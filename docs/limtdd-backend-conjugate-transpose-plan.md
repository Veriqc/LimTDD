# LimTDD 后端 —— Conjugate / Transpose / MatrixMultiply 紧凑化分析与实现计划

> **Audience:** 下一次接手本后端的 agent(可以是 Claude 或人)。
> **Date:** 2026-08-17
> **状态:** 分析已完成,**尚未实现**。这是把「矩阵代数」从密集 O(4^n) 降到次指数的最后一块,但比门构造/状态构造/Kron 那几项更难,需要单独一轮带交叉验证的实现。

---

## 0. 背景:什么已紧凑、什么还密集

截至 2026-08-17 已完成的紧凑化(均验证通过,见 `limtdd-backend-implementation-status.md`):

- **门构造**:`MkSingleQubitGateOnN`/`MkCNOT`/`MkCCNOT`/`MkSwap`/`MkiSwap`/`MkCP` = 小张量(arity ≤ 3)+ `cont`。
- **状态构造**:`MkBasisVector`/`NoDistinctionNode` = 逐 qubit 张量 + `cont` 张量积(O(n),`kMaxCompactQubits=256`)。
- **`KroneckerProduct`** = `cont(a, shiftKeys(b, na))`(key 平移助手)。

**仍为密集实现(本文目标)**:`MatrixMultiply`(O(8^n))、`Conjugate`(O(2^n) 向量 / O(4^n) 矩阵)、`Transpose`(O(4^n) 矩阵)。三者都只对 n≤12 可用。

---

## 1. 依赖的核心原语(已确认 public 可用)

| 原语 | 调用方式 | 说明 |
|---|---|---|
| 切片 | `limtdd::backendPackage().backendSlice(e, v, c)` | public wrapper(见 `Package.hpp`),带零边守卫 `if (e.w == zero) return e;`,内部调私有 `Slicing`。`Slicing` 会把父边权重/相位乘进子边 |
| 建节点 | `limtdd::backendPackage().makeDDNode(v, edges, cached)` | public,`normalize` + unique-table 插入 |
| 权重 | `limtdd::backendPackage().cn.lookup(re, im)` + `cn.incRef`/`decRef` | `lookup` 返回 **refcount 0** 的 ComplexTable 项,要 `incRef` 才归 DD 拥有 |
| 收缩 | `limtdd::backendPackage().cont(tdd1, tdd2)` | 部分缩并(留输出索引)正确;**全缩并标量路径有 bug**(见 `limtdd-innerproduct-bug-report.md`) |

`backendSlice`/`backendAdd`/`backendSlice2` 是当初为适配层加的 public 转发(对应私有 `Slicing`/`T_add2`/`Slicing2`)。

---

## 2. Conjugate —— 可行路径(递归走树)

**核心思路**:`Slicing` 会把父边权重 `e.w` 和父 map `e.map` 乘进子边(`temp.w = mulCached(child.w, e.w)`、`temp.map = mapmul(e.map, child.map)`)。所以**父边的权重/相位在切片时就"消费"了**——递归共轭子边会自动共轭父边的相位,**不需要手动造共轭 map**。

```
Edge conjEdge(Edge e):
    if e.p == nullptr: return e                 # 零边
    if e.p->v == -1:                            # 终端
        r = e
        r.w = conjWeight(e.w)                   # 只共轭权重;map 默认 header(无相位)
        return r
    c0 = conjEdge(backendSlice(e, e.p->v, 0))
    c1 = conjEdge(backendSlice(e, e.p->v, 1))
    return makeDDNode(e.p->v, {c0, c1}, true)   # normalize 会重建 map
```

这样走的是 `renormalize` 的同款模式(见 `Package.hpp` 的 `renormalize`,它也是 Slicing → 递归 → makeDDNode)。

### 必须处理的坑

1. **权重共轭 + refcount**:`conjWeight(w)` = 用 `CTEntry::val(w.r/i)` 取 double,再 `cn.lookup(re, -im)`。返回 refcount 0,要 `incRef` 才归新 DD 拥有——和 `DDTypes.hpp` 里 `operator*`(scalar·DD)同款坑,之前踩过「ref==0 before decref」。

2. **`the_maps` 的 `extra_phase` 是可变字段**(`Maps.hpp` + `Package.hpp` 的 `append_new_map`):`append_new_map` 缓存只按 `(level, x, rotate)` 去重、**不看 `extra_phase`**,且 `self->next` 是全局共享缓存。**不要用 `append_new_map` 手动造共轭 map**(会污染共享节点)。走上面的 Slicing 路径就绕开了。

3. **终端权重可能是非 header 的 map?** 实际中 `Edge::terminal(w)` 的 map 是默认 header。但为稳妥,终端分支要保留 map 字段(即 `r = e` 再改 `r.w`),不要假设能丢掉 map。

### 验证策略

- 用 dense 版(`matrixToDense`/`GetNonZeroAmplitudes` + 重建)当 oracle,对随机小 n(2~6)状态/矩阵逐用例交叉验证。
- 已有回归用例:`Conjugate((|0>+i|1>)/√2)`、`Conjugate(S)`、`Conjugate`/`Transpose` 的 dense 断言在 `test_ddmatrix.cpp`。
- 相位门(`MkPhaseShift`/`MkU3`/`MkArbitrary`)要特别测——它们的相位走精确 double 权重(不进 map),共轭要正确取反虚部。

---

## 3. Transpose —— 比 Conjugate 硬(节点重排)

- 向量 Transpose 是恒等(已实现)。
- 矩阵 Transpose = 交换行/列 = 交换 `"o{q}"↔"q{q}"` 的**语义角色**,不是简单换 key label。
- **为什么换 label 不行**:`varOrder["o{q}"]=2q`、`varOrder["q{q}"]=2q+1`(交错)。单纯把 `key_2_index`/`index_set` 里的 o↔q,会让「节点变量按 varOrder 有序」不变式破裂(树结构没变但语义反了),`cont` 会错。
- **需要重排节点树**:o/q 变量在树里的父子位置互换。没有现成原语,是最硬的部分。
- 建议:若语义层不常在**矩阵**上做大 n Transpose,押后;先做 Conjugate(向量+矩阵都更常用,resetall 依赖向量共轭)。

---

## 4. MatrixMultiply(顺带,同批剩余 dense 项)

- 原理上是单次索引收缩:C[i][k] = Σ_j A[i][j] B[j][k],即收缩 A 的 col `"q{j}"` 与 B 的 row `"o{j}"`。
- 卡点:当前只有两族 key(`"o"`/`"q"`),A 的 col 和 B 的 row 撞 key,需要**第三族「中间索引」(如 `"m{q}"`)**或通用 key 重命名原语。
- `detail::shiftQubitKeys`(本次为 Kron 加的)是通用重命名的一块积木,但不支持 o↔q 交换 + 中间族。
- 同样建议押后:语义层若只在 gate·vector 上跑(顺序施加门),就不需要矩阵-矩阵乘法。

---

## 5. 建议推进顺序(风险从低到高)

1. **Conjugate 向量**(递归走树,拿 dense 交叉验证)——收益最大,resetall 依赖。
2. **Conjugate 矩阵**(同构,只是变量数是 2n)。
3. **Transpose 矩阵**(节点重排,最硬)。
4. **MatrixMultiply**(中间索引族)。

每一步都单独验证后再进下一步。
