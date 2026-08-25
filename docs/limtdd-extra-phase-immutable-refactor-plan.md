# LimTDD `extra_phase` 重构方案：方向 2（拆出 pending phase；方向 1 已证伪）

> **日期**: 2026-08-23
> **状态**: 方向 1（phase 纳入节点身份）**已实施并证伪**；本文档给出正确方向 2 的设计方案（尚未实施）。
> **目标**: 在**保持相对容差**、**保留 map memoization（性能不输 TDD）** 的前提下，根除 `334time.qasm` 的非确定性节点爆炸。
> **关联**: 现象与根因见 [`limtdd-cliffordt-334time-nondeterminism.md`](./limtdd-cliffordt-334time-nondeterminism.md)；`extra_phase` 的历史语义实验见仓库根 `prompts_expr.md §12`。

---

## 0. 结论摘要（更新）

**方向 1（把 `extra_phase` 纳入 `append_new_map` 缓存键、使 map 节点按 phase 区分）是错误的，已实施并证伪**。它使 `ae_10` 正确，但使 `dj_60` 从 121 节点爆炸到 13GB 内存（Clifford 抵消被破坏）。根本原因见 §1。

**正确方案是方向 2：把「pending phase」从 map 结构中彻底拆出**，让 `the_maps` 只表示结构（`level/x/rotate` 链，不可变、按 `(level,x,rotate)` 去重），phase 作为**独立值**由边/调用方持有、在 `cos(phase·π/4)` 处消费。这与「carry outside」协议一致，且使 memoization 真正变纯。

---

## 1. 关键发现：`extra_phase` 不是 map 的结构身份，而是「边上的 pending phase」

这是方向 1 失败的根因，也是理解整个 map/phase 语义的关键。

### 1.1 两个正交的量

- **`rotate` / `x` / `level`（结构）**：`the_maps` 链上的局部置换与旋转相位。`mapmul`/`mapdiv` **只读**这些字段做组合（`Package.hpp:1262-1268 / 1404-1411`）。
- **`extra_phase`（pending phase）**：一个**待施加的全局相位**，属于**边（权重）**，不是 map 结构的一部分。它由 `mapmul`/`mapdiv` 的结构组合「溢出」产生，或由 `normalize` 的权重相位 `rot` 叠加，最终在消费点以 `cos(extra_phase·π/4)` **乘进权重**。

### 1.2 一个正确的不变量：`mapmul`/`mapdiv` 从不读输入 `extra_phase`

`mapmul`/`mapdiv` 组合**结构**，忽略输入的 `extra_phase`。这是**正确**的，因为 pending phase 由调用方在 `cos` 处消费（"carry outside" 协议），而不是 map 结构的一部分。

### 1.3 旧代码靠「别名」让该不变量成立，方向 1 打破了它

旧代码里 `append_new_map` 键只看 `(level,x,rotate)`（不含 phase），于是：

```
同结构 map ⇒ 同一个节点 ⇒ 同一个 extra_phase（共享可变字段）
```

`mapdiv` 的 `self == other` 快捷路径、以及「忽略输入 phase」都**依赖这个隐式不变量**：两个同结构的 map 相位必然相同，所以「忽略输入 phase」等价于「相位差恒为 0」。

**方向 1（phase 纳入键）把「同结构不同 phase」拆成了不同节点**，打破了这条不变量：`mapdiv`/`mapmul` 对两个「同结构不同 phase」的 map 运算时，**丢掉了本应保留的 phase 差**（`self == other` 不命中 → 递归 → 返回 header 而非 phase 差）→ Clifford 抵消失效 → 节点爆炸。

**实测证据**（prefix 二分）：

| prefix | max_nodes |
|---|---|
| 120 | 121（正常，Bell-pair 制备结束） |
| 126 | 3700 |
| 132 | 208942 |

爆炸点恰在 Bell-pair 制备结束后、U^†/V（H 门）开始作用处——即 Clifford 抵消开始失效处。

### 1.4 结论

「把 phase 纳入节点身份」在概念上就是错的：**phase 是 pending edge phase，不是 map 结构身份**。正确做法是把二者**分开**，而不是把它们**合并进同一个身份**。

---

## 2. 方向 2 设计：拆出 pending phase

### 2.1 核心改动

1. **`the_maps` 去掉 `extra_phase` 字段**（或保留但恒为 0），只表示结构。
2. **`append_new_map` 键回到 `(level, x, rotate)`**（不含 phase），结构节点不可变、按结构去重。
3. **`mapmul`/`mapdiv` 返回 `(结构, phase)`**，其中 phase 是「结构组合产生的溢出」，输入 phase 依旧不读（保持 carry-outside 协议）。
4. **phase 由调用方持有并显式合并**：`normalize`/`cont2` 里 `extra_phase + rot` 等操作从「改写共享节点字段」变成「独立 int 运算」。
5. **memoization 变纯**：`ComputeTable3` 缓存 `(self, other) → (result_structure, result_phase)`，命中/未命中都返回同一 `(结构, phase)`，无 writeback、无别名。

### 2.2 phase 的存放位置（已定稿，无需 `Edge.phase`）

进一步深读确认：**phase 总是由 `mapmul`/`mapdiv` 产生、并在同一调用点消费**（折进权重），只有两处「跨行 carry」——`normalize` 的晋升（可用 local 变量）和 `find_remain_map` 的 remain_phase（用 `comm_maps.remain_phase`）。因此 **不需要给 `Edge` 加 `phase` 字段**。

`mapmul`/`mapdiv` 返回 `struct map_res { the_maps* map; int phase; }`，phase 由调用点立即折进权重或短暂 carry。

**完整逐行实施清单见 [`limtdd-extra-phase-direction2-implementation.md`](./limtdd-extra-phase-direction2-implementation.md)**（含 phase 生命周期分类、逐函数改动、实施顺序与回归）。

### 2.3 逐处语义对照（方向 2 视角）

| 位置 | 旧（可变字段） | 方向 2 |
|---|---|---|
| `append_new_map` | 键 `level_x_rotate`，节点建后 `extra_phase=0` | 键 `level_x_rotate`，无 phase 字段 |
| `mapmul`/`mapdiv` 递归 | `res->extra_phase = r->extra_phase + local` | 返回 `phase = r_phase + local`（结构单独返回） |
| base case | `other->extra_phase = 0; return other` | 返回 `(other 结构, phase=0)` |
| `ComputeTable3::lookup` | writeback `entry.result->extra_phase = entry.extra_phase` | 命中返回 `(entry.result, entry.extra_phase)`，无 writeback |
| `normalize` 709/778 | `map->extra_phase += rot` | `phase = mapdiv(...).phase + rot`（存边） |
| `normalize` 849 | `append_new_map(..., child.map->extra_phase)` | `append_new_map(..., child.phase)` |
| `find_remain_map` | `remain_map->extra_phase += ...` | `comm_maps.remain_phase += ...`（已在方向 1 中验证此拆分可行） |
| `cont2` 2856/2870 | `cos(e.map->extra_phase)` / `cos(extra_phase)` | `cos(e.phase)` / `cos(remain_phase)` |
| 各收缩点 1755/1771/…/2342/3641 | `cos(e.map->extra_phase)` | `cos(e.phase)` |

### 2.4 风险与注意

1. **`Edge` 改动影响面大**：`Edge<mNode>` 是 DD 核心结构，加 `phase` 字段会波及相等/哈希/序列化/唯一表。需谨慎。
2. **phase 合并的协议要重新核对**：`cont2` 里 `e.phase`、`remain_phase`、`mapmul 溢出` 三者的合并顺序（谁先谁后、是否重复计数）需逐行确认——这正是方向 1 里没有吃透的部分。
3. **不重走方向 1**：不要再把 phase 塞进 `append_new_map` 键。结构去重键必须保持 `(level,x,rotate)`。

---

## 3. 已证伪：方向 1（记录留档）

方向 1 = 把 `extra_phase` 纳入 `append_new_map` 键 + 移除 writeback + base 查表 + `find_remain_map` carry 拆分 + `normalize` 查表。

结果：

| sanity | 结果 |
|---|---|
| `ae_10` e=0 / e=1 | ✅ fidelity 1 / 0 |
| `dj_60` e=0 | ❌ 121 节点 → 13GB 内存 / >6min 未完成（Clifford 抵消失效） |

**根因**：见 §1.3。phase 是 pending edge phase，方向 1 把它当成了结构身份，打破了「同结构 ⇒ 同 phase」的隐式不变量，使 `mapdiv`/`mapmul` 对「同结构不同 phase」丢掉相位差。

所有改动已 `git checkout` 回滚，树干净。

---

## 4. 验证协议（方向 2 实施时复用）

```bash
cmake -S LimTDD -B LimTDD/build && cmake --build LimTDD/build -j2
# 正确性 sanity（任何核心改动后必跑）
./LimTDD/build/test/test_fidelity expr/qasm/ae_10.qasm 0 1   # 期望 fidelity 1
./LimTDD/build/test/test_fidelity expr/qasm/ae_10.qasm 1 1   # 期望 fidelity 0
./LimTDD/build/test/test_fidelity expr/qasm/dj_60.qasm 0 1    # 期望 fidelity 1，且快速完成（max_nodes≈121）
# 非确定性验证（连跑 ≥5 次，MAX node 应收敛）
LIMTDD_TN_PREFIX=550 ./LimTDD/build/test/test_data LimTDD/Benchmark/cliffordT0.02big20/334time.qasm 00000000000000000000
```

> ⚠️ `dj_60` 必须显式加入 sanity 回归：方向 1 的教训是「`ae_10` 过不代表 Clifford 抵消没坏」，而 334time 的非确定性正是 Clifford+T 电路的相位问题。
