# LimTDD 后端实现现状(供 QReach agent 对接参考)

> **Audience:** QReach 侧的 agent,用于感知 LimTDD 后端目前的实现进度与对接约定。
> **Status:** 实现中 —— 全部 API 已实现、58 个单元测试通过;门构造(小张量 + `cont`)、状态构造(`MkBasisVector`/`NoDistinctionNode`)、`KroneckerProduct` 均已紧凑化;仍为密集实现的只剩 `MatrixMultiply`/`Conjugate`/`Transpose`。
> **Last updated:** 2026-08-17
> **对应文件:** `LimTDD/DDPackage/dd/backend/{DDTypes,DDVector,DDMatrix}.hpp`

## 0. 总体状态

`DDVector` 与 `DDMatrix` 两个命名空间的**全部 API 已实现**,配套 58 个单元测试(35 DDVector + 23 DDMatrix)全部通过。实现是 header-only,位于 `LimTDD/DDPackage/dd/backend/` 下:

- `DDTypes.hpp` — 核心类型 `DD`/`DDComplex` + 包单例 + `Initialize`
- `DDVector.hpp` — 状态向量
- `DDMatrix.hpp` — 门矩阵 + 矩阵代数

测试:`LimTDD/test/test_ddvector.cpp`、`test_ddmatrix.cpp`(已注册进 `test/CMakeLists.txt`)。

## 1. 已实现的 API 全貌

**`DDVector`(契约 §2,13 个函数全部实现)**:
`Initialize`、`MkBasisVector(level,index)`、`MkBasisVector(level,bitstring)`、`NoDistinctionNode`、`InitializeWithAmplitudes(qnum,amps)`、`VectorToMatrixInterleaved`(no-op)、`GetLevel`、`IsApproximatelyZero`、`ExtractSingleAmplitude`、`GetNonZeroAmplitudes`、`Normalize`、`InnerProduct`、`VectorPrintColumnHead`。

**`DDMatrix`(契约 §3,23 个函数全部实现)**:
`Initialize`、9 个单 qubit 门(`MkIdRelation/MkWalsh/MkNegation/MkPauliY/MkPauliZ/MkSGate/MkPhaseShift/MkU3/MkArbitrary`)、5 个多 qubit 门(`MkCNOT/MkCCNOT/MkSwap/MkiSwap/MkCP`)、`MkSingleQubitGateOnN` + `WithParam`/`WithParamVec`、`KroneckerProduct`、`MatrixMultiply`、`MatrixMultiplyWithVector`(热路径)、`Conjugate`、`Transpose`。

## 2. 与设计文档不一致 / QReach 侧必须知道的点

1. **类型不是裸 LimTDD 类型,而是包装类型**:
   - `DD = limtdd::DD`,内部包裹 `dd::TDD`(node 边 + index_set + key_2_index)。
   - `DDComplex = limtdd::DDComplex`,内部包裹 `std::complex<double>`,完整满足契约 §6.3 运算符集(含与整数 0/1 比较、`(re,im)`/`double` 构造)。
   - 命名空间是 **`limtdd::`**,不是 `qreach::`。QReach 侧在 `dd_backend.hpp` 里 `using DD = limtdd::DD; using DDComplex = limtdd::DDComplex;` 即可。

2. **`DD` 是 RAII 引用计数,不是值语义深拷贝**:LimTDD 内部用显式 `incRef`/`decRef`/`garbageCollect` 管理 DAG 生命周期。`DD` 包装成 RAII——拷贝=incRef、析构=decRef、移动=转移所有权。语义层**自由拷贝 `DD` 是安全的**。

3. **✅ 门构造器已改为"紧凑小张量 + `cont`"(不再密集 O(4^n))**。`MkSingleQubitGateOnN`(及 param 变体)、`MkCNOT`、`MkCCNOT`、`MkSwap`、`MkiSwap`、`MkCP` 现在都构建一个只作用于被触及 qubit 的小张量(2×2 / 4×4 / 8×8,arity ≤ 3),由 `MatrixMultiplyWithVector = cont(gate, vec) + rename` 直接收缩到整个态上(即 `test_data.cpp` 的 tensor-network 仿真方式)。构造代价 O(4^arity),**与系统大小 n 无关**;实测 X on 30q 门仅 3 个节点。已加 16-qubit 大 n 回归测试。
   - **✅ 状态构造也紧凑化了**:`MkBasisVector`/`NoDistinctionNode` 现改为「逐 qubit 张量 + `cont` 张量积」,O(n) 而非 O(2^n);实测 32-qubit basis/全 1 态仅 ~33 节点(旧 dense 会 OOM)。`kMaxCompactQubits=256` 单独放宽了上限。
   - **✅ `KroneckerProduct` 已改为 `cont`**:disjoint 张量积 = `cont(a, shiftKeys(b, na))`(一个「key 平移」助手,不改节点变量号)。已实测正确。
   - **仍为密集实现的部分**:矩阵代数 `MatrixMultiply`、`Conjugate`(矩阵)、`Transpose`(矩阵)仍是"枚举/密集 → `to_tdd`",只对 n≤12 可用。这是下一步(DD 级共轭/转置、矩阵乘法)的优化点。
   - 原因见 §4。

4. **`Transpose` 对向量返回自身(恒等)**:契约 §6.8 里 `resetall`(Conjugate∘Transpose 作用于向量)的精确语义**尚未验证**,需要 dense 交叉校验确认。

5. **`VectorToMatrixInterleaved` 是 no-op**(契约允许,因为后端直接对向量做 `cont`)。

6. **精度是 `double`**(不是 CFLOBDD 的 100 位),对应之前定下的 route a。

7. **⚠️ `level` 参数已改为 `n`(真实量子数)**,2026-08-17 QReach 对接改动:QReach 不再把 qubit 数 padding 到 2 的幂,直接把 n 传进原 `level` 参数位。因此:
   - `MkBasisVector(n,…)`、`NoDistinctionNode(n,…)`:首参即 n(不再是 `1<<level`)。
   - `MkSwap(n,i,j)`、`MkiSwap(n,i,j)`、`MkCP(n,ctrl,tgt,θ)`:首参即 n(不再是 `1<<(level-1)`)。
   - `GetLevel(c)` 现在返回真实 n(= `index_set.size()`),不再 `ceil(log2(n))`。
   - `MkCNOT`/`MkCCNOT`/`MkSingleQubitGateOnN`/`InitializeWithAmplitudes` 本来就接受 n,无需改。
   - 支持任意 n(含非 2 幂),n=0 返回标量(已处理 `basisState`/`allOnesState` 的 n=0 分支)。

## 3. 我补充的设计决策(文档里没有、对接需要知道)

1. **Index 命名约定(内部)**:状态向量 qubit q 用 key `"q{q}"`(idx 0);门矩阵 **行(输出)** 用 `"o{q}"`、**列(输入)** 用 `"q{q}"`(idx 0)。这样 `cont(gate, vec)` 收缩列的 `"q{q}"`、留下行的 `"o{q}"`。

2. **`MatrixMultiplyWithVector` = `cont(gate, vec)` + 把结果行 key `"o{q}"` 改回 `"q{q}"`**。改名只改 index_set/key_2_index 元数据,不改节点变量号(安全,因为行/列 key 排序一致)。

3. **⚠️ `Initialize()` 预注册"交错"变量序 `"o0","q0","o1","q1",…`**。这是**正确性必需**:分组序(`o0..oN,q0..qN`)会让 `cont` 产生 factor-2 错误;交错序才正确(对应契约 §6.2 提到的 VOC12 交错序)。预注册保证了无论先建向量还是先建门,变量序都确定。

4. **大端 basis 约定**(契约 §6.2):qubit 0 = MSB;稠密数组的 row-major 线性下标 == 大端基态整数。

5. **角度以 π 为单位**(契约 §6.4):用 `cospi/sinpi`。

6. **`root_of_unit = 8` 相位 map 限制**:map 只能精确表示 π/4 整数倍的相位。任意角门(PhaseShift/U3/Arbitrary)的相位走精确 double 权重(不进 map),正确性无影响,但这些门享受不到 map 压缩。

7. **改了核心文件 `Package.hpp`**:加了 3 个 public 转发包装 `backendAdd`(=私有 `T_add2`)、`backendSlice`(=私有 `Slicing`)、`backendSlice2`(=私有 `Slicing2`),供适配层做加法和振幅切片。纯改可见性,无逻辑改动。

## 4. 已知限制与根因(重要)

1. **门构造/状态构造/Kron 已紧凑化(2026-08-17)**,但矩阵乘法的密集根因仍在:
   - **变量序必须交错**:分组序会 factor-2(已用交错预注册解决)。
   - **`cont` 对 disjoint 张量做 Kron 后再次 `cont` 的"畸形 DD(v=0→v=0)"实为分组序时代的旧症状**——交错预注册后已消失,`KroneckerProduct = cont(a, shiftKeys(b, na))` 实测正确(含三步 cont 再收缩)。真正仍未绕开的是 `MatrixMultiply`(矩阵-矩阵乘法需要第三族"中间索引"或通用 key 重命名原语),它仍走密集构造。这是 `cont2` 的 key-mapping 深坑(与 fidelity 那轮 bug 同源)。

2. **`GetNonZeroAmplitudes`/`Normalize` 枚举时有有界 cache 泄漏**(切片产生的临时权重未 `returnToCache`),不 crash,ComplexCache 会回收。

3. **`Initialize` 预注册上限 `kMaxKey=256`**:≥256 qubit 的 key 会退回自然插入序,可能破坏交错不变式。

4. **`InnerProduct` / `dot`(Gram-Schmidt 热路径)性能**:`InnerProduct(a,b)` 现为**稀疏点积**(两个向量枚举非零振幅 + 按 big-endian 基索引 hash join),O(#非零振幅)。之前是 `cont(conj(a), b)` + 密集共轭,但 `cont` 的**全缩并(标量)路径有数据依赖的缩放 bug**(见 `limtdd-innerproduct-bug-report.md`):残留一个 `v=0` 节点、根权重被 `2^(n-1)` cache 缩放 + map 归一化再乘上一个数据相关因子(如 `<v1|v3>` 的 `-3`)。已改为稀疏点积绕开 `cont` 标量路径。**DD 级共轭/收缩**(修复 `cont` 标量路径,或逐边权重共轭 + map `rotate`/`extra_phase` 取反)仍是下一步把大 n 降到次指数的优化方向,`the_maps` 结构已摸清(`x` 实数不变,`rotate`/`extra_phase` 取反即可)。

## 5. 对接要点(给 QReach agent 的 actionable 清单)

1. `dd_backend.hpp`:`namespace qreach { using DD = limtdd::DD; using DDComplex = limtdd::DDComplex; }`,并 include `dd/backend/DDVector.hpp`、`dd/backend/DDMatrix.hpp`。
2. `DDVector::Initialize()` 和 `DDMatrix::Initialize()` 均幂等,内部都调 `limtdd::Initialize()`。
3. `DD` 默认构造 = 零向量/零矩阵;`DD + DD`、`DDComplex * DD`、`==` 都已定义(§1 要求)。
4. `MkSingleQubitGateOnN` 的 `gate1q` 函数指针会被以 `gate1q(1)` 调用,返回 2×2 再提升。
5. **门·态热路径(`Mk*` + `MatrixMultiplyWithVector`)、状态构造、`KroneckerProduct` 均已可上大 n**(小张量/逐 qubit 张量 + cont,已测 16/30/32 qubit);但**矩阵代数 `MatrixMultiply`/`Conjugate`/`Transpose` 仍是密集 O(4^n)**,若语义层在这些上做大 n 会失效。DD 级共轭/转置的分析与实现计划见 `limtdd-backend-conjugate-transpose-plan.md`。
