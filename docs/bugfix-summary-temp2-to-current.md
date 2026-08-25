# LimTDDexpr 实验:temp2 → 当前分支 Bug Fix 汇总

> **用途**:展示 / 汇报用。每个 bug 按 **表象 → 发生点 → 原因 → 解决** 四段给出。
> **范围**:把实验从旧 `temp2` 基线搬到当前分支、并让 fidelity / Clifford+T 两条主线跑通所做的全部修复,以及当前仍待处理的 open issue。
> **详细过程**见 `HISTORY.md`、`LimTDD/docs/*`、`debug_tdd_c_fidelity.md`。

---

## 总览(一页速览)

| # | Bug | 后端 | 状态 |
|---|---|---|---|
| 1 | 绝对容差碰撞 → `fidelity > 1` | LimTDD | ✅ 已修复 |
| 2 | `u3` 参数顺序(theta/lambda 交换) | LimTDD | ✅ 已修复 |
| 3 | `sx` 门不支持 | LimTDD | ✅ 已修复 |
| 4 | 逆电路门名 `Sdg`/`Tdg` 查找失败 | LimTDD | ✅ 已修复 |
| 5 | macOS 移植编译失败 | TDD_C | ✅ 已修复 |
| 6 | Bell-pair 节点指数爆炸 | TDD_C | ✅ 已修复 |
| 7 | `fidelity ≈ 2^39` | TDD_C | ✅ 已修复 |
| 8 | 334time 非确定性节点爆炸 | LimTDD | 🔴 根因定位,方向 2 待实施 |

> ⚠️ **贯穿线索**:修复 Bug #1 的相对容差改动(`7fce6f4`)**恰好暴露**了 Bug #8 这个潜伏已久的 map memoization 非确定性。

---

## 一、LimTDD 后端 —— fidelity 相关

### Bug #1:ComplexTable 绝对容差碰撞(最长的 debug)

| 维度 | 内容 |
|---|---|
| **表象** | `fidelity > 1`(非物理)。`dj_60` ≈ 2;满足 `fidelity(n) = 2^(2n-81)`(n ≥ 41)。任意 ≥41 qubit 电路在最后一个 evaluation gate 触发,与门类型、oracle 结构无关 |
| **发生点** | `LimTDD/DDPackage/dd/ComplexTable.hpp` 的 `approximatelyEquals` |
| **原因** | 用**绝对**容差 `TOLERANCE = ε·1024 ≈ 2.27e-13`。n ≥ 41 时收缩权重衰减到 ~`2^(-41)` 以下,不同的小权重被「容差内判等」→ 返回 stale 缓存指针;`approximatelyZero` 又把合法小权重清零 |
| **解决** | `approximatelyEquals` 改**相对**容差 `TOLERANCE · max(|a|,|b|)`;`approximatelyZero` 保持绝对容差不变 |

**证伪过的方向(勿重走)**:扩大 `ifContract`、禁用 mapdiv writeback、尾部手工乘 `1/√2`、hard-guard `extra_phase=0`、`Slicing2` 兄弟 rotate 传播(破坏 `ae_10`)。

### Bug #2:`u3` 参数顺序错误

| 维度 | 内容 |
|---|---|
| **表象** | `grover-noancilla_7` e=0 → fidelity 0.14;`qwalk-noancilla_7` e=0 → 7.5e-5 |
| **发生点** | `LimTDD/DDPackage/Cir_import.h` |
| **原因** | `U3mat(lambda, phi, theta)` 被以 `(theta, phi, lambda)` 调用,theta 与 lambda 互换 |
| **解决** | 单行:`U3mat(parameters[2], parameters[1], parameters[0])`。最小复现 = 单 qubit 单 `u3(π/4, 5π/8, -π/2)` 门 |

### Bug #3:`sx` 门不支持

| 维度 | 内容 |
|---|---|
| **表象** | `portfoliovqe_8`(含 82 个 `sx` 门)运行失败 |
| **发生点** | `LimTDD/DDPackage/Cir_import.h` 的 `supportGate` |
| **原因** | 门分发表缺少 `sx` / `sxdg` |
| **解决** | 增加 `{"sx", SXmat}`、`{"sxdg", SXdagmat}` |

### Bug #4:逆电路门名查找失败

| 维度 | 内容 |
|---|---|
| **表象** | 逆电路(如 `U^dagger` / `B^dagger` 里)的门解析错误 |
| **发生点** | `LimTDD/DDPackage/Cir_import.h` 门分派 |
| **原因** | 用 `op->getName()`(可变显示名),逆电路的 `Sdg`/`Tdg` 查找失败或错配 |
| **解决** | 改用 canonical 类型名 `qc::toString(op->getType())` |

**结果**:2026-07-05 全量 **18/18 benchmarks 通过**(e=0 → fidelity 1,e=1 → fidelity 0)。

---

## 二、TDD_C 后端(对照)

### Bug #5:macOS 移植编译失败

| 维度 | 内容 |
|---|---|
| **表象** | 移植到 macOS 无法编译 / 运行 |
| **发生点** | build 缓存 + `ComplexNumbers` 等 |
| **原因** | stale Linux build 缓存;GNU 专用 `std::exception_ptr::__cxa_exception_type`;依赖系统级 xtensor/xtl;Debug 构建在 `ComplexNumbers::mul` 触发断言 |
| **解决** | 清缓存、去掉 GNU 专用类型、改用 workspace 级 vendored `include/xtensor*`、改 Release 构建 |

### Bug #6:Bell-pair 节点指数爆炸

| 维度 | 内容 |
|---|---|
| **表象** | 制备 Bell 对时节点数指数爆炸(约 `3·2^n - 2`) |
| **发生点** | 默认变量序 `get_var_order` |
| **原因** | 把 Bell 伙伴排得太远,DD 无法共享 |
| **解决** | 实验性 `TDD_C_BELL_ORDER=1` 使 Bell-only 线性化(`empty_60` 降到 181 节点) |

### Bug #7:`fidelity ≈ 2^39`

| 维度 | 内容 |
|---|---|
| **表象** | fidelity 量级 ~2^39,非物理 |
| **发生点** | `ComplexTable.hpp`(与 LimTDD Bug #1 同根因) |
| **原因** | 同样绝对容差碰撞;但 TDD_C **无 edge-map 托底**,权重无下限保护 |
| **解决** | 两步:①`approximatelyEquals` 改相对容差;②`test_fidelity.cpp` 里 per-qubit 初始权重 ×2 缩放(抬到容差之上),最后除 `2^(eval_qubits)` 恢复 |

**结果**:`ae_10`、`empty_20/60`、`dj_60`、`ghz_60`、`graphstate_30` 全部 fidelity = 1。

---

## 三、mqt-limdd(结论,非 bug)

`ddsim_simple --pv` 的 statevector 导出在最小 QASM 上即语义可疑(如 `x q[0]` 导出仍像 `|0⟩`),问题在 parser / simulate / export 某一层。**未形成可信 baseline,不建议纳入对照表**。

---

## 四、Clifford+T:334time 非确定性(🔴 进行中)

| 维度 | 内容 |
|---|---|
| **表象** | `334time.qasm`(20q,600 gate)间歇性耗时过长,节点数每次随机爆炸:prefix=560 三次 419k / 17.5k / 18.8k |
| **发生点** | `Package.hpp` 的 `mapmul`/`mapdiv` + `ComputeTable.hpp` 的 `ComputeTable3` |
| **原因** | map memoization 缓存 direct-mapped + 指针哈希 + 从不 clear;命中时把 `extra_phase` **写回共享 map 节点**,而 `append_new_map` 键只看 `(level,x,rotate)` → 别名腐蚀;ASLR 使碰撞/命中模式每次不同 → 结果非确定 |
| **状态** | 根因已定位。**方向 1**(phase 纳入节点身份)已实施并**证伪**(`dj_60` 121→13GB)。**方向 2**(把 phase 从 `the_maps` 结构拆出、`mapmul`/`mapdiv` 返回 `(结构, phase)`)已设计,尚未实施 |

**已证伪方向**:确定性/内容哈希、线性探测、周期 clear、禁用 writeback、删 base 相位清零、方向 1。
**决定性证据**:禁用 map memoization → 节点收敛 ±0.7%(但 ~5x 慢,输给 TDD)。

**为什么 LimTDD 在 Clifford+T bad case 上比TDD更弱**

On the hardest Clifford+T instances LimTDD is somewhat slower than the baseline TDD, not merely tied. This is a deliberate design trade-off rather than a bug. LimTDD encodes phases in a linked map structure keyed by (level, x, rotate) and memoizes map composition; this yields the ~10× advantage whenever phase structures are shared across edges. But its per-operation cost is strictly higher than the baseline's inline complex-phase multiplication — each composition walks a parent chain, constructs a string key, and performs an ordered-tree lookup — and this cost is amortized only while memoization hits. On high-T-density circuits where contraction produces largely distinct phase structures, sharing collapses and the amortization vanishes, exposing the higher constant factor. A second, pathological effect compounds this: the direct-mapped, never-flushed cache saturates, and its write-back of the pending phase into structurally shared map nodes aliases phases that should differ, further inflating the node count. The baseline, which stores phases directly in edge weights and keeps no map cache, is immune to both effects.

---

## 五、当前现状总结

**Fidelity(已稳)**
- **LimTDD**:18/18 benchmark 通过,`ae_10`/`dj_60` e=0 → 1、e=1 → 0,当前目标是保持健康基线不回退。
- **TDD_C**:macOS 已跑通,`ae_10`/`empty_20/60`/`dj_60`/`ghz_60`/`graphstate_30` 全 fidelity = 1,但 trace/debug 工具链比 LimTDD 少。
- **mqt-limdd**:暂无可信 baseline。

**Clifford+T(还差最后一刀)**
- **健康基线已恢复**:Clifford `1time`/`600time` = 11/11 节点;`cliffordT0.02big20` 的 `600time`/`713time`/`655time` 显著恢复。
- **剩余坏例**:`334time`(65589 节点)+ `490time`(24588 节点),第一跳都在 contraction 侧。
- **334time 的非确定性已根因定位**,正确修法是 **方向 2(拆出 pending phase,使 map memoization 变纯)**,设计与逐处对照已备好,待实施;关键教训是 `dj_60` 必须加入 sanity 回归——`ae_10` 过不代表 Clifford 抵消没坏。
