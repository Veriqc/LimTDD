# LimTDDexpr 工作历史：修改与 Debug 总结

> 本文档汇总 LimTDDexpr 项目历史上的主要修改与 debug 情况，供后续接手者快速建立全局认知。
> 详细过程见仓库根目录的 `prompts.md`、`prompts_expr.md`、`prompts_expr2.md`、`prompts_prepare.md`、`debug_tdd_c_fidelity.md`、`diagnosis_usage.md` 以及 `LimTDD/docs/*`。

---

## 一、项目定位与总目标

LimTDDexpr 是一个量子电路实验项目，核心是 **LimTDD**（Tensor Decision Diagram，张量决策图，C++20 后端），辅以两个对照后端 **TDD_C** 和 **mqt-limdd**，外加一套 Python 编排层（`expr/`）。

工作沿三条主线推进，时间上大致先后：

1. **节点数 / 运行时调试**（Clifford 与低 T 比例 Clifford+T 电路）——最早的工作面。
2. **Jamiolkowski fidelity 实验**（LimTDD / TDD_C / mqt-limdd 三后端）——工作量最大、debug 最深的一条线。
3. **QReach 后端替换**（用 LimTDD 的 DD 实现替换 CFLOBDD）——最新一条线。

---

## 二、主线一：Clifford / Clifford+T 节点数调试

早期目标不是「泛跑 benchmark」，而是**恢复健康基线、区分「实现 bug」与「表示能力边界」**。

**已恢复的健康基线：**

- Clifford `1time` / `600time`：`11/11` 节点。
- Clifford+T（`cliffordT0.02big20`）`600time` / `713time` / `655time` 显著恢复。

**已明确的关键修复方向：**

- 早期大回归与 map/phase 边界处理错误有关。
- `Slicing` / `Slicing2` 旧语义残留曾导致健康 Clifford 样例爆炸，已修。
- `mapmul` / `mapdiv` 返回 phaseful map 并进入 cache 是低 T 坏例的重要来源，做过结构性修复。
- `normalize` 中 `abs` 与 `std::abs` 的差异是可信的编译器敏感点。

**剩余坏例（可能接近表示能力边界，而非明确 bug）：**

- `490time.qasm`：24588 节点。
- `334time.qasm`：65589 节点。
- 二者的第一跳都出现在 contraction 侧，而非张量构造侧。

---

## 三、主线二：Fidelity 实验（最核心、debug 最深）

### 3.1 理论规约

目标量是 Jamiolkowski fidelity：

$$F_J(U,V)=\frac{|\mathrm{tr}(U^\dagger V)|^2}{2^{2n}}$$

通过最大纠缠态恒等式规约成 2n-qubit 的 overlap：

$$B^\dagger(I\otimes U^\dagger V)B\,|0^{2n}\rangle$$

只需读全零振幅。LimTDD 与 TDD_C 都实现了这条路径，且已避免「单门扩张成全系统大矩阵」，但仍是**顺序逐 gate contraction**（非全局最优 TN 收缩计划）。

### 3.2 LimTDD 后端 —— 修复的 bug

#### Bug #1：ComplexTable 绝对容差碰撞（`LimTDD/DDPackage/dd/ComplexTable.hpp`）

这是最长、最难的一次 debug。

- **症状**：`fidelity > 1`（非物理）。`dj_60` ≈ 2，且满足 $fidelity(n)=2^{2n-81}$（n ≥ 41）。任何 ≥41 qubit 电路在最后一个 evaluation gate 都会触发，与 gate 类型、oracle 结构无关。
- **根因**：`approximatelyEquals` 用绝对容差 `TOLERANCE = epsilon*1024 ≈ 2.27e-13`。n ≥ 41 时收缩权重衰减到 ~$2^{-41}$ 以下，不同的小权重被「容差内判等」（返回 stale 缓存指针），`approximatelyZero` 又把合法小权重清零。
- **最终修复**：`approximatelyEquals` 改为**相对容差**（`TOLERANCE * max(|a|,|b|)`）；`approximatelyZero` 保持绝对容差不变。

**debug 过程中证伪的方向（重要，避免重走）：**

- 扩大 `ifContract(float k)` 判定；
- 禁用 `mapdiv` lookup writeback；
- 尾部手工乘 `1/sqrt(2)`；
- hard guard `the_maps_header()->extra_phase = 0`；
- `Slicing2` 兄弟 rotate 传播（会破坏 `ae_10`）。

debug 主复现器从完整 `dj_60` 换成更快的 `/tmp/dj_pattern_41.qasm`（step 405 / prefix 406），沿 `cont2` → `eq_contract`（`newk1==newk2==39.5`）→ `Slicing2`（`map_level=-1` vs `map_level=0, rot=4`）逐层下探，最终收敛到容差碰撞这一根因。

#### Bug #2：u3 参数顺序错误（`LimTDD/DDPackage/Cir_import.h`）

- **症状**：`grover-noancilla_7`（7q）e=0 时 fidelity=0.14，`qwalk-noancilla_7` e=0 时 7.5e-5。
- **根因**：`U3mat(lambda, phi, theta)` 被以 `(theta, phi, lambda)` 调用，theta 与 lambda 交换。
- **修复**：`U3mat(parameters[2], parameters[1], parameters[0])`，单行修复。
- **最小复现**：单 qubit 单 `u3(pi/4, 5*pi/8, -pi/2)` 门即可触发。

#### Bug #3：sx 门不支持（`LimTDD/DDPackage/Cir_import.h`）

- **症状**：`portfoliovqe_8` 含 82 个 `sx` 门，之前失败。
- **修复**：`supportGate` 增加 `{"sx", SXmat}` 和 `{"sxdg", SXdagmat}`。

#### 更早已修：gate 名查找

`Cir_import.h` 从 `op->getName()` 改为 canonical `qc::toString(op->getType())`，否则逆电路里的 `Sdg` / `Tdg` 会查找失败或错配。

#### 结果

2026-07-05 全量 **18/18 benchmarks 通过**（e=0 → fidelity=1，e=1 → fidelity=0）。

> 注意：e=1 → 0 是数学精确正确的（单 Pauli 无迹），要得到非零 fidelity 需要 error_count ≥ 2。

### 3.3 TDD_C 后端

- **macOS 移植修复**：清除 stale Linux build cache、替换 GNU 专用 `std::exception_ptr::__cxa_exception_type`、改用 workspace 级 vendored xtensor/xtl 头、使用 Release 构建（Debug 会在 `ComplexNumbers::mul` 触发断言）。
- **Bell-pair 节点爆炸**：根因是默认变量序 `get_var_order` 把 Bell 伙伴排得太远，节点数约 $3\cdot 2^n - 2$ 指数爆炸。实验性 `TDD_C_BELL_ORDER=1` 使 Bell-only 线性化（empty_60 只需 181 节点）。
- **fidelity ≈ 2^39 bug**：根因与 LimTDD 相同（ComplexTable 绝对容差），但 TDD_C **没有 edge-map** 托底，所以修复是两步：
  1. `approximatelyEquals` 改相对容差（同 LimTDD）；
  2. `test_fidelity.cpp` 里 per-qubit 初始权重 ×2 缩放（把权重抬到容差之上），最后除以 $2^{eval\_qubits}$ 恢复。
- 验证通过：ae_10、empty_20/60、dj_60、ghz_60、graphstate_30 全部 fidelity=1。

### 3.4 mqt-limdd 探索

只做了机制性探索，**未形成可信 baseline**：`ddsim_simple --pv` 的 statevector 导出在最小 QASM 上就语义可疑（如 `x q[0]` 导出仍像 `|0⟩`），问题在 parser / simulate / state export 的某一层。结论：不要把它纳入论文对照表，除非先单独验证 `ddsim_simple` 的正确性。

---

## 四、主线三：QReach 后端替换（最新）

把 LimTDD 的 DD 包装成 QReach 的 drop-in 后端，替换 CFLOBDD。实现是 header-only，位于 `LimTDD/DDPackage/dd/backend/`。

**已完成：**

- `DDVector`（13 个函数）+ `DDMatrix`（23 个函数）全部实现，58 个单元测试通过。
- 关键紧凑化（从密集 O(4^n) 降到与 n 无关）：门构造（小张量 + `cont`）、状态构造（逐 qubit 张量）、`KroneckerProduct`（`cont` + key 平移）。
- **约定变更（2026-08-17）**：`level` 参数改为真实量子数 `n`（不再 pad 到 2 的幂）。

**修复 / 规避的 bug：**

- **InnerProduct 缩放 bug**：`cont` 全缩并标量路径残留 `v=0` 节点，`sumRemaining` 叠加出 `2^(n-1)` 缩放（且叠加项非均匀，如 `<v1|v3>` 多出 -3 因子），破坏 Gram-Schmidt。**规避**：改用稀疏点积（枚举非零振幅 + big-endian 基索引 hash join），绕开 `cont` 标量路径。

**仍未完成（有分析计划 `limtdd-backend-conjugate-transpose-plan.md`）：**

- `Conjugate` / `Transpose`（矩阵）/ `MatrixMultiply` 仍是密集实现，只对 n ≤ 12 可用。

---

## 五、当前状态一句话总结

- **LimTDD fidelity**：已修复并稳定（`ae_10`、`dj_60` e=0 → 1、e=1 → 0，全 18 benchmark 通过），当前目标是保持健康基线不回退。
- **TDD_C**：已能跑通 fidelity（macOS 已适配），但 trace / debug 工具链比 LimTDD 少。
- **mqt-limdd**：暂无可信 baseline。
- **QReach 后端**：门 / 态 / Kron 已紧凑化，InnerProduct 已绕开标量收缩 bug，矩阵代数紧凑化是下一块硬骨头。

---

## 六、贯穿始终的方法论提醒

多份文档反复强调：

- **节点数下降 ≠ 语义正确**。改核心 DD 逻辑后必须重跑小样例语义对照（`expr/compare_cpp_state_dirac.py`）和 fidelity sanity 用例（`ae_10` / `dj_60`）。
- **不要重走已证伪的修复方向**（见 §3.2 Bug #1 的证伪清单）。
- 实验性环境变量（`LIMTDD_*` / `TDD_C_*`）是诊断探针，不是默认修复手段。

---

## 七、334time 非确定性收缩 bug（2026-08，最新）

`LimTDD/Benchmark/cliffordT0.02big20/334time.qasm`（20 qubit，600 gate）间歇性「耗时过长」，节点数每次运行随机爆炸（prefix=560 三次 419k / 17.5k / 18.8k），并非「固定难」而是**非确定**。

**根因（已定位，未修复）**：map 运算（`mapmul`/`mapdiv`）的 memoization 缓存（`ComputeTable3`，direct-mapped + 指针哈希 + 从不 clear）命中时把 `extra_phase` **写回共享 map 节点**，而 `append_new_map` 缓存键只看 `(level,x,rotate)` 不含 phase → 同一骨架节点被不同 phase 复用 → 别名腐蚀。ASLR 让指针哈希的碰撞/命中模式每次不同 → 结果非确定。

- memoization 是 **2023-10 原始设计**（`temp2` 也有），不是后来引入；**暴露点**是 2026-05-28 的相对容差修复（`7fce6f4`）。
- **TDD_C** 用同样相对容差但无 `the_maps`/map memoization，故无此问题。

**已证伪的方向**（详见 `LimTDD/docs/limtdd-cliffordt-334time-nondeterminism.md`）：确定性/内容哈希、线性探测、周期 clear、禁用 writeback、删除 base 相位清零（打坏 `ae_10`）；禁用 memoization 能收敛但 ~5x 慢（输给 TDD）；**方向 1（把 `extra_phase` 纳入 map 节点身份）**——`ae_10` 过但 `dj_60` 从 121 节点爆炸到 13GB（Clifford 抵消失效）。

**当前结论**：`extra_phase` 是「边上的 pending phase」，不是 map 结构身份；`mapmul`/`mapdiv` 从不读输入 phase 是**正确不变量**，旧代码靠别名（同结构⇒同节点⇒同 phase）让其成立。正确修法是**方向 2：把 phase 从 `the_maps` 结构拆出**（结构键保持 `(level,x,rotate)`，`mapmul`/`mapdiv` 返回 `(结构, phase)`，phase 由边/调用方持有并在 `cos` 处消费）。设计与逐处语义对照见 [`LimTDD/docs/limtdd-extra-phase-immutable-refactor-plan.md`](docs/limtdd-extra-phase-immutable-refactor-plan.md) §2。教训：`dj_60` 必须加入 sanity 回归——`ae_10` 过不代表 Clifford 抵消没坏。

## 八、非确定性根治：两处 ASLR 依赖 → 确定性 ID（2026-08-30）

方向 2 落地并验证正确后，334time 仍有残留非确定（prefix=500 三次 13707/13638/13638）。追查发现真正根因是**两处「按指针地址比较/哈希」的 ASLR 依赖**，均已用「确定性创建序号 ID」替换：

1. **`T_add2` 的 `if (x.p > y.p)`**（DD 加法操作数顺序）→ 指针地址随机 → 浮点求和顺序随机 → 1-ULP 权重差。这是 §8.4(d) 一直没找到的「1-ULP 最终来源」、也是 §8.5 观察「ASan 使其确定」的原因。改为 `x.p->id > y.p->id`（`mNode.id`，`UniqueTable::getNode` 分配序号）。
2. **`hash<Complex>` 的 `reinterpret_cast(指针)`**（桶分布 → 缓存命中/未命中模式）→ 改为表项 ID（`ComplexTable::Entry.id`，`getEntry` 分配序号）。

关键认识：之前的「内容哈希」（`llround(值/tolerance)`）虽然确定，但用**绝对容差**与 `approximatelyEquals` 的**相对容差**不一致 → 大权重时哈希分得过细 → 唯一表去重漏掉 → 节点爆炸。**entry ID 既确定又和去重一致**（近似相等 ⇒ 同一表项 ⇒ 同一 ID），因此同时拿到「确定」与「不爆炸」。

**结果**：334time prefix=450/500/550 全部收敛为单一值（1183 / 5067 / 402089，连跑全同）；Clifford `1time` = 11；`ae_10`/`dj_60` fidelity 正确。450/500 节点数已恢复；550 的 402089 是**确定但膨胀**，属 §8.6 的 ±1/-i 全局相位非规范性（独立于非确定性）。详见 [`LimTDD/docs/limtdd-aslr-nondeterminism-fixed.md`](docs/limtdd-aslr-nondeterminism-fixed.md)。

