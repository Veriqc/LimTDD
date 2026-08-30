# LimTDD Clifford+T 334time 非确定性收缩 bug 报告

> **日期**: 2026-08-22（更新 2026-08-26）
> **状态**: 方向 2（拆出 pending phase）**已实施并验证正确**；底层指针 hash **已全部改为内容哈希**（Clifford 已确定）。**根因已精确定位**：334time 的非确定是 **argmax 在模长相等处选了不同子节点 → 全局相位（±1/-i）歧义 → DD 失去 canonical**。checkpoint 工具把它定位到 step 294。修复方向是**全局相位规范化**（尚未实施）。详见 §8.6。
> **关键结论**: 非确定性由 map 运算的 memoization 缓存（`ComputeTable3`，direct-mapped + 指针哈希 + 从不 clear）导致；该缓存是 **2023-10 就有的原始设计**（`temp2` 分支也有），**不是后来引入的**。真正"暴露"它的是 2026-05-28 的相对容差修复（`7fce6f4 repair abs bug`，fidelity 那轮）。TDD_C 用同样的相对容差但**没有** map memoization（无 `the_maps` 结构），因此无此问题。
> **相关 badcase**: `LimTDD/Benchmark/cliffordT0.02big20/334time.qasm`（20 qubit，600 tensor/gate）
> **相关历史**: `prompts.md` 记录该例为 CliffT0.02 坏例（65589/65501 节点）；本报告是"耗时过长"这一症状的根因追踪。

---

## 1. 现象

`334time.qasm` 运行时**间歇性耗时过长**，且节点数随机爆炸。同一个 `prefix=560` 连跑三次：

| 运行 | 耗时 | MAX node |
|---|---|---|
| A | 111.7s | 419243 |
| B | 16.6s | 17544 |
| C | 20.6s | 18808 |

全量 600 gate 更会跑到 5 分钟以上 + 8GB 内存。**这不是"固定难"的电路，而是每次运行会随机撞进节点爆炸路径的非确定性。**

复现命令：

```bash
./LimTDD/build/test/test_data LimTDD/Benchmark/cliffordT0.02big20/334time.qasm 00000000000000000000
LIMTDD_TN_PREFIX=560 ./LimTDD/build/test/test_data LimTDD/Benchmark/cliffordT0.02big20/334time.qasm 00000000000000000000
```

---

## 2. 非确定性 onset 定位

用 `LIMTDD_TN_PREFIX` 前缀扫描，观察 MAX node 的确定性：

| prefix | 3 次 MAX node | 确定性 |
|---|---|---|
| 300 | 102 / 102 / 102 | ✅ 确定 |
| 400 | 1035 / 1035 / 1035 | ✅ 确定 |
| 450 | 1171 / 1137 / 1171 | ⚠️ 开始漂移 |
| 500 | 1135 / 1554 / 1900 | ❌ 发散 |
| 550 | 17489 / 18090 / 28303 | ❌ 发散 |

非确定性在 **gate 400–450（T-gate 区）** 起爆，之后是随机游走式放大。

---

## 3. 差分定位：第一个发散点

用 `LIMTDD_REGRESSION_DIAG=1` + 逐步 trace（`step_cont[i]` 计数器指纹）对两个 run 逐步 diff（prefix=450，窗口 380–449）：

```
=== FIRST DIVERGENCE at step 381 ===
  nodes: run1 224->224   run2 224->224      ← 节点结构完全一致！
  mapdiv.lookup_phase_overwrite            run1=11  run2=10
  mapdiv.lookup_phase_overwrite_non_header run1=11  run2=10
  mapdiv.result_phaseful                   run1=29  run2=30
  mapmul.lookup_phaseful                   run1=7   run2=6
  mapmul.result_phaseful                   run1=45  run2=40
```

**结论**：首个发散在 step 381；此刻**节点结构完全相同**（都是 224->224），但 **mapdiv/mapmul 的 phase 相关计数器已经不同**。非确定性最早进入的是 map 的 phase 处理，而非节点结构。

---

## 4. 根因

`mapmulTable` / `mapdivTable` 的 memoization 表（`ComputeTable3`）是 **direct-mapped 哈希表**，且键是 **map 指针地址**。

### 4.1 键是 ASLR 随机的指针地址

`ComputeTable.hpp:235-240`：

```cpp
static std::size_t hash(const LeftOperandType& leftOperand, const RightOperandType& rightOperand) {
    const auto h1 = std::hash<LeftOperandType>{}(leftOperand);   // the_maps* → 地址
    const auto h2 = std::hash<RightOperandType>{}(rightOperand);
    return combineHash(h1, h2) & MASK;
}
```

`std::hash<the_maps*>` 是默认指针哈希（地址），每次进程 ASLR 不同。

### 4.2 冲突静默覆盖（direct-mapped）

`ComputeTable.hpp:262-267`（`insert`）：

```cpp
table[key] = { leftOperand, rightOperand, result, c};   // 单槽直接覆盖，无链表/探测
```

16384 个槽、每槽 1 个 entry。两个不同 `(map1, map2)` 撞到同一槽，前一个被静默挤掉。`findEntry` 发现 `leftOperand != ...` 就返回 miss。

### 4.3 关键：memoization 有副作用（不纯）

`ComputeTable.hpp:288`（`lookup` 命中时）与 `Package.hpp:1366`（`mapdiv` 命中时）：

```cpp
entry.result->extra_phase = entry.extra_phase;   // writeback，改写共享 map 的 extra_phase
```

而 miss 时走重算路径 `res->extra_phase = r->extra_phase`（`Package.hpp:1375` 等）。**hit 与 miss 产生不同的 phase**。

### 4.4 完整链条

```
ASLR → 指针地址 → direct-mapped 碰撞/挤掉模式（每次 run 不同）
  → mapdiv/mapmul 的 cache hit/miss 不同
  → extra_phase writeback 结果不同
  → 共享 map 的 phase 被非确定地改写
  → 收缩结果非确定 → 节点爆炸
```

### 4.5 为什么 map 表会填满（加剧问题）

`mapmulTable`/`mapdivTable` **从不 clear**。`Package.hpp:902-903` 只有 `addTable.clear()` / `contTable.clear()`（这两个表周期清空），而 map 表在整个 run 中持续累积。

---

## 5. 已尝试并排除的假设

| 假设 | 改动 | 结果 |
|---|---|---|
| 权重指针哈希 | `std::hash<Complex>` 改基于值 | ❌ 非确定性依旧（17k→476k） |
| 边相等模糊性 | `Edge::operator==` 从 `approximatelyEquals` 改精确指针 `w==other.w` | ❌ 非确定性依旧 |
| 异常被吞 | `catch(...)` 加打印 | ❌ 无 `STEP_EXCEPTION`（无异常） |
| 碰撞安全（线性探测） | `ComputeTable3` 改线性探测 | ❌ 节点数正确但 O(N) 慢（map 表填满） |
| 确定性内容哈希 | `std::hash<the_maps*>` 改哈希 (level,x,rotate) 链 | ⚠️ O(1) 恢复但节点数暴增 474k，且仍 ~4% 非确定 |
| 禁用 writeback | `LIMTDD_DISABLE_MAPDIV_LOOKUP_WRITEBACK=1` | ❌ 节点爆炸（prefix=300 从 102 → 103803），writeback 是**必要**的，不是根因 |
| **禁用 map memoization** | 注释掉 `mapmulTable.insert` / `mapdivTable.insert`（`Package.hpp:1287/1428`） | ✅ **节点数收敛**（prefix=550 三次 17166/17145/17054，±0.7%）；但 prefix=550 从 ~11s 变 ~60s |
| **删除 base 相位清零** | 注释掉 `mapmul`/`mapdiv` 的 4 处 `extra_phase = 0`（`1206/1222/1324/1338`） | ❌ **打坏 `ae_10`**（返回旧 carry，与 `cont2` 单独施加的 `extra_phase` 重复计数）→ 证明清零是「输出 carry 从 0 起点」的语义，非可删的脏数据 |

**所有改动均已回滚，代码树干净。**

关键教训：

1. **线性探测不可行**：map 表从不 clear，填满后 miss 变成 O(N) 全表扫描（prefix=400 从 0.56s 退化到 8s）。
2. **单改哈希不可行**：确定性内容哈希分布差导致碰撞变多→节点暴增；且 `contTable`/`addTable`（`Edge.hpp:96/99` 的 `murmur64(e.p)`/`murmur64(e.map)`）仍是指针哈希，是残留非确定性源。
3. **禁用 writeback 不可行**：writeback 是为了在命中时恢复 map 的 phase，禁用它会导致 phase 失配、节点爆炸。writeback 是"可变 `extra_phase`"设计的必要补丁，不是根因。
4. **禁用 memoization 收敛了节点数**：这是决定性证据——非确定性正是由 map 的 memoization 缓存引入的。
5. **base 相位清零是语义必要、不可直接删除**：`mapmul(header, other)` 返回「`other` 结构 + carry 0」，这是调用方「先捕获 carry、再单独施加」协议的组成部分；删掉会重复计数打坏 `ae_10`。正确改法是用查表（返回 phase=0 骨架）取代突变，见重构方案 §3。

---

## 6. 时间线：原始设计 + 暴露点（git 查证）

### 6.1 memoization 是原始设计（2023-10），不是后来引入的

> ⚠️ 修正：上一版曾把 `ComputeTable.hpp` 里的 `"=================我加的======================="` 标记误读成"后来 bug-fix 加的"。经 git 查证，该标记实为"LimTDD 对 MQT DD 库的扩展"。

`LimTDD/` 是独立 git 仓库（`LimTDD/.git`）。关键提交（`git log`）：

| 符号 | 首次引入提交 | 日期 |
|---|---|---|
| `ComputeTable3`（map memoization 表） | `fc7c6e0 update` | **2023-10-09** |
| `mapmulTable` / `mapdivTable` | `c23c8e1 make a new version` | **2023-10-17** |

且 `git show temp2:DDPackage/dd/ComputeTable.hpp` 确认 **temp2 分支也包含 ComputeTable3**。所以 map memoization 从 2023 年就在，Windows/temp2 实验里也有。

### 6.2 真正的暴露点：相对容差修复（temp2 → 当前唯一语义差异）

对比 `temp2` 与当前分支的 `ComplexTable.hpp`，`approximatelyEquals` 的**唯一语义差异**是：

```cpp
// temp2（旧，绝对容差）
return left == right || std::abs(left - right) <= TOLERANCE;

// 当前（新，相对容差，7fce6f4 "repair abs bug" 引入）
const auto scale = std::max(std::abs(left), std::abs(right));
return std::abs(left - right) <= TOLERANCE * scale;
```

map 运算本身（mapmul/mapdiv）在 temp2 与当前**功能等价**（temp2 用 `lookup`，当前用 `findEntry`+手动 writeback，都是同一个 writeback）。

**假说**：旧绝对容差更激进（小权重更易合并），掩盖了非确定性；新相对容差更紧（小权重不再合并）→ 更多不同权重 → 更多不同 `arg(c)` → 更多不同 `extra_phase` → 更多不同 map → direct-mapped 表碰撞更频繁 → 指针哈希的非确定性被放大显现。这对应"浮点数精度问题掩盖了此问题"的推测。

### 6.3 TDD_C 对照：相对容差有、但无 map memoization

TDD_C 的 `ComplexTable.hpp:88-101` **也是相对容差**（与 LimTDD 当前一致）。但 grep `TDD_C/include/dd/` 中 `the_maps` / `mapmul` / `mapdiv` / `mapmulTable` / `mapdivTable` / `extra_phase` **全部为空**；`TDD_C/include/dd/ComputeTable.hpp` 只有 `ComputeTable` 和 `ComputeTable2`，**没有 `ComputeTable3`**。

即 TDD_C 的相位直接以 Complex 权重存于边上（无 `the_maps` 结构、无 map memoization），所以 TDD_C 用相对容差却**没有这个非确定性来源**。这为"在保持相对容差下修复 LimTDD"提供了现成参照。

### 6.4 两层问题（不变）

1. **map 的 memoization 缓存不纯**：命中时 `extra_phase` writeback 副作用（`ComputeTable.hpp:288`、`Package.hpp:1366`），导致 hit/miss 产生不同结果。memoization 应当是纯缓存。
2. **`extra_phase` 是共享 map 上的可变字段**：被 `mapmul`/`mapdiv`/`normalize`（`Package.hpp:689/709/778/1206/1222/1244`）反复改写，而 `append_new_map`（`Package.hpp:1161`）缓存键只看 `(level,x,rotate)`、不看 `extra_phase`。这正是 `prompts_expr.md` §12 早就警告的 singleton hazard。

> 注：writeback 是"可变 `extra_phase`"设计的**必要补丁**（命中时把 map 的 phase 恢复到缓存值，避免失配），不是根因。禁用它（§5）会导致 phase 失配、节点爆炸。

---

## 7. 建议的修复方向（均可与相对容差共存）

非确定性的直接来源是 map 的 memoization 缓存（direct-mapped + 指针哈希 + 从不 clear），与相对容差正交。选项：

1. **移除 memoization**（最简、最确定，等价 TDD_C 的无 memoization 设计）：删除 `mapmulTable.insert`/`mapdivTable.insert` 及对应 `lookup`/`findEntry`，让 map 运算始终重算。代价：prefix=550 从 ~11s 变 ~60s（约 5x 慢），但结果完全确定。

2. **周期性 clear**（最小改动，2 行）：在 `Package.hpp:902-903` 的 `addTable.clear()`/`contTable.clear()` 旁加 `mapmulTable.clear()`/`mapdivTable.clear()`，避免填满 → 降低碰撞 → **降低**非确定性（不根除，有残留）。

3. **碰撞安全 + 确定性哈希**（保留提速、根治）：map 表改链地址/探测，哈希改成确定性（给 map 分配稳定整数 ID）。改动较多。

4. **根除可变 `extra_phase`**（最彻底、风险最高）：让 `extra_phase` 不再作为共享 map 的可变字段（拆出按边存储），使 memoization 变纯、hit/miss 结果一致。这是 `prompts_expr.md` §12 的 singleton hazard，动核心语义，需配合 `ae_10`/`dj_60` sanity 回归。

> **方向 4（根除可变 `extra_phase`）已选为方向，但具体落地需按「方向 2」而非「方向 1」**。方向 1（phase 纳入 `append_new_map` 键）已实施并证伪：`ae_10` 过但 `dj_60` 从 121 节点爆炸到 13GB（Clifford 抵消失效）。根因是 `extra_phase` 是「边上的 pending phase」，不是 map 结构身份；`mapmul`/`mapdiv` 从不读输入 phase 是**正确不变量**，旧代码靠别名（同结构⇒同节点⇒同 phase）让其成立，phase-in-key 打破了它。
>
> **正确方向 2**：把 phase 从 `the_maps` 结构拆出（结构键保持 `(level,x,rotate)`），`mapmul`/`mapdiv` 返回 `(结构, phase)`，phase 由边/调用方持有并在 `cos` 处消费。完整设计与逐处语义对照见 [`limtdd-extra-phase-immutable-refactor-plan.md`](./limtdd-extra-phase-immutable-refactor-plan.md) §2。教训：`dj_60` 必须加入 sanity 回归——`ae_10` 过不代表 Clifford 抵消没坏。

> 已证伪：禁用 writeback（会爆炸）、线性探测（O(N) 慢）、单改哈希（分布差/不完整）。

---

## 8. 阶段结论（2026-08-26）：方向 2 已实施，指针哈希已根除，残留为 T 门特有非确定

### 8.1 方向 2（拆出 pending phase）已实施并验证正确

按 [`limtdd-extra-phase-direction2-implementation.md`](./limtdd-extra-phase-direction2-implementation.md) 全量落地：

- `the_maps` 删除 `extra_phase` 字段，新增 `struct map_res { the_maps* map; int phase; }`。
- `mapmul`/`mapdiv` 返回 `map_res`（结构 + 溢出 phase），**不再写回节点、不再读输入 phase**（保持「结构-only」不变量）。
- memoization 变纯：命中用 `findEntry` 读 `(result, extra_phase)`，**无 writeback**；`ComputeTable3::lookup` 的 writeback 行删除。
- 调用点全改：9 处 `sliceStateEdge`（`Slicing`/`Slicing2`/`T_add2`/`cont2` + 3 个测试文件的副本）、`normalize`（`promoted_phase` 局部变量）、`find_remain_map`（`remain_phase`）、`cont2`。

验证（本轮 macOS checkout）：

```text
ae_10 e=0 -> 1；ae_10 e=1 -> 0
dj_60 e=0 -> 1（120 评估量子位，max_nodes≈121，快速）
6-qubit Clifford+T（11 个 T 门）与 Qiskit 吻合到 2.5e-16
```

### 8.2 指针哈希深挖：全部改为内容哈希

把 DD 相关全部指针哈希改成基于值/内容（消除 ASLR 依赖）：

| 文件 | 改动 |
|---|---|
| `Complex.hpp` | `hash<Complex>` 从 `murmur64(指针)` → `round(值/tolerance)` |
| `Maps.hpp` | 新增 `hash<the_maps*>`：沿 `(level,x,rotate)` 链哈希 |
| `Node.hpp` | `mNode` 加 `hash` 字段（内容哈希，跨运行确定） |
| `Edge.hpp` | `hash<Edge>`/`hash<CachedEdge>` 用 `p->hash` + 值权重 + 内容 map |
| `Package.hpp` | `makeDDNode` 计算并写入 `p->hash` |

验证：Clifford `1time.qasm` 连跑 5 次均 **11 节点**，完全确定。证明 `Complex`/`Edge`/`the_maps` 的指针 hash 非确定已根除。

### 8.3 残留：334time（Clifford+T）非确定，是 T 门特有的另一条线

数据（`LIMTDD_TN_PREFIX=500` 六次）：`2451 / 1563 / 88467 / 2451 / 4177 / 20269`。

关键诊断：

- `prefix=480` 两次运行 MAX node 都 1598、回归计数器指纹一致（仅耗时不同）。
- `prefix=500` 三次分别 2451 / 1254 / 7568，计数器（`normalize.child_phase_adds` 等）**开始发散**。

结论：

1. **指针 hash 这条线已挖干净**：Clifford（相位 ±1/±i 精确）已确定，说明 `Complex`/`Edge`/`the_maps` 的指针 hash 非确定已根除。
2. **334time 的残留非确定是「T 门特有的、更底层的另一条线」**，不是指针 hash——Clifford 已确定，而 T 门（相位 π/4，权重含 √2 无理数）才非确定。

### 8.4 残留非确定的精确定位（2026-08-26 追查）

用 `setContStageTrace` + `LIMTDD_STEP_TRACE_START/END` 逐层 trace 后，把残留非确定收敛到以下事实：

**（a）指针 hash 已彻底根除。** step 440 的逐层 trace（map 链、权重、phase、`ifContract` 判定）在 8 次运行中**内容完全一致**，只有 ASLR 地址不同。Clifford 已确定。

**（b）残留是「罕见浮点边界事件」，onset 在 step 430–440 的收缩路径**：

- 10 次运行 trace step 428–442 的**状态节点数全部一致**（step 440 恒为 1029），说明该点非确定概率 < 1/10。
- 但 prefix=500 的 3 次运行 max node `7568 / 2451 / 1254` **全不同** —— 非确定在 442→500 之间从「罕见」累积成「必然」。

**（c）机制定位到 `normalize` 的相位取整边界**（`Package.hpp:718`）：

```cpp
int rot = round(angle / rotate_angle);   // rotate_angle = π/4
```

T 门组合（T·H / H·T 等）产生相位 **π/8、3π/8** 的权重，使 `angle/rotate_angle ≈ 0.5/1.5/2.5…`，恰落在 `round()` 的**半整数边界**。权重差 1 ULP → `atan2` 差 ~1e-16 → `round()` 翻转 → `rot` 差 1 → `promoted_phase` 差 1 → map 结构变 → 节点数分叉。这是「微小浮点差在边界被放大」的经典签名。

**（d）未闭合的一环**：1-ULP 权重差的**最终来源**仍未定位。所有 hash 已是内容哈希、map/node 已去重规范化，遍历顺序理论上确定。

**（e）已证伪：DD 缓存（`addTable`/`contTable`）静默驱逐不是来源。** 给 `addTable` 加了 `LIMTDD_DISABLE_ADD_CACHE` 开关后禁用，prefix=500 仍非确定（`27055/69744/27194/27194`，且节点数普遍更大）—— 说明非确定在**核心 DD 运算**（`T_add2`/`normalize`/`makeDDNode`）里，缓存只是放大器，不是根因。`contTable` 已有 `LIMTDD_DISABLE_CONT_CACHE` 开关，禁用后慢到超时（contraction 失去 memoization）。

**（f）`round()` 半整数边界已实证存在，但对相同输入是确定的。** 新增 `LIMTDD_TRACE_ROUND_BOUNDARY` 诊断，抓到的真实例子（step 454）：

```text
ratio 0.49999999999998862  rot 0  angle 0.3926990816987152
ratio 2.4999999999999853   rot 2  angle 1.9634954084936092
```

`ratio` 距离 0.5 仅 ~1.1e-14，正是 T·H 产生的 π/8 相位落在 `round()` 边界上。但 5 次运行这些 `ratio` 的**精确比特完全一致** —— 说明边界本身是确定的，真正非确定的是**到达边界的权重**（1-ULP 差来自更上游）。

已排除的候选（本轮验证）：

- `cir_2_tn` 张量构造：`step_tensor[i]` 各步节点数确定（tensor 构建确定，非确定只在 `step_cont[i]` 收缩）。
- `ComplexCache` 的 `complexMap`（`std::hash<fp*>` 指针哈希）：其唯一读者 `isInCache` 在热路径中被注释掉，`insert`/`erase` 不影响结果。
- `ComplexTable` 本身：值哈希 + 链地址 + 有序桶，确定。
- `addTable`/`contTable` 缓存驱逐：禁用不改变非确定（见 e）。

下一步（待办）：定位 1-ULP 权重的来源。两个方向：

1. **sanitizer 构建**：`-fsanitize=undefined,address`（或 Clang `-ftrivial-auto-var-init=pattern`）抓未初始化内存读 / UB —— 罕见非确定最可能是某个默认构造的 `Edge`（`p`/`w` 无默认初始值）在极少数分支被未初始化读。
2. **在更早 step 抓罕见事件**：非确定 onset 在 step 430–440 且 <1/10 概率，需大量运行 + 逐层 trace diff 才能定位第一个 1-ULP 分叉的确切操作。

新增诊断开关（均已落地，env 门控、默认不改变行为）：`LIMTDD_DISABLE_ADD_CACHE`、`LIMTDD_TRACE_ROUND_BOUNDARY`（配合 `LIMTDD_REGRESSION_DIAG=1` + `LIMTDD_STEP_TRACE_START/END`）、`LIMTDD_DISABLE_GC`、`LIMTDD_DISABLE_CONT_CACHE`（已有）。

### 8.5 sanitizer 构建结果（2026-08-26）

做了 `-fsanitize=address,undefined -ftrivial-auto-var-init=pattern` 构建，逐项隔离后结果：

| 配置 | prefix=500 MAX node | 结论 |
|---|---|---|
| baseline（`-g -O0`） | `7030/1531/1531/13707/3586/13638` | 非确定 |
| `-ftrivial-auto-var-init=pattern` | `13638/13707/14047/3586/13638` | 仍非确定 → **非未初始化栈** |
| `MallocPreScribble=1` / `MallocScribble=1` | 仍非确定 | → **非未初始化堆** |
| `-fsanitize=undefined`（UBSan） | `13638/3586/13638/13707` | 仍非确定 |
| `-ffp-contract=off -fno-omit-frame-pointer` | `1322/1287/1322/1466/1322` | 仍非确定（节点数变小）→ **非 FMA** |
| `LIMTDD_DISABLE_GC=1`（禁 GC） | `1337/6691/9482/9482/48127` | 仍非确定 → **非 GC/缓存悬挂指针** |
| **`-fsanitize=address`（ASan）** | **`98190 × 4` 全同** | **确定！** |

**结论**：非确定**不是**未初始化内存（栈/堆）、不是 GC/缓存、不是 FMA、不是 UBSan 能抓的 UB。但 **ASan 使其确定**（地址 + redzone 布局变化）。这说明根子是**地址（ASLR）依赖**的某个操作——只是不在 hash 里（hash 已内容化），而是某种更隐蔽的、ASan 的 allocator 恰好中和掉的地址依赖。UBSan 也顺带暴露并已修复一个 hash UB：`static_cast<size_t>(负 double)`（`Complex.hpp`/`ComplexValue.hpp` 的 `round`→`llround` 修复），但它不是非确定源（修后仍非确定）。

下一步：用**逐 step 的 DD 状态 checkpoint**（把每步的节点/权重/map 按规范顺序 dump + 跨 run diff）确定性地定位第一个发散的地址依赖操作。

### 8.6 checkpoint 逐层定位：第一个分叉是 argmax 翻转（2026-08-26）

**checkpoint 工具**（`test_data.cpp` 新增 `LIMTDD_CHECKPOINT`，利用已有的 `mNode::hash` 内容哈希打印每步根边指纹）把第一个发散从「靠节点数看到的 step 440」提前到 **step 294**。

在 step 294 内做逐层诊断（`frm`=find_remain_map、`contk`=cont2 分支、`mmul/mdiv`=mapmul/mapdiv、`norm`/`normout`=normalize、`slc`=Slicing、`argmax`），跨 run diff 结果：

| 层级 | 结果 |
|---|---|
| `frm`（`newk1/newk2/ifContract`） | 452 行全同 ✅ |
| `contk`（cont2 的 `newk1/newk2/ifc1/ifc2`） | 71 行全同 ✅ |
| `mmul`/`mdiv`（mapmul/mapdiv 输入 map + 相位） | 分叉点前全同 ✅ |
| `norm`/`normout`/`slc`/`argmax` | **发散** |

**第一个确切分叉**是 `argmax`（depth 10, var 7）选了不同的子节点：

```text
run1(非发散):  maxIdx=0  maxmag2=1      c1w=-0.292893-0.292893i  (模长 0.414)
run11(发散):   maxIdx=1  maxmag2≈5.83   c1w= 1.70711+1.70711i    (模长 2.414)
```

子 0 完全相同（权重 1），子 1 权重**完全不同的模长**（0.414 vs 2.414 = `(1-1/√2)` vs `(1+1/√2)`），导致 argmax 的 `mag - max_mag2 > tolerance/2` 翻转，进而 `add_x`、`promoted_phase`、map 结构、下游 Slicing 的 `cos(rotate·π/4)` 全部连锁分叉。

**完整因果链（已闭环）**：

```
权重差 1 ULP（种子，在 T_add2 求和 / ComplexNumbers::mul·div 的边界）
  → angle=arg() 差 ~1e-16
  → rot=round(angle/rotate_angle) 在 π/8 半整数处翻转
  → promoted_phase 差
  → normalize 的 map 结构变
  → Slicing 读 e.map->rotate，cos(rotate·π/4) 符号翻转
  → 子权重差 ~5.83 倍（被放大成结构差异）
  → argmax 翻转（maxIdx 0↔1）
  → DD 结构继续分叉 → 非 canonical → 节点爆炸
```

**为什么影响效率**：DD 失去唯一性（canonicity）→ 语义等价但结构不同的节点无法被 uniqueTable 去重 → 节点数增长 → direct-mapped 缓存（`addTable`/`contTable`）content-hash 因结构不同而大量 miss → 反复重算 → GC 频繁 → 节点爆到 8 万+、耗时 5 分钟+。

**一个被证伪的修复尝试**：给 `round(angle/rotate_angle)` 加「半整数（π/8）吸附」——跑 prefix=500 仍非确定（`74017/1396/58284/1396/70154/73889/1396/73889`），**已回退**。这说明分叉**不在 π/8 半整数边界**，而在**整数边界（因子 -1 或 -i 的全局相位歧义）**，是一个自洽闭环：

```
child.w 差 -1（全局相位 π）→ c 差 -1 → angle 差 π → rot 差 4 → promoted_phase 差 4
→ map.rotate 差 4 → Slicing cos 差 -1 → 子权重再差 -1 → ……（±1 自洽闭环）
```

这个 -1/-i 是**物理不可观测的全局相位**，但 DD 把它当结构存下来，导致非 canonical。

**下一方向（未做）**：不再追 1 ULP，而是做**全局相位规范化**——让 normalize 的 `res.w`（根边权重）强制落入一个固定相位范围（如 `arg(res.w) ∈ [0, π/4)`），使同一物理态总是映射到同一个 DD。这才是「消除连锁反应」的正解。

