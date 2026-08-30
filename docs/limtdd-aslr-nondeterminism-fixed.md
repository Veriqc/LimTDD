# LimTDD ASLR 非确定性：根因与修复（节点 ID + 表项 ID）

> **日期**: 2026-08-30
> **状态**: 非确定性已根治；节点数在 prefix≤500 恢复；prefix≥550 仍有 ±1/-i 非规范性膨胀（独立于非确定性，见 [`limtdd-cliffordt-334time-nondeterminism.md`](./limtdd-cliffordt-334time-nondeterminism.md) §8.6）。

---

## 0. 结论一句话

`334time.qasm` 的非确定性由**两处 ASLR 依赖**引起，均已用「确定性创建序号 ID」替换指针地址而修复：

| 源 | 位置 | 修法 |
|---|---|---|
| DD 加法操作数顺序 | `Package.hpp` `T_add2` 的 `if (x.p > y.p)` | `mNode.id`（`UniqueTable::getNode` 的分配序号） |
| 权重哈希桶分布 | `hash<Complex>` 的 `reinterpret_cast(指针)` | `ComplexTable::Entry.id`（`getEntry` 的分配序号） |

关键认识：**指针地址是 ASLR 随机的**，任何「按地址比较/哈希」的代码都会让结果跨运行漂移。之前的「内容哈希」方向（`llround(值/tolerance)`）虽然确定了，但用**绝对容差**和 `approximatelyEquals` 的**相对容差**不一致，导致大权重时哈希分得过细 → 唯一表去重漏掉 → 节点爆炸。**entry ID 是既确定又和去重一致**的（近似相等的值共享同一个 ComplexTable 表项 → 同一个 ID），因此同时拿到「确定」与「不爆炸」。

---

## 1. 两个 ASLR 依赖源

### 1.1 `T_add2` 的指针比较（操作数顺序 → 1-ULP 权重差）

`Package.hpp` 的 DD 加法 `T_add2` 原本用：

```cpp
if (x.p > y.p) {
    return T_add2(y, x);   // 规范化操作数顺序
}
```

`x.p`/`y.p` 是节点指针，地址 ASLR 随机。跨运行时 `x.p > y.p` 翻转 → 操作数顺序不同 → **浮点求和顺序不同 → 1-ULP 权重差**。这正是 §8.4(d) 一直找不到的「1-ULP 最终来源」，也是 §8.5 观察到「ASan 使其确定（改变 allocator 恰好中和了地址依赖）」的原因。

### 1.2 `hash<Complex>` 的指针哈希（桶分布 → 缓存命中/未命中模式）

`hash<Complex>` 原本 `murmur64(reinterpret_cast<std::size_t>(c.r))`，即对表项**指针地址**哈希。桶分布随 ASLR 漂移 → `addTable`/`contTable`/唯一表的缓存命中/未命中模式不同 → 缓存值与重算值差 1-ULP → 非确定。

---

## 2. 为什么「内容哈希」方向会爆炸（已被证实并绕开）

之前的内容哈希把 `hash<Complex>` 改成：

```cpp
std::llround(dd::CTEntry::val(c.r) / tolerance())
```

这是**绝对容差网格**。但 `approximatelyEquals` 用的是**相对容差** `TOLERANCE * max(|a|,|b|)`。两者不一致：

- 大权重（|w|>1）时相对容差比绝对网格宽 → 两个「近似相等」的值落在不同网格点 → 哈希不同 → 唯一表跨桶 → 去重漏掉 → **节点爆炸**（prefix=500 从 ~3500 涨到 13638，prefix=550 涨到 ~411k）。
- 小权重时则反过来（哈希太粗，但等式检查会兜住假命中）。

所以「内容哈希」在 334time 上换来的是「确定但爆炸」。这不是可选权衡，而是**哈希与相等性不一致**这个真实缺陷。

---

## 3. 修复：确定性的创建序号 ID

思路与「给 `mNode` 分配序号」完全同构，只是作用对象换成 `ComplexTable` 表项：

```cpp
// ComplexTable::Entry 增加 id
struct Entry {
  fp value{};
  Entry* next{};
  RefCount refCount{};
  std::size_t id = 0;   // 唯一创建序号
};

// getEntry() 里分配（两个路径统一）
entry->id = ++nextEntryId;
```

`hash<Complex>` 改为用表项 ID（并保留符号位，因为 +v 与 −v 是不同值）：

```cpp
auto* er = dd::CTEntry::getAlignedPointer(c.r);   // 解码符号位
auto h1 = dd::murmur64(static_cast<std::size_t>(er->id));
h1 = dd::combineHash(h1, dd::CTEntry::isNegativePointer(c.r) ? 1U : 0U);
// i 同理
```

配套改动：

- `hash<Edge>`/`hash<CachedEdge>`：`murmur64(e.p)`（指针）→ `e.p->id`（节点序号）。
- `hash<the_maps*>`：保留内容哈希（`(level,x,rotate)` 链是**精确**的，非容差取整，本身一致且确定）。
- `T_add2`：`x.p > y.p` → `x.p->id > y.p->id`。

为何一致：`ComplexTable::lookup` 把近似相等的值去重到**同一个表项**，所以「近似相等 ⇒ 同一表项 ⇒ 同一 ID ⇒ 同一哈希」，与 `approximatelyEquals` 完全对齐；同时 ID 是创建顺序，跨运行确定。

---

## 4. 验证

### 4.1 正确性

| 用例 | 结果 |
|---|---|
| `ae_10` e=0 / e=1 | fidelity 1 / 0 ✅ |
| `dj_60` e=0 | fidelity 1，max_nodes 121 ✅ |
| Clifford `1time` | 11 节点 ✅ |

### 4.2 确定性（`test_data` 连跑，MAX node）

| prefix | 原始（非确定） | 内容哈希（确定但爆炸） | **本修复（确定且低）** |
|---|---|---|---|
| 450 | 1171/1137/1171 | 1183 | **1183 ×5 全同** |
| 500 | 1135/1554/1900 | 13638 | **5067 ×5 全同** |
| 550 | 17489/18090/28303 | 411854 | **402089 ×2 全同**（确定但膨胀） |

450/500 已恢复「确定 + 低节点数」。550 的 402089 是**确定的**（非非确定），是 §8.6 的 ±1/-i 非规范性膨胀，见下。

---

## 5. 遗留：±1/-i 全局相位非规范性（独立问题）

prefix≥550 的节点膨胀是 §8.6 描述的「±1/-i 全局相位歧义」：DD 把物理不可观测的全局相位（`map.rotate=4` 即 −1）当**结构**存下来，导致同一物理态映射到不同 DD、无法去重、节点膨胀。这是**正确但非最小表示**，与本次修复的非确定性正交。

尝试过「全局相位规范化」（把相位折进/折出 map）均打坏 fidelity（`dj_60` → 0 或 3.3e35），已回退。正确修法需要区分 map 的结构相位与标量全局相位，是下一个独立课题。
