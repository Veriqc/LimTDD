# LimTDD 方向 2 详细实施文档：拆出 pending phase

> **日期**: 2026-08-23（实施完成 2026-08-26）
> **状态**: **已实施并验证正确**。方向 1 已证伪，见 [`limtdd-extra-phase-immutable-refactor-plan.md`](./limtdd-extra-phase-immutable-refactor-plan.md)。实施结果与残留非确定见 [`limtdd-cliffordt-334time-nondeterminism.md`](./limtdd-cliffordt-334time-nondeterminism.md) §8。
> **结论先行**: phase 是「边上的 pending phase」，`mapmul`/`mapdiv` 从不读输入 phase 是**正确不变量**。方向 2 把它从 `the_maps` 结构拆出，作为 `mapmul`/`mapdiv` 的**返回值**，由调用点立即消费（折进权重）或短暂 carry（local 变量 / `comm_maps.remain_phase`）。**不需要给 `Edge` 加 `phase` 字段。**

---

## 1. phase 生命周期（本次深读的关键结论）

全库 grep `->extra_phase`，把所有读点按「是否紧跟 `mapmul`/`mapdiv`」分类：

### 1.1 模式 A：产生即消费（折进权重）——占绝大多数

读点**全部**紧跟一个 `mapmul`/`mapdiv` 调用，phase 是「结构组合产生的溢出」，立即折进权重：

| 读点 | 前一条 map 操作 | 函数 |
|---|---|---|
| 1755/1771/1794 | `temp.map = mapmul(...)` | `Slicing`（供 add） |
| 1886/1915/1952 | `temp.map = mapmul(...)` | `Slicing2`（供 contract） |
| 2206 | `yCopy.map = mapdiv(...)` (2192) | `T_add2` |
| 2253 | `temp_map = mapmul(...)` (2250) | `T_add2` |
| 2342 | `e.map = mapmul(...)` (2338) | `T_add2` |
| 2856 | `e.map = mapmul(...)` (2853) | `cont2`（cache hit） |
| 3641 | `r.map = mapmul(...)` (3638) | `cont2`（cache miss） |

### 1.2 模式 B：产生后 carry 一小段，再消费——只有 2 处

1. **`normalize`**：`mapdiv` 溢出 + 权重相位 `rot` → 存到 `res.p->e[1].map->extra_phase`（778/709）→ 在 849 行**晋升为父节点 rotate**。这是唯一「跨多行 carry」的 phase。
2. **`find_remain_map`**：`remain_map->extra_phase` 累加（2433）→ `cont2` 2722 捕获 → 2870/3645 折进权重。已在方向 1 探索中改为 `comm_maps.remain_phase`，本次直接采用。

### 1.3 关键推论

- **`Edge` 上的 map 的 `extra_phase` 从不被「下游直接读」**——下游总是把它喂给 `mapmul`/`mapdiv`（后者忽略输入 phase），读到的是**新产生**的溢出。所以「carry」只在模式 B 的两处发生。
- 因此 phase 可以用 **`mapmul`/`mapdiv` 的返回值 + local 变量 + `comm_maps.remain_phase`** 承载，**无需 `Edge.phase`**。

---

## 2. 设计

### 2.1 数据结构

`Maps.hpp`：
```cpp
struct the_maps {
    short level;
    bool  x;
    int   rotate;                       // 结构相位，[0,8)
    std::map<std::string, the_maps*> next;
    the_maps* father;
    // ❌ 删除 int extra_phase;
    static the_maps the_maps_header_element;
    static constexpr the_maps* the_maps_header() { return &the_maps_header_element; }
};

// mapmul/mapdiv 的返回：结构 + 溢出 phase
struct map_res {
    the_maps* map;
    int phase;
};
```

`the_maps_header_element` 初始化为 `{ -1, 0, 0, {}, nullptr }`（去掉第 4 个字段）。

### 2.2 `mapmul`/`mapdiv` 返回 `map_res`

- 结构部分照旧（`append_new_map` 键回到 `(level,x,rotate)`）。
- phase 部分 = 原 `res->extra_phase` 的计算值（递归里 `r_phase + local_phase`），**作为返回值返回**，不再写回节点。
- **不再读输入 phase**（保持「结构-only」不变量）；输入 phase 由调用方负责（模式 B 里 carry）。

### 2.3 memoization 变纯

`ComputeTable3` 已经用 `Entry { leftOperand, rightOperand, result, extra_phase }` 存了 phase。方向 2：

- `mapmul`/`mapdiv` 命中时用 `findEntry` 读 `entry->result`（结构）+ `entry->extra_phase`（phase），**返回 `map_res`，无 writeback**。
- `ComputeTable3::lookup`（含 writeback 的旧方法）删除或停用；`mapmul` 从 `lookup` 改用 `findEntry`。
- `insert` 存 `(self, other, structure, phase)`（phase 归一化到 `[0,8)`）。

---

## 3. 逐文件改动清单

### 3.1 `Maps.hpp`

- `the_maps` 删 `int extra_phase;`。
- 新增 `struct map_res { the_maps* map; int phase; };`。
- `comm_maps` 保留 `int remain_phase;`（方向 1 探索已加，本次为正式设计）。

### 3.2 `Maps.cpp`

- `the_maps_header_element` 改为 `{ -1, 0, 0, {}, nullptr }`（4 字段）。
- `phase_carrier` 删除（方向 1 的产物，方向 2 不需要——phase 不再需要「level -1 带 phase」的载体，直接作为 int 返回）。

### 3.3 `ComputeTable.hpp`

- `ComputeTable3::lookup` 删除 writeback 行（或整个方法停用）。
- `findEntry` 保持（返回 `const Entry*`，`mapmul`/`mapdiv` 用它读 `result` + `extra_phase`）。

### 3.4 `Package.hpp` — `append_new_map`（~1153）

```cpp
the_maps* append_new_map(the_maps* self, short level, bool x, int rotate) {
    rotate = (rotate % root_of_unit + root_of_unit) % root_of_unit;
    if (x == 0 && rotate == 0) return self;
    std::string new_key = std::to_string(level) + "_" + std::to_string(x) + "_" + std::to_string(rotate);
    auto it = self->next.find(new_key);
    if (it != self->next.end()) return it->second;
    self->next[new_key] = new the_maps{ level, x, rotate, {}, self };
    return self->next[new_key];
}
```

### 3.5 `Package.hpp` — `mapmul`（~1179）

```cpp
map_res mapmul(the_maps* self, the_maps* other) {
    // (trace 不变，去掉读 extra_phase 的诊断)
    if (self->level == -1) {              // base：结构=other，phase=0
        return { other, 0 };
    }
    if (other->level == -1) {
        return { self, 0 };
    }
    if (const auto* entry = mapmulTable.findEntry(self, other); entry != nullptr) {
        return { entry->result, entry->extra_phase };   // 纯命中，无 writeback
    }
    the_maps* res;
    int phase = 0;
    if (self->level > other->level) {
        auto r = mapmul(self->father, other);
        res = append_new_map(r.map, self->level, self->x, self->rotate);
        phase = r.phase;
    }
    else if (self->level < other->level) {
        auto r = mapmul(self, other->father);
        res = append_new_map(r.map, other->level, other->x, other->rotate);
        phase = r.phase;
    }
    else {
        auto r = mapmul(self->father, other->father);
        auto rotate = other->x ? (other->rotate - self->rotate) : (other->rotate + self->rotate);
        int nr = (rotate % root_of_unit + root_of_unit) % root_of_unit;
        int local = (other->x ? self->rotate : 0);
        if ((self->x + other->x) % 2 == 0 && nr == 0) {
            // 结构抵消：只有 phase 增量，挂在父结构上
            res = r.map;
            phase = r.phase + local;
        } else {
            res = append_new_map(r.map, self->level, (self->x + other->x) % 2, nr);
            phase = r.phase + local;
        }
    }
    phase = (phase % root_of_unit + root_of_unit) % root_of_unit;
    mapmulTable.insert(self, other, res, phase);
    return { res, phase };
}
```

> ⚠️ 注意「结构抵消」分支（`x==0 && nr==0`）：旧代码 `append_new_map` 会早退返回 `r`，phase 记在 `r` 上。方向 2 里结构仍返回 `r.map`，phase 单独 `r.phase + local`。这一支在方向 1 里曾漏掉、导致 `dj_60` 爆炸，必须保留。

### 3.6 `Package.hpp` — `mapdiv`（~1297）

与 `mapmul` 对称，phase 公式用减法：

```cpp
map_res mapdiv(the_maps* self, the_maps* other) {
    if (other->level == -1) {
        return { self, 0 };
    }
    if (self == other) {                    // 同指针 → 同结构同 phase
        return { the_maps::the_maps_header(), 0 };
    }
    if (const auto* entry = mapdivTable.findEntry(self, other); entry != nullptr) {
        return { entry->result, entry->extra_phase };   // 纯命中
    }
    // 递归分支：结构用 append_new_map，phase = r.phase + local
    //   self->level >  other->level: phase = r.phase
    //   self->level <  other->level: other->x==0 → phase=r.phase；other->x==1 → phase=r.phase-other->rotate
    //   equal: x=(sx+ox)%2；nr；local=(x==1 ? -other->rotate : 0)
    //          x==0&&nr==0 → res=r.map, phase=r.phase+local；else append_new_map(...)
    ...
    mapdivTable.insert(self, other, res, phase);
    return { res, phase };
}
```

### 3.7 `Package.hpp` — 模式 A 调用点（折进权重）

统一把 `X.map = mapmul(...)` + `cn.mul(..., cos(X.map->extra_phase*rotate_angle))` 改为：

```cpp
auto mr = mapmul(a, b);
X.map = mr.map;
if (X.w != Complex::zero)
    cn.mul(X.w, X.w, cn.getTemporary(cos(mr.phase*rotate_angle), sin(mr.phase*rotate_angle)));
```

涉及：`Slicing` 1751/1755、1767/1771、1790/1794；`Slicing2` 1885/1886、1909/1915、1950/1952；`T_add2` 2250/2253、2338/2342、`mapdiv` 2192/2206；`cont2` 2853/2856、3638/3641。

> ⚠️ **还有 3 个测试文件的 `sliceStateEdge` 副本**（`Slicing` 的复制品）也调 `mapmul` + 折 `extra_phase`，方向 2 必须同步改：
> - `test/test_fidelity.cpp`：289/290、299/300、311/312（3 处 `mapmul` + `cos`）
> - `test/test_state_output.cpp`：119/120、129/130、141/142（3 处）
> - `test/test_fidelity_trace_tn.cpp`：223/224、233/234、245/246（3 处）
>
> 全库 grep 确认：**没有**「直接读存量 edge 的 map->extra_phase、无前置 mapmul/mapdiv」的读点。`edge.map->rotate`（结构 rotate）读点不受影响。

> 注：`Slicing`/`Slicing2` 里若 `c==1` 或 `c==0` 还有 `cos(e.map->rotate*rotate_angle)`（读的是**结构 rotate**，不是 extra_phase），这部分不变。

### 3.8 `Package.hpp` — `normalize`（模式 B，carry → 晋升）

核心：用 local 变量 `promoted_phase` 承载非 max 子边的 phase，849 行晋升。

```cpp
int promoted_phase = 0;
for (i = 0; i < 2; ++i) {
    if (i == maxArgIndex) {
        res.map = res.p->e[i].map;
        res.p->e[i].w = Complex::one;
        res.p->e[i].map = the_maps::the_maps_header();      // phase 0
    } else {
        if (isZero[i]) { /* 零子边：phase 0，promoted_phase 不变 */ }
        else {
            // c = child.w / max_w; rot/angle 计算不变
            auto mr = mapdiv(res.p->e[i].map, res.map);
            res.p->e[i].map = mr.map;
            promoted_phase = mr.phase + rot;   // (mode==2 用 int(arg/rotate_angle) 替代 rot)
            res.p->e[i].w = ...;               // 残差权重不变
        }
    }
}
res.map = append_new_map(res.map, res.p->v, add_x, promoted_phase);   // 原 849
```

- 删除 689 的 `map->extra_phase = 0`（零子边 phase 恒 0）。
- 709/778 的 `map->extra_phase += ...` 变成 `promoted_phase = mr.phase + ...`。
- 849 读 `promoted_phase` 而非 `res.p->e[1].map->extra_phase`。

> ⚠️ 务必确认 `add_x`/swap 之后 `res.p->e[1]` 恒为「非 max 子边」——经核对（§1 之前读码），`add_x = (maxArgIndex>0)`，swap 后 `res.p->e[1]` 始终落到非 max 子边，故 `promoted_phase` 就是它。

### 3.9 `Package.hpp` — `find_remain_map`（模式 B，carry）

沿用方向 1 探索已定稿的拆分（本次正式采用）：

- `comm_maps{ header, header, header, 0 }`（`remain_phase=0`）。
- 2433 `res->remain_phase += map2->rotate`。
- 删 `temp_pahse` 保存/恢复（2381-2396 的 6 行）。
- `cont2` 2722 `extra_phase = r_maps->remain_phase`。

### 3.10 `Package.hpp` — `cont2`（模式 A + remain_phase）

- 2853/2856：`auto mr = mapmul(remain_map, e.map); e.map = mr.map; cn.mul(e.w, cos(mr.phase))`。
- 2870/3645 的 `cos(extra_phase)` 不变（`extra_phase = remain_phase`，已是 int）。
- 3638/3641 同 2853/2856。

---

## 4. 实施顺序与回归

建议顺序（每步可编译、可回归）：

1. `Maps.hpp`/`Maps.cpp`：删 `extra_phase`、加 `map_res`、改 header 初始化。
2. `Package.hpp`：`append_new_map` 回退 + `mapmul`/`mapdiv` 改返回 `map_res`（含 `ComputeTable3` 命中改 `findEntry`）。
3. `Package.hpp`：改所有模式 A 调用点（Slicing/Slicing2/T_add2/cont2）。
4. `Package.hpp`：改 `normalize`（`promoted_phase`）。
5. `Package.hpp`：改 `find_remain_map`（`remain_phase`）。
6. `ComputeTable.hpp`：删 writeback。

每步后回归（务必含 `dj_60`，方向 1 的教训）：

```bash
cmake --build LimTDD/build -j2
./LimTDD/build/test/test_fidelity expr/qasm/ae_10.qasm 0 1   # 1
./LimTDD/build/test/test_fidelity expr/qasm/ae_10.qasm 1 1   # 0
./LimTDD/build/test/test_fidelity expr/qasm/dj_60.qasm 0 1    # 1，且 max_nodes≈121、快速
# 最终：334time 确定性（连跑 5 次 max node 收敛）
LIMTDD_TN_PREFIX=550 ./LimTDD/build/test/test_data LimTDD/Benchmark/cliffordT0.02big20/334time.qasm 00000000000000000000
```

---

## 5. 风险与待验证点

1. **`normalize` 的 `promoted_phase`**：需确认 swap 后 `res.p->e[1]` 恒为非 max 子边（§3.8 已核对，实施时用 `ae_10`/`dj_60` 回归确认）。
2. **`mapmul`/`mapdiv` 命中/未命中返回的 phase 一致性**：命中返回 `entry->extra_phase`（mod 8），未命中归一化 `phase % 8`，两者一致（cos 周期 8）。实施时用 `dj_60` 回归确认无意外依赖未归约 phase。
3. **`Slicing`/`Slicing2` 里 `e.map->rotate` 的读**（结构 rotate）与 extra_phase 无关，不要误改。
4. **`addTable`/`contTable` 的 `CachedEdge` 含 `map` 指针**：map 结构去重键不变，`CachedEdge` 无需改动。

> 本方案改动约 20 处、跨 3 个文件，全部是「把 phase 从节点字段改成返回值/int」的机械变换，不改变相位代数。方向 1 已证伪的核心教训——`extra_phase` 是 pending phase、结构去重键必须保持 `(level,x,rotate)`——在本方案中成立。
