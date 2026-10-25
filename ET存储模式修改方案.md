# ET（Esoteric Twist）存储模式修改方案

## 1. 文档目标

本文只描述把当前双数组传播方式改造成 Esoteric Twist（ET）单数组原位传播的实施方案，不直接修改源码。

改造目标：

1. 只保存一套D2Q9分布函数。
2. 把碰撞与传播合并到一次格点更新中。
3. 使用对向分布的“扭转写回”完成传播。
4. 保持当前 `main()` 时间步骨架、颗粒结构体数组和物理公式。
5. 正确适配固定壁面、压力边界、多个移动椭圆和新生流体格点。

## 2. 参考依据

ET的核心是：碰撞前从一组本地/邻居地址读取分布，碰撞后把结果写回“读取到的反方向分布所在地址”，最后交换对向方向数组的指针。这样不需要第二个格点数组，也不需要单独传播整场数据。

主要参考：

1. Geier, Schönherr, *Esoteric Twist: An Efficient in-Place Streaming Algorithmus for the Lattice Boltzmann Method on Massively Parallel Hardware*, 2017, DOI: [10.3390/computation5020019](https://doi.org/10.3390/computation5020019)。
2. Lehmann, *Esoteric Pull and Esoteric Push: Two Simple In-Place Streaming Schemes for the Lattice Boltzmann Method on GPUs*, 2022, DOI: [10.3390/computation10060092](https://doi.org/10.3390/computation10060092)。其附录给出了D2Q9 ET的具体地址形式。

## 3. 为什么ET比当前模式复杂

当前代码使用AoS形式：

```cpp
double (*Dist)[DY][DQ];
double (*Temp)[DY][DQ];
```

同一个格点的9个方向连续存放，传播结果写入第二套数组。

ET最清晰的实现方式是SoA形式：

```cpp
double * Dist[DQ];
```

每个 `Dist[f]` 是方向 `f` 的整张二维场。这样一个时间步结束后，可以直接交换对向方向的指针：

```cpp
Swap( Dist[f], Dist[Re[f]] );
```

因此ET方案比AA方案对数据布局的修改更大，但每个时间步可以使用同一套更新逻辑，不必维护两个完全不同的奇偶碰撞内核。

## 4. 推荐的数据布局

### 4.1 一维展开

建议增加：

```cpp
inline int Id( int i, int j )
{
    return i*DY+j;
}
```

分布函数改为：

```cpp
double * Dist[DQ];
```

初始化时：

```cpp
for( int f=0; f<DQ; ++f )
    Dist[f] = new double [DX*DY];
```

访问方式由：

```cpp
Dist[i][j][f]
```

改为：

```cpp
Dist[f][Id(i,j)]
```

### 4.2 不再分配 `Temp`

ET正式版本不需要 `Temp`。但实施初期建议保留原数组版本作为对照，待结果一致后再删除。

## 5. 对向方向分组

当前D2Q9使用以下对向组：

| 正向代表 `f` | 方向 | 反方向 `Re[f]` | 方向 |
| ---: | --- | ---: | --- |
| 1 | 东 | 3 | 西 |
| 2 | 北 | 4 | 南 |
| 5 | 东北 | 7 | 西南 |
| 6 | 西北 | 8 | 东南 |

定义：

```cpp
const int HalfQ[4] = { 1, 2, 5, 6 };
```

每一组只处理一次，静止方向0单独处理。

## 6. ET单步的基本算法

以下是适合当前D2Q9编号的、保持原代码风格的概念流程。

对格点 `(i,j)`：

```cpp
F[0] = Dist[0][Id(i,j)];

for( int h=0; h<4; ++h )
{
    int f = HalfQ[h];
    int r = Re[f];
    int ii = i+Ex[f];
    int jj = j+Ey[f];

    F[f] = Dist[f][Id(i,j)];
    F[r] = Dist[r][Id(ii,jj)];
}

// 使用F计算宏观量、平衡分布和碰撞后结果Post；

Dist[0][Id(i,j)] = Post[0];

for( int h=0; h<4; ++h )
{
    int f = HalfQ[h];
    int r = Re[f];
    int ii = i+Ex[f];
    int jj = j+Ey[f];

    Dist[r][Id(ii,jj)] = Post[f];
    Dist[f][Id(i,j)]   = Post[r];
}
```

完成整场更新后交换对向方向指针：

```cpp
Swap( Dist[1], Dist[3] );
Swap( Dist[2], Dist[4] );
Swap( Dist[5], Dist[7] );
Swap( Dist[6], Dist[8] );
```

这里的“扭转”就是：

- 沿正方向的碰撞结果写到邻居的反方向数组。
- 沿反方向的碰撞结果写到本地的正方向数组。
- 最后交换两个方向数组的名字。

传播被包含在读取位置、写回位置和指针交换中。

## 7. 线程安全原则

ET能够原位并行的关键不是“随便覆盖”，而是：

> 一个格点更新只写回本次更新已经读取过的那些存储位置，并且这些位置不会由另一个格点同时负责。

实现时必须逐个方向证明：

1. 每个物理分布在一个时间步只被读取一次。
2. 每个物理存储位置只被一个线程写入。
3. 对向指针交换只能在所有格点更新完成以后执行。

OpenMP循环结束形成隐式同步点。不能在并行循环内部交换方向指针。

## 8. 建议新增的辅助函数

为避免边界和颗粒代码到处重复ET地址公式，建议增加普通C风格函数：

```cpp
int Id( int i, int j );
void LoadNodeET( int i, int j, double F[DQ] );
void StoreNodeET( int i, int j, const double Post[DQ] );
double GetDistET( int i, int j, int f );
void SetDistET( int i, int j, int f, double Value );
void SwapEtDirections();
void GlobalCollideET();
```

其中：

- `LoadNodeET()` 按ET规则从本地和指定邻居装入9个物理方向。
- `StoreNodeET()` 把碰撞后结果扭转写回。
- `GetDistET()/SetDistET()` 为边界和颗粒模块提供逻辑方向访问。
- `SwapEtDirections()` 在整场更新完成后交换四组对向指针。

## 9. 各文件修改规划

### 9.1 `common.h`

计划增加：

- `HalfQ[4]` 对向组代表方向。
- `Id(i,j)` 一维索引函数。
- 必要的周期下标函数。

保持不变：

- `Ex、Ey、Re、Alpha`。
- `GridData、ParticleData、LinkData`。
- `Feq()`。

### 9.2 `Migration.cpp`

计划修改：

1. 将 `Dist/Temp` 改成9个方向数组指针。
2. 按方向分配、清零和释放内存。
3. `SetFlowField()` 使用SoA地址初始化分布。
4. `GlobalCollide()` 内部改为ET读、碰撞、扭转写回和指针交换，但保留函数名。
5. `BoundaryCollide()` 改用ET逻辑访问函数。
6. `InletAndOutlet()` 先装载逻辑分布到局部数组，修改后再编码写回。
7. `MacroCalculate()` 使用 `LoadNodeET()` 取得当前格点的9个物理分布。
8. `LatticeRemesh()` 对9个方向数组分别执行原有x方向移动。

### 9.3 `Particle.cpp`

计划修改：

1. `ParticleCollide()` 使用 `GetDistET()/SetDistET()` 执行插值反弹。
2. `ParticleForce()` 从ET逻辑分布中读取动量交换所需值。
3. `ParticleGeometry()` 通过局部数组重建新生流体节点，再用ET写回函数编码。
4. 清除颗粒内部节点时，不能假定9个物理方向都位于同一格点地址。
5. 多颗粒编号、椭圆内部判断和 `GetDelta()` 几何计算保持不变。

## 10. 主程序骨架保持方式

主循环仍保持：

```cpp
GlobalCollide();
BoundaryCollide();
InletAndOutlet();
ParticleCollide();
ParticleForce();
ParticleGeometry();
LatticeRemesh();
MacroCalculate();
```

`GlobalCollide()` 内部完成ET碰撞传播和对向指针交换。后续函数都通过ET逻辑访问层解释当前分布。

不要把边界或颗粒物理公式直接混入 `main()`。

## 11. ET与边界条件

ET改造最困难的不是BGK碰撞，而是边界地址映射。

### 11.1 固定壁面

论文中的简单反弹可以在装载阶段把指向壁面的读取重定向到本地反方向槽位。

本项目还支持BFL和LL2插值，需要读取多个相邻位置，因此建议第一版保留 `BoundaryCollide()` 独立函数，通过ET逻辑访问层翻译每次读取和写入。

完成正确性验证后，再考虑把HBB直接合并到 `LoadNodeET()`。

### 11.2 压力入口出口

压力边界公式不变，但 `Dist[0][j]` 不再是一段包含9个物理方向的连续数据。

需要：

```text
LoadNodeET
    ↓
在局部F[9]中计算入口/出口未知分布
    ↓
逐方向SetDistET
```

### 11.3 移动椭圆边界

`LinkData.f` 继续表示物理方向，不能改成ET物理槽位编号。

`ParticleCollide()` 中的二次插值公式和移动壁面速度修正保持不变，只替换数据访问方式。

`ParticleForce()` 必须在边界处理完成后读取同一个逻辑时间层的分布函数。

## 12. 新生流体节点和移动几何

当前新生流体节点会重建9个方向，并同时写入 `Dist` 和 `Temp`。

ET模式下建议：

1. 先把外推结果保存在局部 `F[DQ]`。
2. 计算 `Den、Vx、Vy`。
3. 通过 `SetDistET()` 写入9个逻辑方向。
4. 完成全部新生节点后再更新格点类型。

ET的物理分布可能存放在邻居的对向数组中。因此直接清空：

```cpp
Dist[f][Id(i,j)] = 0;
```

不一定等于清空格点 `(i,j)` 的逻辑方向 `f`。清空颗粒节点也必须经过地址转换。

## 13. 网格重映射

当前 `LatticeRemesh()` 直接复制整行 `Dist[i]`。改成SoA后，应对每个方向分别移动：

```cpp
for( int f=0; f<DQ; ++f )
{
    // 按原来的x方向范围移动 Dist[f]；
}
```

注意：

- 重映射前后不能执行对向指针交换。
- 所有方向必须移动相同的x距离。
- 颗粒中心和 `Link[n].i` 继续按原方法同步平移。

## 14. 初始化、输出和释放

### 初始化

初始时每个方向数组直接写入 `Feq()`。第一次ET更新前，方向指针处于标准配对状态。

### 输出

当前颗粒轨迹输出不直接依赖分布布局，可以保留。

如果输出完整流场，必须由 `LoadNodeET()` 解码后再写文件。

### 释放

方向指针在运行中会交换，但始终仍是最初9块内存的一个排列。结束时可以逐个释放当前 `Dist[f]`，前提是每个指针仍然唯一且没有另存重复别名。

更稳妥的方案是另存：

```cpp
double * DistMemory[DQ];
```

只用 `DistMemory` 负责释放，`Dist` 负责运行时换向。

## 15. 实施顺序

建议分八个提交实施：

1. **增加SoA方向数组和索引函数**，暂不替换原算法。
2. **增加ET装载、写回和指针交换函数**。
3. **实现无边界的ET碰撞传播**。
4. **适配上下壁面和压力入口出口**。
5. **适配宏观量计算**。
6. **适配多颗粒边界、受力和新生流体**。
7. **适配网格重映射和输出**。
8. **对照验证通过后删除 `Temp` 和旧传播路径**。

## 16. 验证方案

### 16.1 地址级测试

在小网格中给每个 `(i,j,f)` 写入唯一编号，只执行一次ET装载、写回和指针交换，检查每个编号是否传播到预期位置。

分别测试：

- 东西方向。
- 南北方向。
- 两组对角方向。
- x周期回绕。
- 上下壁面重定向。

### 16.2 无颗粒泊肃叶流

比较ET和原双数组版本的：

- 密度场。
- 速度场。
- 总质量。
- 解析解误差。

### 16.3 固定椭圆

令颗粒不运动，比较 `LinkNum、d、Fx、Fy、ParTor`。

### 16.4 运动单颗粒和多颗粒

逐步比较每个颗粒的位置、速度、方向角、受力和力矩，最后比较完整输出文件。

### 16.5 并行一致性

分别使用1线程和多线程运行，确认结果可重复，排除原位写入竞争。

## 17. 验收标准

1. 只保留一套D2Q9分布数据。
2. 一个时间步只读取和写入每个有效分布一次。
3. 对向方向指针只在全场更新完成后交换。
4. 单线程和OpenMP多线程结果一致。
5. HBB、BFL、LL2和压力边界均可运行。
6. 移动椭圆和新生流体节点不产生NaN。
7. 多颗粒受力与轨迹不会串位。
8. 与原双数组版本在允许浮点误差内一致。

## 18. 主要风险

| 风险 | 后果 | 控制方法 |
| --- | --- | --- |
| 对向指针交换时机错误 | 整场方向解释颠倒 | 只允许 `GlobalCollide()` 末尾统一交换 |
| SoA迁移遗漏某个函数 | 局部结果错误或访问旧布局 | 全仓库搜索全部 `Dist` 直接访问 |
| 对角方向映射错误 | 产生非物理横向速度 | 对四组方向分别做唯一编号传播测试 |
| 边界公式使用物理槽位错误 | 壁面滑移或密度发散 | 边界统一使用ET逻辑访问层 |
| 新生流体写错地址 | 颗粒尾部产生NaN | 局部重建后统一编码写回 |
| OpenMP写入冲突 | 结果随机变化 | 地址级证明和1线程/多线程对照 |
| 指针交换后错误释放 | 内存泄漏或重复释放 | 使用独立 `DistMemory` 保存原始地址 |

## 19. 与AA方案的选择建议

| 项目 | AA Pattern | Esoteric Twist |
| --- | --- | --- |
| 主数据布局 | 可以继续AoS | 推荐改成SoA |
| 时间步实现 | 奇偶两种内核 | 同一内核加对向指针交换 |
| 邻居访问 | 每隔一步访问全部方向邻居 | 每步访问选定半组邻居 |
| 对当前代码改动量 | 较小 | 较大 |
| 边界适配难度 | 高 | 更高 |
| 后续GPU/高带宽优化潜力 | 高 | 高 |

对于当前包含移动椭圆、新生流体和三种插值边界的项目，建议先实现AA方案并建立单数组验证体系，再实施ET。ET不应与AA同时改造，否则出现数值错误时很难判断是物理边界问题还是存储映射问题。

## 20. 预期收益

ET将分布函数主体内存从：

```text
2 × DX × DY × DQ × sizeof(double)
```

降低为：

```text
1 × DX × DY × DQ × sizeof(double)
```

同时消除独立全场传播数组。ET论文指出它在原位传播和间接寻址中可以减少内存占用与部分邻居索引读取；但本项目使用规则直角网格和CPU/OpenMP，最终性能必须在目标主机上实测。
