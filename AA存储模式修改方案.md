# AA Pattern 存储模式修改方案

## 1. 文档目标

本文只描述把当前双数组传播方式改造成 AA Pattern 单数组原位传播的实施方案，不直接修改源码。

改造目标：

1. 删除长期使用的第二套分布函数数组 `Temp`。
2. 使用一个 `Dist` 数组完成碰撞和传播。
3. 保持 `main()` 中现有时间步调用顺序和函数骨架。
4. 保持D2Q9方向编号、椭圆颗粒结构体数组和流固耦合公式不变。
5. 让固定壁面、压力边界、移动颗粒边界、宏观量和新生流体节点都能识别AA存储状态。

## 2. 参考依据

AA Pattern是一种线程安全的原位传播方法。它只保存一套分布函数，并在两种存储状态之间交替：

- 一个时间步只在本格点内读写，但把碰撞后分布写入反方向槽位。
- 另一个时间步从邻居的反方向槽位读取，并把碰撞结果写回对应的传播位置。

主要参考：

1. Wittmann, Zeiser, Hager, Wellein, *Comparison of different propagation steps for lattice Boltzmann methods*, 2013, DOI: [10.1016/j.camwa.2012.05.002](https://doi.org/10.1016/j.camwa.2012.05.002)。
2. Lehmann, *Esoteric Pull and Esoteric Push: Two Simple In-Place Streaming Schemes for the Lattice Boltzmann Method on GPUs*, 2022, DOI: [10.3390/computation10060092](https://doi.org/10.3390/computation10060092)。其附录给出了AA Pattern的奇偶步读写形式。

## 3. 当前程序的数据方式

当前程序使用：

```cpp
double (*Dist)[DY][DQ];
double (*Temp)[DY][DQ];
```

`GlobalCollide()` 的工作方式是：

```text
从 Dist[i][j][f] 读取
        ↓
完成BGK碰撞
        ↓
把结果传播到 Temp[ii][jj][f]
        ↓
交换 Dist 和 Temp 指针
```

优点是所有函数都可以认为 `Dist[i][j][f]` 永远表示格点 `(i,j)` 上物理方向 `f` 的分布函数。缺点是需要两套同样大的数组。

## 4. AA Pattern 的核心状态

为了让方案容易验证，定义两种明确状态。

### 4.1 状态 A：正常状态

```text
Dist[i][j][f]
```

直接表示格点 `(i,j)` 上沿方向 `f` 运动的分布函数。

### 4.2 状态 B：扭转状态

逻辑上属于格点 `(i,j)`、方向 `f` 的分布函数，实际存放在：

```text
Dist[i-Ex[f]][j-Ey[f]][Re[f]]
```

也就是说，传播和反方向换槽被编码进了存储位置。

程序每执行一个时间步就在状态A和状态B之间切换。

## 5. 两种时间步的读写规则

### 5.1 A到B：本地碰撞并写入反方向

对每个流体格点：

```cpp
for( int f=0; f<DQ; ++f )
    F[f] = Dist[i][j][f];

// 使用F计算密度、速度和平衡分布，再完成碰撞；

for( int f=0; f<DQ; ++f )
    Dist[i][j][Re[f]] = Post[f];
```

这一时间步只访问本格点，不访问邻居。

碰撞后沿方向 `f` 的结果写入本格点的反方向槽位。结合状态B的解释，这相当于完成了一次传播。

### 5.2 B到A：从邻居读取并恢复正常状态

对每个流体格点：

```cpp
for( int f=0; f<DQ; ++f )
{
    int ii = i-Ex[f];
    int jj = j-Ey[f];
    F[f] = Dist[ii][jj][Re[f]];
}

// 使用F完成碰撞；

for( int f=0; f<DQ; ++f )
{
    int ii = i+Ex[f];
    int jj = j+Ey[f];
    Dist[ii][jj][f] = Post[f];
}
```

完成后重新回到状态A。

每个线程写入的位置是唯一的，因此不需要额外的 `Temp` 数组。

## 6. D2Q9方向保持不变

当前方向编号继续使用：

| `f` | 方向 | `Ex[f]` | `Ey[f]` | `Re[f]` |
| ---: | --- | ---: | ---: | ---: |
| 0 | 静止 | 0 | 0 | 0 |
| 1 | 东 | 1 | 0 | 3 |
| 2 | 北 | 0 | 1 | 4 |
| 3 | 西 | -1 | 0 | 1 |
| 4 | 南 | 0 | -1 | 2 |
| 5 | 东北 | 1 | 1 | 7 |
| 6 | 西北 | -1 | 1 | 8 |
| 7 | 西南 | -1 | -1 | 5 |
| 8 | 东南 | 1 | -1 | 6 |

AA实现必须始终通过 `Re[f]` 获取反方向，不应假设对向方向编号相邻。

## 7. 建议新增的状态和辅助函数

在保持C风格的前提下，建议增加：

```cpp
const int AA_NORMAL = 0;
const int AA_TWISTED = 1;
int AaState;
```

建议增加以下普通函数，不引入类或模板封装：

```cpp
void LoadNodeAA( int i, int j, double F[DQ] );
void StoreNodeAA( int i, int j, const double F[DQ] );
double GetDistAA( int i, int j, int f );
void SetDistAA( int i, int j, int f, double Value );
void GlobalCollideAA();
```

职责：

- `LoadNodeAA()`：按照 `AaState` 把一个格点的9个物理分布读取到局部数组。
- `StoreNodeAA()`：按照当前阶段把碰撞后分布写入正确物理地址。
- `GetDistAA()`：让边界和受力函数读取一个逻辑分布函数，而不关心实际存储位置。
- `SetDistAA()`：让边界和新生流体重建写入逻辑分布函数。
- `GlobalCollideAA()`：替换原 `GlobalCollide()` 的双数组实现。

## 8. 各文件修改规划

### 8.1 `common.h`

计划增加：

- AA状态常量。
- 必要时增加周期下标辅助函数。

保持不变：

- D2Q9方向数组。
- `GridData`、`ParticleData` 和 `LinkData`。
- `Feq()` 公式。

### 8.2 `Migration.cpp`

计划修改：

1. 初始化时只分配 `Dist`。
2. 删除正式版本中的 `Temp` 分配、清零、交换和释放。
3. 初始化结束时令 `AaState=AA_NORMAL`。
4. 将 `GlobalCollide()` 内部改成AA奇偶状态分支，但保留函数名和主循环位置。
5. `MacroCalculate()` 改用 `LoadNodeAA()` 获取逻辑分布。
6. `BoundaryCollide()` 和 `InletAndOutlet()` 改用 `GetDistAA()/SetDistAA()`。
7. `LatticeRemesh()` 只移动一套分布数据，并保持 `AaState` 不变。

### 8.3 `Particle.cpp`

计划修改：

1. `ParticleCollide()` 的所有 `Dist` 直接访问替换为AA逻辑访问函数。
2. `ParticleForce()` 通过AA逻辑访问函数读取入射和反射分布。
3. `ParticleGeometry()` 重建新生流体时先在局部 `double F[DQ]` 中计算，再通过 `SetDistAA()` 写入。
4. 颗粒内部节点清零时，按当前 `AaState` 清除对应的逻辑分布。
5. 颗粒结构体数组、连线颗粒编号和椭圆几何公式保持不变。

## 9. 主循环骨架

主循环顺序不改变：

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

改变的是这些函数访问 `Dist` 的方法。

在每个时间步结束时，所有函数必须看到同一个 `AaState`。禁止某个函数私自切换状态。

## 10. 边界条件改造重点

### 10.1 周期下标

AA的B到A步骤会访问邻居，应沿用当前x方向周期回绕规则。y方向不能简单周期回绕，因为上下边界是固定壁面。

### 10.2 上下固定壁面

原 `BoundaryCollide()` 中的HBB、BFL和LL2公式可以保留，但参与公式的每个分布值都必须通过AA逻辑访问函数获取。

不允许直接假设：

```cpp
Dist[i][j][f]
```

就是当前物理方向 `f`。

### 10.3 压力入口出口

入口出口公式本身不变。需要把左端和右端的9个逻辑分布读入局部数组，重建未知方向后再按当前AA状态写回。

### 10.4 移动颗粒边界

`LinkData` 中的 `i、j、f、d、p` 都可以保留。

需要保证：

- `f` 始终表示物理方向，而不是数组当前槽位。
- `r=Re[f]` 始终表示物理反方向。
- `ParticleCollide()` 的读写经过AA地址转换。
- `ParticleForce()` 读取的是完成边界处理后的逻辑分布。

## 11. 新生流体格点

这是AA改造中风险最高的部分之一。

当前代码同时写：

```cpp
Dist[i][j][f]
Temp[i][j][f]
```

AA模式只剩一个数组，而且逻辑上属于同一格点的9个方向在扭转状态下可能分散在邻居槽位中。

建议流程：

1. 在局部数组 `F[DQ]` 中完成外推平均。
2. 由 `F` 计算新生节点密度和速度。
3. 使用 `SetDistAA(i,j,f,F[f])` 逐方向写入当前AA状态对应的位置。
4. 完成全部新生节点写入后，再把 `MARKED` 改回 `FLUID`。

不能直接对 `Dist[i][j]` 连续写9个值。

## 12. 初始化、输出和重启

### 初始化

初始分布采用状态A：

```text
AaState = AA_NORMAL
Dist[i][j][f] = Feq(...)
```

### 输出

`SaveData()` 当前只输出颗粒数据，不读取完整分布函数，因此基本不受影响。

如果以后输出完整流场，必须通过 `LoadNodeAA()` 解码。

### 重启文件

如果以后保存重启数据，必须额外保存 `AaState`。只保存原始 `Dist` 而不保存状态，会导致恢复后所有方向解释错误。

## 13. 实施顺序

建议分七个提交实施：

1. **增加AA逻辑访问层**：保留 `Temp`，只验证访问器。
2. **实现A到B本地时间步**。
3. **实现B到A邻居时间步**。
4. **适配固定壁面和压力边界**。
5. **适配颗粒碰撞、受力和新生流体节点**。
6. **适配宏观量、重映射和输出**。
7. **结果一致后删除 `Temp`**。

不要在第一个提交中直接删除 `Temp`，否则出现数值差异时难以定位。

## 14. 验证方案

### 14.1 无颗粒通道流

暂时关闭 `MIGRATION`，比较原模式和AA模式：

- 密度场。
- 速度场。
- 总质量。
- 泊肃叶解析解误差。

### 14.2 固定单颗粒

令颗粒速度和角速度保持0，比较：

- `LinkNum`。
- 每条连线的 `d`。
- 颗粒表面附近密度和速度。
- 颗粒受力。

### 14.3 运动单颗粒

比较颗粒的：

- `Sx、Sy`。
- `Vx、Vy`。
- `ParFai、Va`。
- `Fx、Fy、ParTor`。

### 14.4 多颗粒

最后比较三个颗粒的完整输出轨迹。

建议先运行偶数步，再运行奇数步，确认两种AA状态都经过验证。

## 15. 验收标准

1. 正式计算只分配一套D2Q9分布函数。
2. OpenMP运行时不存在同一存储位置的并发写入。
3. 奇偶步的质量和速度没有周期性跳变。
4. 边界格式HBB、BFL、LL2都能单独运行。
5. 新生流体节点不会产生NaN或 `Num==0` 问题。
6. 与原双数组版本在允许浮点误差内一致。
7. 多颗粒编号、受力和轨迹不会串位。

## 16. 主要风险

| 风险 | 后果 | 控制方法 |
| --- | --- | --- |
| 奇偶状态解释错误 | 速度方向颠倒或结果每步振荡 | 为两种状态分别编写小规模地址测试 |
| 边界仍直接访问数组 | 壁面和颗粒附近数值错误 | 所有边界统一经过逻辑访问层 |
| 新生流体写错物理位置 | 颗粒尾部出现NaN | 局部数组重建后统一编码写回 |
| OpenMP地址冲突 | 结果不稳定且不可重复 | 审核每个方向的唯一写入关系 |
| 输出或重启忽略状态 | 文件解读错误 | 输出前解码，重启时保存 `AaState` |

## 17. 预期收益

当前分布函数内存约为：

```text
2 × DX × DY × DQ × sizeof(double)
```

AA正式版本约为：

```text
1 × DX × DY × DQ × sizeof(double)
```

分布函数主体内存可降低约一半，同时省去每步交换两套数据所有权的双数组需求。实际速度收益需要在目标主机上测量，不能只根据内存数量推断。
