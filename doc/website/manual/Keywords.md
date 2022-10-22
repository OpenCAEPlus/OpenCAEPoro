---
head:
  - - link
    - rel: stylesheet
      href: https://cdnjs.cloudflare.com/ajax/libs/KaTeX/0.5.1/katex.min.css
title: Introduction to OpenCAEPoro Keywords
---

# Introduction to OpenCAEPoro Keywords

OpenCAEPoro 目前用一套关键字体系来构建模拟选项，控制模拟过程，大部分关键字格式与内容与斯伦贝谢的 Eclipse 保持一致，对于这些关键字，将在其后用 **(e)** 或 **(e300)** 标出，**(e300)** 表示仅和适用于 E300 的选项匹配。有些关键字后面还会用 **(/)** 标出，这表示必须在其内容结束的**下一行**使用 `/` 表示其已结束，例如

```text
PVTW (e)(/)
```

表示关键字`PVTW`和E100兼容，需用`/`标记输入结束，具体使用方式如下：

```text
PVTW
4014.7      1.0     3E-6       0.3100    0.0
/
```

文本中的 `'` 和 `,` 将会自动被替换成空格。**(TODO)** 标识了目前尚不完善的地方。

注意：目前，从关键字读入的数据只适用于**英制**单位。

## 注释

在 OpenCAEPoro 中，有以下几种注释方式

* 以 `#` 或 `--` 开头的行，例如

```text
  # comments
     # comments
  -- comments 
     -- comments
```

* 在任何 `/` 的背后，例如

```text
  DIMENS
  6  6  6  / comments
```

注意：不可在关键字后面注释 ([EQUALS](#_EQUALS) 下的子关键字除外)，否则该关键字将被忽略，例如

```text
DIMENS  /  comments
```

## BLACKOIL<span id=_BLACKOIL></span> (e)

BLACKOIL 用来定义黑油模型，示例：

```text
BLACKOIL
```

## OIL<span id=_OIL></span> (e)

OIL 表明在模拟过程中油相可能存在，用于黑油模型，在组分模型中将会被忽略。示例：

```text
OIL
```

## GAS<span id=_GAS></span> (e)

GAS 表明在模拟过程中气相可能存在，用于黑油模型，，在组分模型中将会被忽略。示例：

```text
GAS
```

## WATER<span id=_WATER></span> (e)

WATER 表明在模拟过程中水相可能存在，注意：当前的组分模型中，水相一定存在。示例：

```text
WATER
```

## DISGAS<span id=_DISGAS></span> (e)

DISGAS 表明在模拟过程中溶解气可能存在于油相中，用于黑油模型，示例：

```text
DISGAS
```

## COMPS<span id=_COMPS></span>

COMPS 激活了组分模型，并输入烃组分数，示例：

```text
COMPS
6
```

## DIMENS<span id=_DIMENS></span> (e)(/)

DIMENS 用来定义储层网格的规模，适用于三维结构网格，用三个数字分别给出了 $x,y,z$ 方向上网格的维度 $N_{x},N_{y},N_{z}$，示例：

```text
DIMENS
3  4  5
```

## TABDIMS<span id=_TABDIMS></span>

TABDIMS 定义了饱和度表格和PVT表格的最大数量，不同的饱和度表格和PVT表格可能用在不同的网格区域，这将通过 [SATNUM](#_SATNUM) 和 [PVTNUM](#_PVTNUM) 来定义。示例：

```text
TABDIMS
-- NTSFUN  NTPVT
      1      2
```

表明储层中至少有一个饱和度计算区域和两个 PVT 计算区域。

## EQUALS<span id=_EQUALS></span> (e)(/)

EQUALS 用来批量地为某些油藏性质赋值，其中的子关键字具有如下的一般格式。示例：

```text
item  var  i1  i2  j1  j2  k1  k2
```

其中 item 表示要赋值的性质名称，var 表示要赋的值，后面六个数表示赋值的网格范围，依次是 $x$ 坐标的范围 $[i_{1},i_{2}]$，$y$ 坐标的范围 $[j_{1},j_{2}]$ 和 $z$ 坐标的范围 $[k_{1},k_{2}]$，例如：

```text
DX  1000  1  3  1  3  1  3
```

表示油藏区域中坐标范围为 $[1,3]\times[1,3]\times[1,3]$ 的网格块的 $x$ 方向的长度为 1000 (英尺)。坐标范围的值还支持简写，假设三维的坐标范围分别为 $[1,N_{x}]\times[1,N_{y}]\times[1,N_{z}]$，则

```text
Dx  1000  6*              / Dx = 1000 in all Region
Dx  1000  1  3  2*  2  4  / DX = 1000 in [1,3]x[1,Ny]x[2,4]
Dx  1000  4* 2  6         / DX = 1000 in [1,Nx]x[1,Ny]x[2,6]
```

目前，在 EQUALS 中，支持如下子关键字的设置

```text
DX        / x 轴方向网格长度，feet
DY        / y 轴方向网格长度，feet
DZ        / z 轴方向网格长度，feet
PORO      / 网格块的孔隙度，dimensionless
NTG       / 网格块的净毛率，dimensionless
PERMX     / 网格块 x 轴方向的绝对渗透率，mD
PERMY     / 网格块 y 轴方向的绝对渗透率，mD
PERMZ     / 网格块 z 轴方向的绝对渗透率，mD
SATNUM    / 网格块所使用的的饱和度表格编号
PVTNUM    / 网格块所使用的的PVT表格编号
TOPS      / 顶层网格块的上表面深度（顶层的 z 坐标为1，因此后两个数字为1），feet
```

## COPY<span id=_COPY></span> (e)(/)

COPY 用于对 EQUALS 中子关键字的内容进行相互拷贝，其坐标格式与 [EQUALS](#_EQUALS) 中一致，例如：

```text
COPY
PERMX  PERMY   1  3  1  3  1  3
```

表示将坐标范围 $[1,3]\times[1,3]\times[1,3]$ 中的网格块的 PERMY 属性的值设为 PERMX 属性的值。

## MULTIPLY<span id=_MULTIPLY></span> (e)(/)

MULTIPLY关键字用于对 EQUALS 中子关键字的内容进行放缩，其坐标格式与 [EQUALS](#_EQUALS) 中一致，例如

```text
MULTIPLY
PERMZ  0.1  1  3  1  3  1  3
```

表示将区域 $[1,3]\times[1,3]\times[1,3]$ 中的网格块的 PERMZ 属性的值乘以 0.1

## DX DY DZ PORO NTG PERMX PERMY PERMZ TOPS (e)(/)

关键字用于定义每个网格块的对应属性的值，见 [EQUALS](#_EQUALS)，网格块依照 x -> y -> z 的坐标字典序。

```text
DX
10 20 30 24*20
/
```

## SWATINIT<span id=_SWATINIT></span> (e)(/)

SWATINIT 用于定义初始水相饱和度，并以此缩放毛管力曲线使得水相饱和度与初始平衡解一致

```text
SWATINIT
0.2  0.22  0.25  24*0.3
/
```

## SATNUM<span id=_SATNUM></span> (e)(/)

SATNUM 用来定义饱和度表格作用区域，即指定网格块使用相应的饱和度表格，饱和度表格的编号按照其输入顺序依次排序，参考 [SWOF](#_SWOF)

```text
SATNUM
2*1   2   24*3
/
```

## ACTNUM<span id=_ACTNUM></span> (e)(/)

ACTNUM 用来定义网格块的状态，目前，0 表示死网格，1表示活网格

```text
ACTNUM
2*1  0   24*1
/
```



## PVTNUM<span id=_PVTNUM></span> (e)(/)
PVTNUM 用来定义PVT区域，即指定网格块使用相应的PVT关系，与 [SATNUM](#_SATNUM) 类似

```text
PVTNUM
2*1   2   24*3
/
```

## COORD<span id=_COORD></span> (e)(/)

COORD 用于角点网格的定义，它给出了角点网格的“骨架”，如果网格规模为 $N_{x} * N_{y} * N_{z}$，那么在 COORD 中就需要给出 $(N_{x} + 1)(N_{y} + 1)$ 条垂直方向的线，按照 $x \rightarrow y$ 的字典序，而每条线由两个三维点表示，顶层的点在前，底层的点在后，往下为正，单位英尺，例如，对于 $N_{x} = 3，N_{y} = 2$

```text
COORD  
0      0     1000    0      0     2000
1000   0     1000    1000   0     2000
2000   0     1000    2000   0     2000
3000   0     1000    3000   0     2000
0      2000  1000    0      2000  1000
1000   2000  1000    1000   2000  1000
2000   2000  1000    2000   2000  1000
3000   2000  1000    3000   2000  1000
0      4000  1000    0      4000  1000
1000   4000  1000    1000   4000  1000
2000   4000  1000    2000   4000  1000
3000   4000  1000    3000   4000  1000
/
```

## ZCORN<span id=_ZCORN></span> (e)(/)

ZCORN 用于角点网格的定义，它通过给出每个网格的8个角点在对应“骨架”(见 [COORD](#_COORD) )上的深度。

```text
 3 -- 4
/|   /|
1 -- 2|
||   ||
| 7-- 8
|/   /
5 -- 6
```

给出的顺序为：

* 网格层： (k=1 -> k=nZ)
* 每一行网格(k固定，j -> nY)：给出各网格块节点 1 和 2 的深度，先给出 j = 1 行，$v_{1,1},v_{2,1},...,v_{1,nX},v_{2,nX}$，其后是 $v_{3,1},v_{4,1},...,v_{3,nX},v_{4,nX}$，再依次给出后面的行
* 上一点结束后，同样的方式给出 $v_{5,1},v_{6,1},...,v_{5,nX},v_{6,nX}$，其后是 $v_{7,1},v_{8,1},...,v_{7,nX},v_{8,nX}$

## RTEMP<span id=_RTEMP></span> (e)

RTEMP 定义恒温油藏的温度，单位为 °F，示例

```text
RTEMP
200
```

## SWOF<span id=_SWOF></span> (e)(/)

SWOF 定义当油相和水相存在时的饱和度表格。如果有多个表格，则需用 `/` 隔开。在SWOF中一共由四列数据，水相的饱和度作为自变量，这四列数据分别表示

```text
SWOF
-- Sw     Krw     Krow     Pcow
 .2200   .0000   1.0000   7.0000
 .3000   .0700   0.4000   4.0000
 .4000   .1500   0.1250   3.0000
 .5000   .2400   0.0649   2.5000
 .6000   .3300   0.0048   2.0000
 .8000   .6500   0.0      1.0000
 .9000   .8300   0.0      .5000
 1.0000  1.0000  0.0      .0000  / table 1
 .18     .00     1.00     0
 .32     .07     0.38     0
 .50     .31     0.05     0
 .60     .38     0.004    0
 .80     .57     0.0      0
 1.00    1.00    0.0      0      / table 2
 /
```

每一列的具体含义如下：

* Sw : 水相的饱和度，dimensionless
* Krw : 水相的相对渗透率，dimensionless
* Krow : 当仅有油相和水相存在时，油相的相对渗透率，dimensionless
* Pcow : 水相的毛管力 Po - Pw，psia

## SWFN<span id=_SWFN></span> (e)(/)

SWFN 定义水相的饱和度表格。如果有多个表格，则需用 `/` 隔开，与 [SWOF](#_SWOF) 类似。在 SWFN 中一共有三列数据，水相的饱和度作为自变量，这三列数据分别表示:

* Sw : 水相的饱和度，dimensionless
* Krw : 对应的水相相对渗透率，dimensionless
* Pcow : 对应的水相的毛管力 Po - Pw，psia

## SGOF<span id=_SGOF></span> (e)(/)

SGOF 定义当气相和水相同时存在时的饱和度表格。如果有多个表格，则需用 `/` 隔开，与 [SWOF](#_SWOF) 类似。在SGOF中一共由四列数据，气相的饱和度作为自变量，这四列数据分别表示：

* Sg : 水的饱和度，dimensionless
* Krg : 对应的气相的相对渗透率，dimensionless
* Krog : 当油相，气相和束缚水存在时，对应的油相相对渗透率，dimensionless
* Pcog : 对应的气相毛管力 Pg - Po，psia

## SGFN<span id=_SGFN></span> (e)(/)

SGFN 定义气相的饱和度表格。如果有多个表格，则需用 / 隔开，与 [SWOF](#_SWOF) 类似。在 SGFN 中一共有三列数据，气相的饱和度作为自变量，这三列数据分别表示：

* Sg : 气相的饱和度，dimensionless
* Krg : 对应的气相相对渗透率，dimensionless
* Pcog : 对应的气相毛管力 Pg - Po，psia

## SOF3<span id=_SOF3></span> (e)(/)

SOF3 用来定义油相的饱和度表格(用于三相模型)。如果有多个表格，则需用 / 隔开，与 [SWOF](#_SWOF) 类似。在 SGFN 中一共有三列数据，油相的饱和度作为自变量，这三列数据分别表示
* So : 油相的饱和度，dimensionless
* Krog : 当仅有油相和水相存在时，对应的油相的相对渗透率，dimensionless
* Krow : 当油相，气相和束缚水存在时，对应的油相相对渗透率，dimensionless

## PVCO<span id=_PVCO></span> (e)(/)

PVCO 适用于黑油模型，它通过表格数据提供了活油（含有溶解气）的 PVT 性质。如果有多个表格，则需用 `/` 隔开，与 [SWOF](#_SWOF) 类似。PVCO一共有六列数据，其中，前两列作为自变量(泡点压力和气油比)，即其余的量随着它们之一的改变而改变。这六列分别表示：

* Pbub : 油相的泡点压力，psia
* Rs : 油相的气油比，Mscf/stb，它等于标准状态下从油相中析出的气体的体积与油相的体积之比
* Bo : 对应的 Pbub 压力下油相的体积系数，rb / stb，它等于相同摩尔数的油相在储层状态下与在标准状态下的体积之比
* Viscosity : 对应的 Pbub 压力下油相的粘性，cP
* Cb : 对应的 Rs 下不饱和油相的压缩性，1/psia，它等于 Bo 对于压力的相对变化率
* Cv : 对应的 Rs 下不饱和油相的粘性的变化率，1/psia，它等于粘性对于压力的相对变化率

## PVDO<span id=_PVDO></span> (e)(/)

PVDO 适用于黑油模型，通过表格数据提供了死油 (不含有溶解气) 的PVT 性质。如果有多个表格，则需用 `/` 隔开，与 [SWOF](#_SWOF) 类似。PVDO 一共有三列数据，其中油相压力为自变量，这三列分别表示：

* Po : 油相压力，psia
* Bo : 对应的油相体积系数，rb/stb
* Viscosity : 对应的油相粘性，cP

## PVDG<span id=_PVDG></span> (e)(/)

PVDG 适用于黑油模型，通过表格数据提供了干燥气的 PVT 性质。如果有多个表格，则需用 `/` 隔开，与 [SWOF](#_SWOF) 类似。PVCO 一共有三列数据，气相的压力作为自变量，这三列数据分别表示：

* Pg : 气相压力，psia
* Bg : 对应的气相体积系数，rb/Mscf
* Viscosity : 对应的气相粘度，cP

## PVTW<span id=_PVTW></span> (e)(/)

PVTW 用于黑油模型或者单独考虑水相的组分模型，通过表格数据提供了水相的 PVT 性质。如果有多个表格，则需用 `/` 隔开，与 [SWOF](#_SWOF) 类似。PVTW 一共有五列数据，水相的压力作为自变量，这五列数据分别表示：

* Pw : 水相的压力，psia
* Bw : 对应的水相体积系数，rb/stb
* Cw : 对应的水相压缩性系数，1/psia，它等于体积系数关于压力的相对变化率
* Viscosity : 对应的水相粘度，cP
* Cv : 对应的水相的粘度关于压力的相对变化率，1/psia

## PBVD<span id=_PBVD></span> (e)(/)

PBVD 通过表格数据给出了初始油藏下，泡点压力关于位置深度的关系。如果有多个表格，则需用 `/` 隔开，与 [SWOF](#_SWOF) 类似。PBVD 一共有两列数据，分别是：

* Depth : 采样点深度，feet
* Pbub : 对应的泡点压力，psia

## MISCIBLE<span id=_MISCIBLE></span> (e300)

在组分模型中，MISCIBLE 关键字开启混溶选项。

## MISCSTR<span id=_MISCSTR></span> (e)

在混溶模型中，对界面的表面张量进行设置。必须先打开 [MISCIBLE](#_MISCIBLE) 选项。需要输入三个值

* 最大混溶表面张量：当表面张量大于次值时，油气不混溶， dynes/cm
* 所期待的最大表面张量：dynes/cm (暂时没有用上)
* 用来放缩输入毛管力曲线的最大表面张量：dynes/cm

## EQUIL<span id=_EQUIL></span> 

EQUIL 用于油藏初始条件的计算，它给出了参考深度处的油藏压力，油气接触面的深度及毛管力，油水接触面的深度及毛管力。它的六列依次为

* Dref : 参考深度，feet
* Pref : 参考深度处的压力，psia
* D1 : 某接触面深度，feet
* P1 : D1 深度的两相毛管力，psia
* D2 : 某接触面深度，feet
* P2 : D2 深度的两相毛管力，psia

根据存在的相的情况依次给出需要的值

* 如果存在油相和水相，则
  - D1 为油水接触面深度，P1 = Po - Pw
  - 如果气相也存在，则 D2 为油气接触面深度，P2 = Pg - Po
* 如果只存在气相和水相，则
  - D1 为气水接触面深度，P1 = Pg - Pw
  - D2 和 P2 无需给出

EQUIL 应与上述的饱和度表格与 PVT 表格一样，在不同的区域有不同的数值，但由于多区域功能目前并未完全实现，因此只读入一行数据。**(TODO)**

## ZMFVD<span id=_ZMFVD></span> (e)(/)

ZMFVD 用于组分模型，输入油藏中烃组分摩尔占比 (dimensionless) 关于深度 (feet) 的变化关系表，输入格式具体可参见 [SWOF](#_SWOF)。如下是一个六种烃组分的示例：

```text
ZMFVD
1000.0   0.5  0.03  0.07  0.2  0.15  0.05
10000.0  0.5  0.03  0.07  0.2  0.15  0.05 
/
```

## ROCK<span id=_ROCK></span> 

ROCK 给出了岩石的可压缩性，一共有两列数据，分别表示

* Pref : 参考压力，psia
* C : 岩石的可压缩性，1/psia，它是岩石体积关于压力的相对变化率。

ROCK 应与上述的饱和度表格与PVT表格一样，在不同的区域有不同的数值，但由于多区域功能目前并未完全实现，因此只读入一行数据。**(TODO)**

## GRAVITY<span id=_GRAVITY></span> 

GRAVITY 用于黑油模型或单独考虑水相的组分模型，用于计算标准状态下流体的参考密度，分别需要键入油，水，气三相的重力因子

* oil：默认值 45.5
* water (with reference to pure water)：默认值 1.0
* gas (with reference to air)：默认值 0.7773

对于组分模型，只需给出水相的信息即可。

GRAVITY 应与上述的饱和度表格与PVT表格一样，在不同的区域有不同的数值，但由于多区域功能目前并未完全实现，因此只读入一行数据。**(TODO)**

## DENSITY<span id=_DENSITY></span>

DENSITY 用于黑油模型或单独考虑水相的组分模型，给出了标准状态下流体的参考密度，分别需要键入油，水，气三相的密度

* oil：默认值 37.457 lb/ft3
* water：默认值 62.366416 lb/ft3
* gas：默认值 0.062428 lb/ft3

对于组分模型，只需给出水相的信息即可

DENSITY 应与上述的饱和度表格与PVT表格一样，在不同的区域有不同的数值，但由于多区域功能目前并未完全实现，因此只读入一行数据。**(TODO)**

## INCLUDE<span id=_INCLUDE></span> (e)

INCLUDE 关键字用于分文件编写输入文件，需要输入被包含文件的  **相对路径**，被包含文件的格式与主输入文件是一致的，示例：

```text
INCLUDE
./SomeFile.inc
```

## METHOD<span id=_METHOD></span>

METHOD 关键字用来确定所使用的的离散方法以及所调用的线性求解器，离散方法包括 IMPEC (隐式压力显式组分)，FIM (全隐式类方法)。线性求解器默认使用FASP (如需使用别的求解器，需要在程序中补充对应的接口)；因此，对应地需要给出求解常量矩阵(IMPEC)或块状矩阵(FIM)的FASP输入文件 (**相对路径**)，例如：

```text
METHOD
FIM ./bsr.fasp
```

模拟器默认使用 IMPEC 方法，而FASP输入文件的优先级为

1. Method 关键字的输入文件中给定
2. (当前目录下) ./csr.fasp (IMPEC)，./bsr.fasp (FIM)
3. (当前目录下) ../conf/csr.fasp (IMPEC)， ../conf/bsr.fasp (FIM)，
4. 内置参数

## TSTEP<span id=_TSTEP></span> (e)(/)

TSTEP 关键字通过时间间隔给出了模拟的关键时间节点 (day)，这是在模拟中会强制到达的时间点。在这些时间节点上，井的控制方式可能会发生改变，油藏状态可能会进入不同的阶段 (根据经验预估)，由此求解参数可能会做出调整。第 0 天始终为第 1 个时间节点，因此无需再输入。TSTEP 关键字支持简写，例如 2*100 表示两个 100 天的间隔。于是下面的关键时间节点为第 0, 5, 15, 45, 145, 245 天。

```text
TSTEP
5  10  30  2*100
/
```

在输入文件中一般包含多个 TSTEP 关键字。对同一对象的控制作用范围为这两个控制关键字中间的时间段，如下，第一个 [TUNING](#_TUNING) 的作用范围为中间的 115 天

```text
TUNING
****** (contents)
/

TSTEP
5 10
/

TSTEP
100
/

TUNING
****** (contents)
/
```

此外，如果激活了 [RPTSCHED](#_RPTSCHED) 关键字，则会在每一个时间节点打印关键字里对应的信息。

## TUNING<span id=_TUNING></span> (/)

TUNING 关键字给出了模拟求解的参数，包括时间步长控制参数和非线性求解控制参数，它是与**时间相关**的控制，见 [TSTEP](#_TSTEP)。TUNING 关键字的内容分为三部分，每一部分内容结束后需用 `/` 结尾，前两部分与时间步长的控制相关，最后一部分与非线性求解相关，主要用于 FIM 类方法。示例：

```text
TUNING
-- Init   max    min    incre   chop    cut
   1      10     0.1    3       0.15    0.3               /
-- dPlim  dSlim   dNlim   dVerrlim
   300    0.5     0.3     0.001                           /
-- itNRmax  NRtol  dPmax  dSmax  dPmin   dSmin   dVerrmax
   10       1E-2   2000   2      1       0.01    0.01     /
/
```

* 第一部分以时间步长控制相关，最重要的参数为前三个
  * Init：每一个关键时间节点开始的最大时间步长，day
  * max：相邻两个关键时间节点中，能取到的最大时间步长，day
  * min：相邻两个关键时间节点中，下一个时间步长的预测下限。在自适应时间步长选取下，这给出了下一个时间步预测的最小值，但有时为了保证解的稳定性和合理性，时间步长会再缩短，因此真实的时间步长是可能会低于预测下限，day
  * incre：在自适应时间步长选取下，下一个时间步长相对于当前时间步长的最大增长因子，dimensionless
  * chop：在自适应时间步长选取下，下一个时间步长相对于当前时间步长的最小缩减因子，dimensionless
  * cut：在自适应时间步长选取下，在 FIM 方法中，当 Newton 方法失败时时间步长的缩短因子，dimensionless
* 第二部分以下一时间步长的预测为主，它是假定变量基于时间线性变化的线性预测，它根据以下变量进行预测，并不保证下一时间步的真实情况也满足条件
  * dPlim：下一时间步每个网格块能接受的最大压力变化，psia
  * dSlim：下一时间步每个网格块每一相能接受的最大饱和度变化，dimensionless
  * dNlim：下一时间步每个网格块每种组分能接的最大摩尔数相对变化，dimensionless
  * dVerrlim：下一时间步每个网格块能接受的孔隙体积与流体体积的相对误差的最大值，dimensionless
* 第三部分与非线性迭代的控制相关，它主要用在 FIM 类方法求解非线性方程的 Newton 迭代中
  * itNRmax：每个时间步最大的 Newton 迭代步数，如果迭代超过该步数仍未收敛，则缩短时间步长重新迭代。
  * NRtol：Newton 迭代的收敛误差，目前在 OpenCAEPoro 中，这个值用于多种误差控制，包括残差的无穷范数、相对于网格块孔隙体积的残差的无穷范数、相对于网格块组分摩尔总数的残差的无穷范数
  * dPmax：每个 Newton 迭代最大的压力变化，psia
  * dSmax：每个 Newton 迭代步最大的饱和度变化，由于从组分摩尔数到相饱和度需要经历复杂的相平衡计算(组分模型)，因此基于程序性能的考量，所使用的的饱和度变化是一个预测值而非真实值
  * dPmin，dSmin：Newton 迭代最小的压力变化和饱和度变化。当一步 Newton 迭代的最大压力变化和最大饱和度变化 (均为真实值) 分别低于 dPmin 和 dSmin 时，则认为 Newton 迭代收敛
  * dVerrmax：每个网格块孔隙体积与流体体积的相对误差，此选项一般用于 IMPEC 类方法，因为 FIM 类方法对此具有较好的保证

## WELSPECS<span id=_WELSPECS></span> (/)

WELSPECS 用于输入井的信息，包括

* 井的名字
* 所属的井组名称 (暂时无用处)
* 井所在网格块的 x 轴坐标
* 井所在网格块的 y 轴坐标
* 井的深度，feet

```text
WELSPECS
INJE1   G   1   1     1*
PROD1   G   10  10    1*
/
```

在井的信息输入中，可能会使用默认值，如上面的 `1*`，表示对应的位置是默认值。如果是 `m*`，则表示对应的 m 个位置的值为默认值。

在 WELSPECS 中，只有井的深度允许为默认值，当它被设为**默认值**或者**负数**时，井的位置位于第 0 个射孔的位置，关于射孔的排序，见 [COMPDAT](#_COMPDAT)。

## COMPDAT<span id=_COMPDAT></span> (/)

COMPDAT 用于定义井的射孔。使用该关键字时，需要将所定义的射孔的信息匹配到对应的井上去，目前有精确匹配和模糊匹配两种匹配方式。

* 精确匹配

  需要准确的井的名字，如

  ```text
  COMPDAT
  INJ1    ... 
  /
  ```

  会将其后定义的信息准确匹配到井 `INJ1`。

* 模糊匹配

  在末尾带上 `*`，则会匹配到名称以前面部分开头的井，如

  ```text
  COMPDAT
  INJ*    ... 
  /
  ```

  则会匹配到所有名字以 `INJ` 开头的井。

在 COMPDAT 中，包括了

* 所需匹配的井的信息
* 射孔所在网格块的 x 轴坐标，默认则为所属井的 x 轴坐标。一般来说，水平井的射孔 $x,y$ 坐标与井并不相同
* 射孔所在网格块的 y 轴坐标，默认则为所属井的 y 轴坐标。一般来说，水平井的射孔 $x,y$ 坐标与井并不相同
* 射孔所在网格块的 z 轴坐标范围，由两个升序的正整数组成
* 射孔与网格块连接处的传播因子(WI)，rb/day/psia，默认或为负则由程序自行计算得到
* 射孔的直径，feet，默认为1  feet
* 射孔与网格块连接处的等效渗透因子(Kh)，mD*ft，默认或为负则由程序自行计算得到
* 射孔的 Skinfactor，默认则为 0
* 射孔的方向：x，y，z。默认为 z，即竖直方向

示例：

```text
COMPDAT
INJE*   2*   1   2     1*   0.5   3*   
PROD1   1   3   3   3     1*   0.5   3*    
/
```

对于名字中以 `INJE` 开头的井，它拥有与井横纵坐标相同的两个射孔，这两个射孔的纵坐标分别为 1 和 2，射孔的直径为 0.5 ft，其余皆为默认值。

对于井 PROD1，它拥有一个坐标位置在 (1,3,3) 的一个射孔，射孔直径为 0.5ft，其余皆为默认值。

注意，在 OpenCAEPoro 中，要求给出的射孔顺序是由高至低的 (这由于射孔间静水压差的计算方式决定)，目前并没有在内部依据射孔深度对射孔进行排序；同时如果井的深度是默认的，则井的深度则为位置最高的射孔的深度。**(TODO)**

## WCONINJE<span id=_WCONINJE></span> (/)

WCONINJE 用于对注入井的控制，这是与时间相关的控制，见 [TSTEP](#_TSTEP)。请注意，在OpenCAEPoro中，一口井是生产井还是注入井，这是一个动态的状态，而不是一开始就确定的结果，这意味着井的类型是可以随时切换的。因此，WCONINJE 不仅赋予了井的一般控制，也赋予了井的类型。

在 WCONINJE 中，同样需要匹配井的名字，方式同 [COMPDAT](#_COMPDAT)。 WCONINJE 中包含如下的控制

1. 所需匹配的井的信息
2. 注入流体的类型：包括水相 `WAT` 或 `WATER`，气相 `GAS`；在组分模型中，也可以是任意的流体，但需要在 [WELLSTRE](#_WELLSTRE) 中给出对应的组分比例
3. 井的状态：开 `OPEN`，关 `CLOSE`
4. 井的控制方式：注入气体的流速控制 `RATE`，注入井的恒压控制(井的参考位置处) `BHP`
5. 流速控制的值或流速上限：注气的单位为 Mscf/day，注液的单位为 stb/day
6. 压力控制的值或压力上限：单位为 psia

如果井的初始控制方式(关键字作用时的控制方式)是 `RATE`，则 5 为流速控制的值，6 为压力上限，当井的压力超过该上限时，则切换为 `BHP` 控制。如果井的当前控制方式是 `BHP`，且流速超过了 5 所指定的值，则切换回 `RATE` 控制。

如果井的控制方式是 `BHP`，且 5 是默认值，则始终保持为 `BHP` 控制。若 5 不是默认值，则为最大注入流速的限制

示例：

```text
WCONINJE
INJE1   GAS   OPEN   RATE   100000.0      10000
INJE2   GAS   OPEN   BHP    1*            10000
/
```

## WCONPROD<span id=_WCONPROD></span> (/)

WCONPROD  用于对生产井的控制，这是与**时间相关**的控制，见 [TSTEP](#_TSTEP)。请注意，在OpenCAEPoro中，一口井是生产井还是注入井，这是一个动态的状态，而不是一开始就确定的结果，这意味着井的类型是可以随时切换的。因此，WCONPROD 不仅赋予了井的一般控制，也赋予了井的类型。

在 WCONPROD 中，需要匹配井的名字，方式同 [COMPDAT](#_COMPDAT)。WCONPROD 中包含如下的控制

1. 所需匹配的井的信息
2. 井的状态：开 `OPEN`，关 `CLOSE`
3. 井的控制方式：产油速率控制 `ORAT`，产气速率控制`GRAT`，产水速率控制`WRAT`，井的恒压控制(井的参考位置处) `BHP`
4. 流速控制的值或流速上限：产出气的单位为 Mscf/day，产出水或油的单位为 stb/day
5. 压力控制的值或压力上限：单位为 psia

如果井的初始控制方式(关键字作用时的控制方式)是 `*RAT`，则 4 为流速控制的值，5 为压力下限，当井的压力低于该下限时，则切换为 `BHP` 控制。如果井的当前控制方式是 `BHP`，且当前流速超过 4 所指定的值，则切换回 `RATE` 控制。

如果井的初始控制方式是 `BHP`，则始终保持为 `BHP` 控制。

示例：

```text
WCONPROD
PROD*   OPEN    ORAT   20000.0     1000
/
```

## WELTARG<span id=_WELTARG></span>(WELLTARG)  (e)(/)

WELTARG 用于改变井的控制方式，这是与**时间相关**的控制，见 [TSTEP](#_TSTEP)。在 WELTARG 中，需要匹配井的名字，方式同 [COMPDAT](#_COMPDAT)。

示例：

```text
WELTARG
PROD*  ORAT  100.0 
/
```

上述示例中，所有名字以 `PROD` 开头的井改为产油控制，速率为 100 stb/day。

## WELLSTRE<span id=_WELLSTRE></span> (e)(/)

WELLSTRE 用于组分模型，定义了注入流体的各组分摩尔占比。如果在 [WCONINJE](#_WCONINJE) 中，给出了某一非基本注入物(即非水)的名称，则需要在 WELLSTRE 中给出各组分的摩尔占比。示例：

```text
WCONINJE
INJ    Solvent    OPEN    RATE     4700     4000 
/

WELLSTRE
Solvent  0.6  0.3  0.05  0.02  0.03
/
```

各组分的摩尔占比和应为 1，末尾的 0 可以省略。

## CNAMES<span id=_CNAMES></span> (e)

CNAMES 用于输入烃组分的名字，需要先输入 NCNP

```text
CNAMES
Meth Ethane C3-C6 C7+
```

## TCRIT<span id=_TCRIT></span> (e)(/)

TCRIT 定义了烃组分的临界温度，单位 $^{\circ}$ R，如有多区域则以 / 区分，区域数量应与 [TABDIMS](#_TABDIMS) 中的 NTPVT 一致，示例：

```text
TCRIT
140 270 450 670  /
100 200 300 400
/
```

## PCRIT<span id=_PCRIT></span> (e)(/)

PCRIT 定义了烃组分的临界压力，单位 psia，格式同 [TCRIT](#_TCRIT)。

## VCRIT<span id=_VCRIT></span> (e)(/)

VCRIT 定义了烃组分的临界摩尔体积，单位 $\mathrm{ft}^{3}/\mathrm{lb\text{-}M}$，格式同 [TCRIT](#_TCRIT)。

## ZCRIT<span id=_ZCRIT></span> (e)(/)

ZCRIT 定义了烃组分的临界Z-factor，格式同 [TCRIT](#_TCRIT)。

## MW<span id=_MW></span> (e)(/)

MW 定义了烃组分的分子质量，单位 $\mathrm{lb/lb\text{-}M}$。格式同 [TCRIT](#_TCRIT)。

## ACF<span id=_ACF></span> (e)(/)

ACF 定义了烃组分的偏心因子，格式同 [TCRIT](#_TCRIT)。

## OMEGAA,OMEGAB<span id=_OMEGAA></span> <span id=_OMEGAB></span> (e)(/)

OMEGAA，OMEGAB 分别定义了用于状态方程计算的系数 $\Omega_{A},\Omega_{B}$，格式同 TCRIT。默认时 $\Omega_{A}=0.457235529,\Omega_{B}=0.077796074$，当前仅限于 PR 方程
，格式同 [TCRIT](#_TCRIT)。

## SSHIFT<span id=_SSHIFT></span> (e)(/)

SSHIFT 定义了烃组分的体积偏移参数，默认时值为 0。格式同 [TCRIT](#_TCRIT)。

## PARACHOR<span id=_PARACHOR></span> (e)(/)

PARACHOR 用于混溶模型的表面张量计算，需要打开 MISCIBLE 选项，单位 $\mathrm{(dynes/cm)}^{1/4}\mathrm{cc}/\mathrm{gm\text{-}M}$ ，格式同 [TCRIT](#_TCRIT)。

## VCRITVIS<span id=_VCRITVIS></span> (e)(/)

VCRITVIS 定义了仅用于粘性计算的临界摩尔体积，单位 $\mathrm{ft}^{3}/\mathrm{lb\text{-}M}$。如果未输入，则使用 ZCRIT 进行计算，若 ZCRIT 也没有输入，则赋值为 VCRIT。格式同 [TCRIT](#_TCRIT)。

## BIC<span id=_BIC></span> (e)(/)

BIC 关键字用于组分模型，定义了组分间的二元相互系数，它是一个实对称矩阵，示例：

```text
BIC
#BIC matrix
0    -.02     .1        .13    .135   0.1277    .1       .1       .1
-.02   0      .036      .05    .08    .1002     .1       .1       .1
.1     .036   0         0      0      .092810   .130663  .130663  .130663
.13    .05    0         0      0      0         .006     .006     .006
.135   .08    0         0      0      0         .006     .006     .006
.1277  .1002  .092810   0      0      0         0        0        0
.1     .1     .130663   .006   .006   0         0        0        0
.1     .1     .130663   .006   .006   0         0        0        0
.1     .1     .130663   .006   .006   0         0        0        0
/
```

也可以只输入下半部分(不含对角线)。

## LBCCOEF<span id=_LBCCOEF></span>

LBCCOEF 定义了使用 Lorentz-Bray-Clark 粘性计算公式时的参数，默认时为 0.1023, 0.023364, 0.058533, -0.040758, 0.0093324，格式同 [TCRIT](#_TCRIT)。

## RR<span id=_RR></span>

RR 关键字用于组分模型，给出了在相分裂计算中用 Newton 法求解 Rachford-Rice 方程的求解参数，包括：

1. 最大迭代步数
2. 残差控制

## SSMSTA<span id=_SSMSTA></span>

SSMSTA 关键字用于组分模型，给出了在相稳定性分析 SSM 的求解参数，包括：

1. 最大迭代步数
2. 残差控制

## NRSTA<span id=_NRSTA></span>

NRSTA 关键字用于组分模型，给出了在相稳定性分析中 Newton 法的求解参数，包括：

1. 最大迭代步数
2. 残差控制

## SSMSP<span id=_SSMSP></span>

SSMSP 关键字用于组分模型，给出了在相分裂计算中 SSM 的求解参数，包括：

1. 最大迭代步数
2. 残差控制

## NRSP<span id=_NRSP></span>

NRSP 关键字用于组分模型，给出了在相分裂计算中 Newton 法的求解参数，包括：

1. 最大迭代步数
2. 残差控制

-----

## SUMMARY<span id=_SUMMARY></span> (/)

SUMMARY 关键字用于控制输出每个时间步的各指标信息，其结果将保存在 `SUMMARY.out` 文件中。

对于一些基本指标会默认输出，包括：

* Time：当前时间点，day
* NRiter：到目前为止花费的总牛顿迭代次数
* LSiter：到目前为止花费的总线性迭代次数

可控制输出的变量如下，这些都是 SUMMARY 关键字下的**子关键字**。对于子关键字一共有三种格式，此处分别记为 S，N，L

* S：只要出现即激活该选项，例如：

  ```text
  SUMMARY
  FPR
  /
  ```

* N：需要给出具体井的名字，如果用 `/`，则表示输出所有井的信息。该类子关键字需要在下一行以 `/` 结束，例如：

  ```text
  SUMMARY
  WOPR
  INJE1  PROD2
  /
  WWPR
  /
  /
  ```

* L：需要给出具体的网格坐标。该类子关键字需要在下一行以 `/` 结束，例如：

  ```text
  SUMMARY
  BPR
  1 1 1
  2 3 4
  /
  /
  ```

属于 **S** 类的子关键字包括：

* FPR：平均油藏压力，psia
* FOPR：当前时间点的总产油速率，stb/day
* FOPT：总产油量，stb
* FGPR：当前时间点的总产气速率，Mscf/day
* FGPT：总产气量，Mscf
* FWPR : 当前时间点的总产水速率，stb/day
* FWPT：总产水量，stb
* FGIR：当前时间点的总注气速率，Mscf/day
* FGIT：总注气量，Mscf
* FWIR：当前时间点的总注水速率，stb/day
* FWIT：总注水量，stb

属于 **N** 类的子关键字包括：

* WBHP：指定井的井底压力(参考深度处的压力)，psia
* WOPR：指定当前时间点各井的产油速率，stb/day
* WGPR：指定当前时间点各井的产气速率，Mscf/day
* WWPR：指定当前时间点各井的产水速率，stb/day
* WOIR：指定当前时间点各井的注油速率，stb/day
* WGIR：指定当前时间点各井的注气速率，Mscf/day
* WWIR：指定当前时间点各井的注水速率，stb/day
* WOPT：指定井的总产油量，stb
* WGPT：指定井的总产气量，Mscf
* WWPT：指定井的总产水量，stb

属于 **L** 类的子关键字包括：

* DG：指定井的各射孔与井处压力的压差，psia
* BPR： 指定网格块的压力，psia
* SOIL： 指定网格块的油相饱和度
* SGAS： 指定网格块的气相饱和度
* SWAT： 指定网格块的水相饱和度

-----

## RPTSCHED<span id=_RPTSCHED></span> (/)

RPTSCHED 用来控制输出各个关键时间节点的详细信息，见 [TSTEP](#_TSTEP)。其结果将打印在`RPT.out`。

可控制的输出信息包括：

* PRE / PRESSURE : 各网格块的压力，psia
* PGAS：各网格块中气相压力，psia
* PWAT：各网格块中水相压力，psia
* DENO：各网格块中油相密度，lb/ft3
* DENG： 各网格块中气相密度，lb/ft3
* DENW：各网格块中水相密度，lb/ft3
* SOIL：各网格块中油相饱和度，dimensionless
* SGAS：各网格块中气相饱和度，dimensionless
* SWAT：各网格块中水相饱和度，dimensionless

## VTKSCHED<span id=_VTKSCHED></span> 

VTKSCHED 会将所指定的信息在每一个 TSTEP 输出成 VTK 的文件格式，可输出的信息参见 [RPTSCHED](#_RPTSCHED)。

示例：

```text
RPTSCHED
PRES   SOIL   SWAT
/
```

-----

## 参考示例 (SPE1)

SPE1 是三相三组分经典黑油模型

```text
--TITLE
    SPE1 Case1 (Fixed BPP)

-- Original size 10x10x3 = 300
DIMENS
 10  10  3  / 
 
-- Black oil model
BLACKOIL

-- oil,gas,water and dissolved gas could exist
OIL
WATER
GAS
DISGAS

-- Field unit used
FIELD

-- Only one saturation region and one PVT region
TABDIMS
1   1

-- grid information

EQUALS
'DX'    1000   6*      /
'DY'    1000   6*      /
'DZ'    20     4* 1 1  /
'DZ'    30     4* 2 2  /
'DZ'    50     4* 3 3  /
'PORO'  0.3    4* 1 1  /
'PORO'  0.3    4* 2 2  /
'PORO'  0.3    4* 3 3  /
'PERMX' 500    4* 1 1  /
'PERMX' 50     4* 2 2  /
'PERMX' 200    4* 3 3  /
'PERMZ' 75     4* 1 1  /
'PERMZ' 35     4* 2 2  /
'PERMZ' 15     4* 3 3  /
'TOPS'  8325   4* 1 1  /
/

COPY
'PERMX' 'PERMY' 4* 1 3 /
/

SWOF 
0.12000    0.00000   1.00000    0.00000
0.18000    0.00001    .85000    0.00000
0.24000     .0732    0.70000    0.00000
0.32000     .1707    0.35000    0.00000
0.37000     .2317    0.20000    0.00000
0.42000     .2927    0.09000    0.00000
0.52000     .4146    0.02100    0.00000
0.57000     .4756    0.01000    0.00000
0.62000     .5366    0.00100    0.00000
0.72000     .6586    0.00010    0.00000
0.75000     .6951    0.00000    0.00000
1.00000    0.9000    0.00000    0.00000
/

SGOF
0.00       0.00000   1.00000     0.00000
0.02       0.00000   0.997       0.00000 
0.05       0.005     0.980       0.00000
0.12       0.025     0.700       0.00000
0.20       0.075     0.350       0.00000
0.25       0.125     0.200       0.00000
0.30       0.190     0.090       0.00000
0.40       0.410     0.021       0.00000
0.45       0.600     0.010       0.00000
0.50       0.720     0.001       0.00000
0.60       0.870     0.0001      0.00000
0.70       0.940     0.00000     0.00000
0.85       0.980     0.00000     0.00000
1.00       1.000     0.00000     0.00000
/

PVCO 
  14.7   0.0010      1.062       1.040       15.1E-6     0.46E-4
 264.7   0.0905      1.150       0.975       15.1E-6     0.46E-4
 514.7   0.1800      1.207       0.910       15.1E-6     0.46E-4
1014.7   0.3710      1.295       0.830       15.1E-6     0.46E-4
2014.7   0.6360      1.435       0.695       15.1E-6     0.46E-4
2514.7   0.7750      1.500       0.641       15.1E-6     0.46E-4
3014.7   0.9300      1.565       0.594       15.1E-6     0.46E-4
4014.7   1.2700      1.695       0.510       15.1E-6     0.46E-4
9014.7   1.3500      1.705       0.500       15.1E-6     0.46E-4
/

PVDG
  14.7   166.67      .0080                                        
 264.7    12.09      .0096                                        
 514.7     6.2741    .0112                                        
1014.7     3.1970    .0140                                        
2014.7     1.6141    .0189                                        
2514.7     1.2940    .0208                                        
3014.7     1.0800    .0228                                        
4014.7      .8110    .0268                                        
5014.7      .6490    .0309                                        
9014.7      .3859    .0470   
/

PVTW
4014.7      1.0     3E-6       0.3100    0.0  /
/

ROCK
4014.7      0.3000E-05    /

GRAVITY
59.53       1.000987      0.792   /

EQUIL
8500  4825.22  8500  0  7000  0  1 /

PBVD
5000    4014.7    
9000    4014.7
/

SUMMARY
EXCEL
FPR
FOPR
FOPT
FGPR
FGPT
FWPR
FWPT
FGIR
FGIT
FWIR
FWIT
FWCT
FWPT
BPR 
1,1,1 /
10,10,3 /
/
WBHP 
/
WPI 
/

-- Well information

WELSPECS
'INJE1'   'G'   1   1     1*    'GAS'   /
'PROD1'   'G'   10  10    1*    'OIL'   /
/

COMPDAT
'INJE**'   2*   1   1     1*   0.5   3*   /
'PROD1'    2*   3   3     1*   0.5   3*   /
/

WCONINJE
'INJE*'   'GAS'   'OPEN'   'RATE'   100000.0      10000    /
/

WCONPROD
'PROD*'   'OPEN'    'ORAT'   20000.0     1000    /
/

TUNING
--  Init     max    min   incre   chop   cut
     0.1     10     0.1     5     0.3    0.3                 /
--  dPlim  dSlim   dNlim   dVerrlim
     300     0.2    0.3     0.001                            /
-- itNRmax  NRtol  dPmax   dSmax  dPmin   dSmin   dVerrmax
     10     1E-3    200     0.2    1E-0    1E-2    0.01      /
/

-- use IMPEC and defaulted linear solver
METHOD
IMPEC
/

-- 10 years simulation
TSTEP
1    3    9    29    8  
/

TSTEP
132.625   182.625   185.625  
/

TSTEP
3*182.625   
/

TSTEP
7*365.25
/
```

## 参考示例 (SPE5)

SPE5 是三相十组分组分模型，其中九个烃组分，先无注入生产两年，再交替注入水/气十八年

```text
------------------------------------------------------------------------
-- SPE 16000
-- "Fifth Comparative Solution Project : Evaluation of Miscible Flood Simulators"
-- J.E. Killough, C.A. Kossack
--
-- The fifth SPE comparison problem , reported by Killough and Kossack
-- ( 9th SPE Symp on Res. Sim., San Antonio, 1987).
-- Dimension 7x7x3
-- 6-component compositional model, run on Ecl300
-- FIELD units
-- The run follows the first of the three suggested production schedules
-- Solvent injection as part of a WAG cycle
------------------------------------------------------------------------

TITLE
   SPE Fifth Comparison Test Problem - Scenario One
   
-- Field units used
FIELD

-- Compositional model used
COMPS
6
/

-- Only one saturation region and one PVT region
TABDIMS
1 1

-- Grid Size
DIMENS
7 7 3 /

-- GRID information  
EQUALS
'DX'     500   6*        /
'DY'     500   6*        /
'DZ'      20   4*  1  1  /
'DZ'      30   4*  2  2  /
'DZ'      50   4*  3  3  /
'PORO'   0.3   4*  1  3  /
'PERMX'  500   4*  1  1  /
'PERMX'   50   4*  2  2  /
'PERMX'  200   4*  3  3  /
'PERMZ'   50   4*  1  2  /
'PERMZ'   25   4*  3  3  /
'TOPS'  8325   4*  1  1  /
/

COPY
'PERMX' 'PERMY' 6*  /
/

-- Components information, see EGOIL.in as follows
INCLUDE
EGOIL.in
/

-- Reservoir temperature
RTEMP
160 
/

SWOF
--SW         KRW     KROW    PCOW
  0.2        0       1       0
  0.2899     0.0022  0.6769  0
  0.3778     0.018   0.4153  0
  0.4667     0.0607  0.2178  0
  0.5556     0.1438  0.0835  0
  0.6444     0.2809  0.0123  0
  0.7        0.4089  0       0
  0.7333     0.4855  0       0
  0.8222     0.7709  0       0
  0.9111     1       0       0
  1          1       0       0
/

SGOF
--SG         KRG     KROg    PCOG
       0.0       0.0   1.00000       0.0
 0.0500000       0.0 0.8800000       0.0
 0.0889000 0.0010000 0.7023000       0.0
 0.1778000 0.0100000 0.4705000       0.0
 0.2667000 0.0300000 0.2963000       0.0
 0.3556000 0.0500000 0.1715000       0.0
 0.4444000 0.1000000 0.0878000       0.0
 0.5333000 0.2000000 0.0370000       0.0
 0.6222000 0.3500000 0.0110000       0.0
 0.6500000 0.3900000       0.0       0.0
 0.7111000 0.5600000       0.0       0.0
 0.8000000   1.00000       0.0       0.0
/

PVTW
14.7   1.00    3.3E-06      0.7     0.00E-01
/

ROCK
3990.30 5E-06
/

GRAVITY
1*       1*           1*   /

EQUIL
8400 4000 9000 0 7000 0 1 1 0  /

PBVD
5000    4014.7    
9000    4014.7
/

ZMFVD
1000.0   0.5  0.03  0.07  0.2  0.15  0.05
10000.0  0.5  0.03  0.07  0.2  0.15  0.05 
/

-- SUMMARY output
SUMMARY
EXCEL
FPR
FOPR
FOPT
FGPR
FGPT
FWPR
FWPT
FGIR
FGIT
FWIR
FWIT
FWCT
FWPT
WBHP 
/
WPI 
/


-- RPTSCHED output
RPTSCHED
PRESSURE
DENO  DENG  DENW
SOIL  SGAS  SWAT
/

SCHEDULE    ==========================================================

TUNING
-- Init   max     min     incre    chop    cut
   1      10      0.1     5        0.3     0.3             /
-- dPlim  dSlim   dNlim   dVerrlim
   300    1       0.3     0.001                            /
-- itNRmax  NRtol  dPmax  dSmax  dPmin   dSmin   dVerrmax
   20       1E-3   200    0.2    1       1E-2    0.01      /
/

-- Use FIM and defaulted linear solver
METHOD
FIM 
/

-- Well information

WELSPECS
--name  group   I   J  depth_ref phase_ref
'PROD1'   'G'   7   7    1*    'OIL'   /
/

COMPDAT
--d
--name   I  J   K1  K2         diameter 
'PROD1'   2*    3   3        1*   0.5   3*   /
/

WCONPROD
--d
'PROD*'   'OPEN'  'ORAT'  12000    1000    /
/

--Start production only 

TSTEP
2*365.25
/

--Define injection well

-- Injected fluid 
WELLSTRE
Solvent 0.77 0.20 0.03 0.0 0.0 /
/

WELSPECS
I Field 1 1 8335 GAS /
/

COMPDAT
I 2* 1 1  1* 0.5 3* /
/

--Start WAG

WCONINJE
--d
--name type   openflag  mode  surface_rate       BHP
I     WATER   OPEN    RATE      12000         10000 /
/

TSTEP
1*365.25 
/

WCONINJE
--d
--name type   openflag  mode  surface_rate       BHP
I     Solvent   OPEN    RATE  12000           10000 /
/

TSTEP
1*365.25 
/

WCONINJE
--d
--name type   openflag  mode  surface_rate       BHP
I     WATER   OPEN    RATE      12000           10000 /
/

TSTEP
1*365.25
/

WCONINJE
--d
--name type   openflag  mode  surface_rate       BHP
I     Solvent   OPEN    RATE  12000           10000 /
/

TSTEP
1*365.25
/

WCONINJE
--d
--name type   openflag  mode  surface_rate       BHP
I     WATER   OPEN    RATE      12000          10000 /
/

TSTEP
1*365.25  
/

WCONINJE
--d
--name type   openflag  mode  surface_rate       BHP
I     Solvent   OPEN    RATE  12000           10000 /
/


TSTEP
1*365.25
/

WCONINJE
--d
--name type   openflag  mode  surface_rate       BHP
I     WATER   OPEN    RATE      12000           10000 /
/

TSTEP
1*365.25 
/

WCONINJE
--d
--name type   openflag  mode  surface_rate       BHP
I     Solvent   OPEN    RATE  12000           10000 /
/

TSTEP
1*365.25
/

WCONINJE
--d
--name type   openflag  mode  surface_rate       BHP
I     WATER   OPEN    RATE      12000           10000 /
/

TSTEP
1*365.25 
/

WCONINJE
--d
--name type   openflag  mode  surface_rate       BHP
I     Solvent   OPEN    RATE  12000           10000 /
/

TSTEP
1*365.25
/

WCONINJE
--d
--name type   openflag  mode  surface_rate       BHP
I     WATER   OPEN    RATE      12000           10000 /
/

TSTEP
1*365.25  
/

WCONINJE
--d
--name type   openflag  mode  surface_rate       BHP
I     Solvent   OPEN    RATE  12000          10000 /
/

TSTEP
1*365.25 
/

WCONINJE
--d
--name type   openflag  mode  surface_rate       BHP
I     WATER   OPEN    RATE      12000          10000 /
/

TSTEP
1*365.25  
/

WCONINJE
--d
--name type   openflag  mode  surface_rate       BHP
I     Solvent   OPEN    RATE  12000           10000 /
/

TSTEP
1*365.25
/

WCONINJE
--d
--name type   openflag  mode  surface_rate       BHP
I     WATER   OPEN    RATE      12000           10000 /
/

TSTEP
1*365.25 
/

WCONINJE
--d
--name type   openflag  mode  surface_rate       BHP
I     Solvent   OPEN    RATE  12000           10000 /
/

TSTEP
1*365.25
/

WCONINJE
--d
--name type   openflag  mode  surface_rate       BHP
I     WATER   OPEN    RATE      12000          10000 /
/

TSTEP
1*365.25  
/

WCONINJE
--d
--name type   openflag  mode  surface_rate       BHP
I     Solvent   OPEN    RATE  12000          10000 /
/

TSTEP
1*365.25 
/


CNAMES
C1
C3
C6
C10
C15
C20
/

TCRIT
 343.0
 665.7
 913.4
1111.8
1270.0
1380.0
/

PCRIT
667.8
616.3
436.9
304.0
200.0
162.0
/

ZCRIT
0.290
0.277
0.264
0.257
0.245
0.235
/

MW
 16.04
 44.10
 86.18
149.29
206.00
282.00
/

ACF
0.013
0.1524
0.3007
0.4885
0.6500
0.8500
/

BIC
   0.0
   0.0    0.0
   0.0    0.0    0.0
   0.05   0.005  0.0     0.0
   0.05   0.005  0.0     0.0     0.0 /



RR
#Rachford Rice equation parameters
#maxit tolerance
   30	   1e-12
/

SSMSTA
#Successive substitution for stability analysis
#maxit  tolerance  eYt(relaxation factor)
  100	   1e-12     1e-8
/

NRSTA
#Newton Raphson for stability analysis
#maxit tolerance 
  55     1e-12
/


SSMSP
#SSM for phase splitting calculation
#maxit toleranceR 
  100	    1E-6
/

NRSP
#Newton Raphson for phase splitting calculation
#maxit toleranceR 
  55      1e-12
/
```

<!-- <script type="text/javascript" src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"></script> -->
<!-- <script type="text/x-mathjax-config">
  MathJax.Hub.Config({ tex2jax: {inlineMath: [['$', '$']]}, messageStyle: "none" });
</script> -->