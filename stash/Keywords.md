# Introduction to OpenCAEPoro Keywords

OpenCAEPoro 目前用一套关键字体系来构建模拟选项，控制模拟过程，大部分关键字格式与内容与斯伦贝谢的 Eclipse 保持一致，对于这些关键字，将在其后用 **(e)** 或 **(e300)** 标出，**(e300)** 表示仅和适用于 E300 的选项匹配。有些关键字后面还会用 **(/)** 标出，这表示必须在其内容结束的**下一行**使用 `/` 表示其已结束，例如

```
PVTW (e)(/)
4014.7      1.0     3E-6       0.3100    0.0
/
```

文本中的 `'` 和 `,` 将会自动替换成空格。**(Todo)** 标识了目前尚不完善的地方。

注意：目前，从关键字读入的数据只适用于**英制**单位。

## 注释

在 OpenCAEPoro 中，有以下几种注释方式

* 以 `#` 或 `--` 开头的行，例如

  ```
  # comments
     # comments
  -- comments 
     -- comments
  ```

* 在任何 `/` 的背后，例如

  ```
  DIMENS
  6  6  6  / comments
  ```

注意：不应在关键字后面注释，否则该关键字将被忽略，例如

```
DIMENS  /  comments
```



## <span id=_BLACKOIL>BLACKOIL</span> (e)

BLACKOIL 用来定义黑油模型，示例

```
BLACKOIL
```



## <span id=_OIL>OIL</span> (e)

OIL 表明在模拟过程中油相可能存在，用于黑油模型，示例

```
OIL
```



## <span id=_GAS>GAS</span> (e)

GAS 表明在模拟过程中气相可能存在，用于黑油模型，示例

```
GAS
```



## <span id=_WATER>WATER</span> (e)

WATER 表明在模拟过程中水相可能存在，当前的组分模型中，水相一定存在，示例

```
WATER
```



## <span id=_DISGAS>DISGAS</span> (e)

DISGAS 表明在模拟过程中溶解气可能存在于油相中，用于黑油模型，示例

```
DISGAS
```



## <span id=_COMPS>COMPS</span>

COMPS 激活了组分模型，可以输入组分数，但并不起作用，烃组分数和最大烃相数目前由 NCNP 输入，示例

```
COMPS
6
```



## <span id=_DIMENS>DIMENS</span> (e)(/)

DIMENS 用来定义储层网格的规模，适用于三维结构网格，用三个数字分别给出了 x，y，z 方向上网格的维度，示例

```
DIMENS
3  4  5
```

注意，向下为正，x，y，z 轴符合左手螺旋定则。



## <span id=_TABDIMS>TABDIMS</span>

TABDIMS 定义了饱和度表格和PVT表格的最大数量，不同的饱和度表格和PVT表格可能用在不同的网格区域，这将通过 [SATNUM](#_SATNUM) 和 [PVTNUM](#_PVTNUM) 来定义。示例

```
TABDIMS
-- NTSFUN  NTPVT
1  2
```

表明储层中至少有一个饱和度计算区域和两个 PVT 计算区域。



## EQUALS<span id=_EQUALS> </span> (e)(/)

EQUALS 用来批量定义油藏的性质，其中的子关键字具有如下的一般格式

```
item  var  i1  i2  j1  j2  k1  k2
```

其中 item 表示要赋值的性质的名称，var 表示要赋值的值，后面六个数表示赋值的网格范围，依次是 x 坐标的范围，y 坐标的范围 和 z 坐标的范围，例如

```
DX  1000  1  3  1  3  1  3
```

表示油藏区域 （i，j，k)，其中 [1,3]x[1,3]x[1,3] 的网格块的 x 方向的长度为 1000 英尺。坐标范围的值还支持简写，假设三维的坐标范围分别为 1-N1，1-N2，1-N3，则

```
Dx  1000  6*              / Dx = 1000 in all Region
Dx  1000  1  3  2*  2  4  / DX = 1000 in [1,3]x[1,N2]x[2,4]
Dx  1000  4* 2  6         / DX = 1000 in [1,N1]x[1,N2]x[2,6]
```

目前，在EQUALS中，支持如下子关键字的设置

```
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
TOPS      / 顶层网格块的上表面的深度，feet。由于顶层的 z 坐标为1，因此后两个数字为1。
```



## <span id=_COPY>COPY</span> (e)(/)

COPY 用于对 [EQUALS](#_EQUALS) 中子关键字的内容进行相互拷贝，其坐标格式与 EQUALS 中一致，例如

```
PERMX  PERMY   1  3  1  3  1  3
```

表示将区域[1,3]x[1,3]x[1,3]中的网格块的PERMY属性的值设为PERMX属性的值。



## <span id=_MULTIPLY>MULTIPLY</span> (e)(/)

MULTIPLY关键字用于对 [EQUALS](#_EQUALS) 中子关键字的内容进行放缩，其坐标格式与 EQUALS 中一致，例如

```
PERMZ  0.1  1  3  1  3  1  3
```

表示将区域[1,3]x[1,3]x[1,3]中的网格块的PERMZ属性的值乘以0.1。



## DX DY DZ PORO NTG PERMX PERMY PERMZ TOPS (e)(/)

关键字用于定义每个网格块的对应属性的值，见 [EQUALS](#_EQUALS)，网格块依照 x -> y -> z 的坐标字典序。

```
DX
10  12  6*20
/
```



## <span id=_SWITINIT>SWITINIT</span> (e)(/)
SWITINIT 用于定义初始水相饱和度，并以此缩放毛管力曲线使得水相饱和度与初始平衡解一致
```
SWATINIT
0.2  0.25  6*0.3
/
```



## <span id=_SATNUM>SATNUM</span> (e)(/)
SATNUM 用来定义饱和度表格作用区域，即指定网格块使用相应的饱和度表格
```
SATNUM
1  2   10*3
/
```



## <span id=_ACTNUM>ACTNUM</span> (e)(/)
ACTNUM 用来定义网格块的状态，目前，0 表示死网格，1表示活网格
```
ACTNUM
1  0   10*1
/
```



## <span id=_PVTNUM>PVTNUM</span> (e)(/)
PVTNUM 用来定义PVT区域，即指定网格块使用相应的PVT关系
```
PVTNUM
1  2   10*3
/
```



## <span id=_COORD>COORD</span> (e)(/)

COORD 用于角点网格的定义，它给出了角点网格的“骨架”，如果网格规模为 Nx * Ny * Nz，那么在 COORD 中就需要给出(Nx + 1)(Ny + 1) 条纵轴方向的线，按照 x -> y 的顺序，而每条线由两个三维点表示，顶层的点在前，底层的点在后，往下为正，单位英尺，例如，对于 Nx = 3，Ny = 2

```
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



## <span id=_ZCORN>ZCORN</span> (e)(/)

ZCORN 用于角点网格的定义，它通过给出每个网格的8个角点在对应“骨架”(见 [COORD](#_COORD) )上的深度。

```
 3 -- 4
/|   /|
1 -- 2|
||   ||
| 7-- 8
|/   /
5 -- 6
```
给出顺序为：
* 网格层： (k=1 -> k=nZ)
* 每一行网格(k固定，j -> nY)：给出各网格块节点 1 和 2 的深度，先给出 j = 1 行，$v_{1,1},v_{2,1},...,v_{1,nX},v_{2,nX}$，其后是 $v_{3,1},v_{4,1},...,v_{3,nX},v_{4,nX}$，再依次给出后面的行
* 上一点结束后，同样的方式给出 $v_{5,1},v_{6,1},...,v_{5,nX},v_{6,nX}$，其后是 $v_{7,1},v_{8,1},...,v_{7,nX},v_{8,nX}$




## <span id=_RTEMP>RTEMP</span> (e)

RTEMP 用来定义恒温油藏的温度，单位为 °F，示例

```
RTEMP
200
```



## <span id=_SWOF>SWOF</span> (e)(/)

SWOF 用来定义当油相和水相存在时的饱和度表格。如果有多个表格，则需用 `/` 隔开。在SWOF中一共由四列数据，水相的饱和度作为自变量，这四列数据分别表示

```
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

每一列的意义如下

* Sw : 水的饱和度，dimensionless
* Krw : 对应的水的相对渗透率，dimensionless
* Krow : 当仅有油相和水相存在时，对应的油相的相对渗透率，dimensionless
* Pcow : 对应的水相的毛管力 Po - Pw，psia



## <span id=_SWFN>SWFN</span> (e)(/)

SWFN 用来定义水相的饱和度表格。如果有多个表格，则需用 / 隔开，与 [SWOF](#_SWOF) 类似。在 SWFN 中一共有三列数据，水相的饱和度作为自变量，这三列数据分别表示
* Sw : 水相的饱和度，dimensionless
* Krw : 对应的水相相对渗透率
* Pcow : 对应的水相的毛管力 Po - Pw，psia



## <span id=_SGOF>SGOF</span> (e)(/)

SGOF 用来定义当气相和水相存在时的饱和度表格。如果有多个表格，则需用 `/` 隔开，与 [SWOF](#_SWOF) 类似。在SGOF中一共由四列数据，气相的饱和度作为自变量，这四列数据分别表示

* Sg : 水的饱和度，dimensionless
* Krg : 对应的气相的相对渗透率，dimensionless
* Krog : 当油相，气相和束缚水存在时，对应的油相相对渗透率，dimensionless
* Pcog : 对应的气相毛管力 Pg - Po，psia



## <span id=_SGFN>SGFN</span> (e)(/)

SGFN 用来定义气相的饱和度表格。如果有多个表格，则需用 / 隔开，与 [SWOF](#_SWOF) 类似。在 SGFN 中一共有三列数据，气相的饱和度作为自变量，这三列数据分别表示
* Sg : 气相的饱和度，dimensionless
* Krg : 对应的气相相对渗透率
* Pcog : 对应的气相毛管力 Pg - Po，psia



## <span id=_SOF3>SOF3</span> (e)(/)

SOF3 用来定义油相的饱和度表格(用于三相模型)。如果有多个表格，则需用 / 隔开，与 [SWOF](#_SWOF) 类似。在 SGFN 中一共有三列数据，油相的饱和度作为自变量，这三列数据分别表示
* So : 油相的饱和度，dimensionless
* Krog : 当仅有油相和水相存在时，对应的油相的相对渗透率，dimensionless
* Krow : 当油相，气相和束缚水存在时，对应的油相相对渗透率，dimensionless



## <span id=_PVCO>PVCO</span> (e)(/)

PVCO 用于黑油模型，通过表格数据提供了活油（含有溶解气）的 PVT 性质。如果有多个表格，则需用 `/` 隔开，与 [SWOF](#_SWOF) 类似。PVCO一共有六列数据，其中，前两列作为自变量(泡点压力和气油比)，即其余的量随着它们之一的改变而改变。这六列分别表示

* Pbub : 油相的泡点压力，psia
* Rs : 油相的气油比，Mscf/stb。它等于标准状态下从油相中析出的气体的体积与油相的体积之比。
* Bo : 对应的 Pbub 压力下油相的体积系数，rb / stb。它等于相同摩尔数的油相在储层状态下与在标准状态下的体积之比。
* Viscosity : 对应的 Pbub 压力下油相的粘性，cP
* Cb : 对应的油气比 Rs 下不饱和油相的压缩性，1/psia，它等于 Bo 对于压力的相对变化率
* Cv : 对应的油气比 Rs 下不饱和油相的粘性的变化率，1/psia，它等于粘性对于压力的相对变化率



## <span id=_PVDO>PVDO</span> (e)(/)

PVDO 用于黑油模型，通过表格数据提供了死油 (不含有溶解气) 的PVT 性质。如果有多个表格，则需用 `/` 隔开，与 [SWOF](#_SWOF) 类似。PVDO 一共有三列数据，其中油相压力为自变量，这三列分别表示

* Po : 油相压力，psia
* Bo : 对应的油相体积系数，rb/stb
* Viscosity : 对应的油相粘性，cP



## <span id=_PVDG>PVDG</span> (e)(/)

PVDG 用于黑油模型，通过表格数据提供了干燥气的 PVT 性质。如果有多个表格，则需用 `/` 隔开，与 [SWOF](#_SWOF) 类似。PVCO 一共有三列数据，气相的压力作为自变量，这三列数据分别表示

* Pg : 气相压力，psia
* Bg : 对应的气相体积系数，rb/Mscf
* Viscosity : 对应的气相粘度，cP



## <span id=_PVTW>PVTW</span> (e)(/)

PVTW 用于黑油模型或者单独考虑水相的组分模型，通过表格数据提供了水相的 PVT 性质。如果有多个表格，则需用 `/` 隔开，与 [SWOF](#_SWOF) 类似。PVTW 一共有五列数据，水相的压力作为自变量，这五列数据分别表示

* Pw : 水相的压力，psia
* Bw : 对应的水相体积系数，rb/stb
* Cw : 对应的水相压缩性系数，1/psia，它等于体积系数关于压力的相对变化率。
* Viscosity : 对应的水相粘度，cP
* Cv : 对应的水相的粘度关于压力的相对变化率，1/psia



## <span id=_PBVD>PBVD</span> (e)(/)

PBVD 通过表格数据给出了初始油藏下，泡点压力关于位置深度的关系。如果有多个表格，则需用 `/` 隔开，与 [SWOF](#_SWOF) 类似。PBVD 一共有两列数据，分别是

* Depth : 采样点深度，feet
* Pbub : 对应的泡点压力，psia。




## <span id=_MISCIBLE>MISCIBLE</span> (e300)
在组分模型中，MISCIBLE 关键字开启了混溶选项。



## <span id=_MISCSTR>MISCSTR</span> (e)
在混溶模型中，对界面的表面张量进行设置。必须先打开 [MISCIBLE](#_MISCIBLE) 选项。需要输入三个值

* 最大混溶表面张量：当表面张量大于次值时，油气不混溶， dynes/cm
* 所期待的最大表面张量：dynes/cm (暂时没有用上)
* 用来放缩输入毛管力曲线的最大表面张量：dynes/cm



## <span id=_EQUIL>EQUIL</span> 

EQUIL 用于油藏初始条件的计算，它给出了参考深度处的油藏压力，油气接触面的深度及毛管力，油水接触面的深度及毛管力。它的六列依次为

* Dref : 参考深度，feet
* Pref : 参考深度处的压力，psia
* D1 : 某接触面深度，feet
* P1 : 对应的两相毛管力，psia
* D2 : 某接触面深度，feet
* P2 : 对应的两相毛管力，psia

根据存在的相的情况依次给出需要的值

* 如果存在油相和水相，则
  * D1 为油水接触面深度，P1 = Po - Pw
  * 如果气相也存在，则 D2 为油气接触面深度，P2 = Pg - Po
* 如果只存在气相和水相，则
  * D1 为气水接触面深度，P1 = Pg - Pw
  * D2 和 P2 无需给出

EQUIL 应与上述的饱和度表格与 PVT 表格一样，在不同的区域有不同的数值，但由于多区域功能目前并未完全实现，因此只读入一行数据。**(Todo)**



## <span id=_ROCK>ROCK</span> 

ROCK 给出了岩石的可压缩性，一共有两列数据，分别表示

* Pref : 参考压力，psia
* C : 岩石的可压缩性，1/psia。它是岩石体积关于关于压力的相对变化率。

ROCK 应与上述的饱和度表格与PVT表格一样，在不同的区域有不同的数值，但由于多区域功能目前并未完全实现，因此只读入一行数据。**(Todo)**



## <span id=_GRAVITY>GRAVITY</span> 

GRAVITY 用于黑油模型或单独考虑水相的组分模型，用于计算流体的密度(重力因子：密度乘以重力常数)，分别需要键入油，水，气三相的重力因子

* oil：默认值 45.5
* water (with reference to pure water)：默认值 1.0
* gas (with reference to air)：默认值 0.7773

对于组分模型，只需给出水相的信息即可

GRAVITY 应与上述的饱和度表格与PVT表格一样，在不同的区域有不同的数值，但由于多区域功能目前并未完全实现，因此只读入一行数据。**(Todo)**



## <span id=_DENSITY>DENSITY</span>

DENSITY 用于黑油模型或单独考虑水相的组分模型，用于计算流体的密度(重力因子：密度乘以重力常数)，分别需要键入油，水，气三相的密度

* oil：默认值 37.457 lb/ft3
* water：默认值 62.366416 lb/ft3
* gas：默认值 0.062428 lb/ft3

对于组分模型，只需给出水相的信息即可

DENSITY 应与上述的饱和度表格与PVT表格一样，在不同的区域有不同的数值，但由于多区域功能目前并未完全实现，因此只读入一行数据。**(Todo)**



## <span id=_INCLUDE>INCLUDE</span> (e)

INCLUDE 关键字用于分文件编写输入文件，需要输入被包含文件的  **相对路径**，被包含文件的格式与主输入文件是一致的，示例

```
INCLUDE
./others.out
```



## <span id=_METHOD>METHOD</span>

METHOD 关键字用来确定所使用的的离散方法以及所调用的线性求解器，离散方法包括 IMPEC (隐式压力显式组分)，FIM (全隐式类方法)，线性求解器默认使用FASP (如需使用别的求解器，需要在程序中补充对应的接口)，因此，对应地需要给出求解常量矩阵或块状矩阵的FASP输入文件 (**相对路径**)，例如

```
METHOD
FIM ./bsr.fasp
```

模拟器默认使用 IMPEC 方法，而FASP输入文件的优先级为

1. Method 中给定
2. ./bsr.fasp (FIM)，./csr.fasp (IMPEC)
3. ../conf/bsr.fasp (FIM)，../conf/csr.fasp (IMPEC)
4. 内置参数



## <span id=_TSTEP>TSTEP</span> (e)(/)

TSTEP 关键字通过时间间隔给出了模拟的关键时间节点 (day)，这是在模拟中会强制到达的时间点。在这些时间节点上，井的控制方式可能会发生改变，油藏状态可能会进入不同的阶段 (根据经验预估)，由此求解参数可能会做出调整。第 0 天始终为第 1 个时间节点，因此无需再输入。TSTEP 关键字支持简写，例如 2*100 表示两个 100 天的间隔。于是下面的关键时间节点为第 0, 5, 15, 45, 145, 245 天。

```
TSTEP
5  10  30  2*100
/
```

在输入文件中一般包含多个 TSTEP 关键字。对同一对象的控制作用范围为这两个控制关键字中间的时间段，如下，第一个 [TUNING](#_TUNING) 的作用范围为中间的 115 天
```
TUNING
******
/

TSTEP
5 10
/

TSTEP
100
/

TUNING
******
/
```
此外，如果激活了 [RPTSCHED](#_RPTSCHED) 关键字，则会在每一个时间节点打印关键字里对应的信息。



## <span id=_TUNING>TUNING</span> (/)

TUNING 关键字给出了模拟求解的参数，包括时间步长控制参数和非线性求解控制参数，它是与**时间相关**的控制，见 [TSTEP](#_TSTEP)。TUNING 关键字的内容分为三部分，每一部分内容结束后需用 `/` 结尾，前两部分与时间步长的控制相关，最后一部分与非线性求解相关，主要用于 FIM 类方法。示例

```
TUNING
-- Init     max    min   incre   chop    cut
     1      10    0.1      3    0.15    0.3 /
--  dPlim  dSlim   dNlim   dVerrlim
     300     0.5    0.3    0.001           /
-- itNRmax  NRtol  dPmax  dSmax  dPmin   dSmin   dVerrmax
       10    1E-2   2000    2    1      0.01    0.01  /
/
```

* 第一部分以时间步长控制相关，最重要的参数为前三个
  * Init：每一个关键时间节点开始的最大时间步长，day。
  * max：相邻两个关键时间节点中，能取到的最大时间步长，day。
  * min：相邻两个关键时间节点中，下一个时间步长的预测下限。在自适应时间步长选取下，这给出了下一个时间步预测的最小值，但有时为了保证解的稳定性和合理性，时间步长会再缩短，因此真实的时间步长是可能会低于预测下限，day。
  * incre：在自适应时间步长选取下，下一个时间步长相对于当前时间步长的最大增长因子，dimensionless。
  * chop：在自适应时间步长选取下，下一个时间步长相对于当前时间步长的最小缩减因子，dimensionless。
  * cut：在自适应时间步长选取下，在 FIM 方法中，当 Newton 方法失败时时间步长的缩短因子，dimensionless。
* 第二部分以下一时间步长的预测为主，它是假定变量基于时间线性变化的线性预测，它根据以下变量进行预测，并不保证下一时间步的真实情况也满足条件
  * dPlim：下一时间步每个网格块能接受的最大压力变化，psia
  * dSlim：下一时间步每个网格块每一相能接受的最大饱和度变化，dimensionless
  * dNlim：下一时间步每个网格块每种组分能接的最大摩尔数相对变化，dimensionless
  * dVerrlim：下一时间步每个网格块能接受的孔隙体积与流体体积的相对误差的最大值，dimensionless
* 第三部分与非线性迭代的控制相关，它主要用在 FIM 类方法求解非线性方程的 Newton 迭代中
  * itNRmax：每个时间步最大的 Newton 迭代步数，如果迭代超过该步数仍未收敛，则缩短时间步长重新迭代。
  * NRtol：Newton 迭代的收敛误差，目前在 OpenCAEPoro 中，这个值用于多种误差控制，包括残差的无穷范数、相对于网格块孔隙体积的残差的无穷范数、相对于网格块组分摩尔总数的残差的无穷范数。
  * dPmax：每个 Newton 迭代最大的压力变化，psia
  * dSmax：每个 Newton 迭代步最大的饱和度变化，由于从组分摩尔数到相饱和度需要经历复杂的相平衡计算(组分模型)，因此基于程序性能的考量，所使用的的饱和度变化是一个预测值而非真实值。
  * dPmin，dSmin：Newton 迭代最小的压力变化和饱和度变化。当一步 Newton 迭代的最大压力变化和最大饱和度变化 (均为真实值) 分别低于 dPmin 和 dSmin 时，则认为 Newton 迭代收敛。
  * dVerrmax：每个网格块孔隙体积与流体体积的相对误差，此选项一般用于 IMPEC 类方法，因为 FIM 类方法对此具有较好的保证。



## <span id=_WELSPECS>WELSPECS</span> (/)

WELSPECS 用来输入井的信息，包括

* 井的名字
* 所属的井组 (暂时无用处)
* 井所在网格块的 x 轴坐标
* 井所在网格块的 y 轴坐标
* 井的深度，feet

```
WELSPECS
INJE1   G   1   1     1*
PROD1   G   10  10    1*
/
```

在井的信息输入中，可能会使用默认值，如上面的 `1*`，表示对应的位置是默认值。如果是 `m*`，则表示对应的 m 个位置的值为默认值。

在 WELSPECS 中，只有井的深度允许为默认值，当它被设为**默认值**或者**负数**时，井的位置位于第 0 个射孔的位置，关于射孔的排序，见 [COMPDAT](#_COMPDAT)



## <span id=_COMPDAT>COMPDAT</span> (/)

COMPDAT 用于定义井的射孔，使用该关键字时，需要将所定义的射孔的信息匹配到对应的井上去，目前有精确匹配和模糊匹配两种匹配方式。

* 精确匹配

  需要准确的井的名字，如

  ```
  COMPDAT
  INJ1    ... 
  /
  ```

  会将其后定义的信息准确匹配到井 `INJ1`。

* 模糊匹配

  在末尾带上 `*`，则会匹配到名称以前面部分开头的井，如

  ```
  COMPDAT
  INJ*    ... 
  /
  ```

  则会匹配到所有名字以 `INJ` 开头的井。

在 COMPDAT 中，包括了

* 所需匹配的井的信息
* 射孔所在网格块的 x 轴坐标，默认则为所属井的 x 轴坐标
* 射孔所在网格块的 y 轴坐标，默认则为所属井的 y 轴坐标
* 射孔所在网格块的 z 轴坐标范围，由两个升序的正整数组成。
* 射孔与网格块连接处的传播因子(WI)，rb/day/psia，默认或为负则由程序自行计算得到。
* 射孔的直径，feet，默认为1  feet
* 射孔与网格块连接处的等效渗透因子(Kh)，mD*ft，默认或为负则由程序自行计算得到。
* 射孔的 Skinfactor，默认则为 0
* 射孔的方向：x，y，z。默认为 z，即竖直方向。

示例

```
COMPDAT
INJE*   2*   1   2     1*   0.5   3*   
PROD1   1   3   3   3     1*   0.5   3*    
/
```

对于名字中以 `INJE` 开头的井，它拥有与井横纵坐标相同的两个射孔，这两个射孔的纵坐标分别为 1 和 2，射孔的直径为 0.5 ft，其余皆为默认值。

对于井 PROD1，它拥有一个坐标位置在 (1,3,3) 的一个射孔，射孔直径为 0.5ft，其余皆为默认值。

在 OpenCAEPoro 中，要求给出的射孔顺序是由高至低的 (这由于射孔间静水压差的计算方式决定)，目前并没有在内部依据射孔深度对射孔进行排序，同时如果井的深度是默认的，则井的深度则为位置最高的射孔的深度。**(Todo)**

现有程序只需要在加上水平井相关参数的计算即可计算水平井的模型，但水平井的摩擦力并未考虑。



## <span id=_WCONINJE>WCONINJE</span> (/)

WCONINJE 用于对注入井的控制，这是与时间相关的控制，见 [TSTEP](#_TSTEP)。请注意，在OpenCAEPoro中，一口井是生产井还是注入井，这只是一个暂时的状态，而不是一开始就确定的结果，这意味着井的类型是可以随时切换的。因此，WCONINJE 不仅赋予了井的一般控制，也赋予了井的类型。

在 WCONINJE 中，同样需要匹配井的名字，方式同 [COMPDAT](#_COMPDAT)。 WCONINJE 中包含如下的控制

1. 所需匹配的井的信息
2. 注入流体的类型：包括水相 `WAT 或 WATER`，气相 `GAS`，在组分模型中，也可以是任意的流体，但需要在 [WELLSTRE](#_WELLSTRE) 中给出对应的组分比例。
3. 井的状态：开 `OPEN`，关 `CLOSE`
4. 井的控制方式：注入气体的流速控制 `RATE`，注入井的恒压控制(井的参考位置处) `BHP`
5. 流速控制的值或流速上限：注入气的单位为 Mscf/day，注入水的单位为 stb/day
6. 压力控制的值或压力上限：单位为 psia

如果井的控制方式是 `RATE`，则 5 为流速控制的值，6 为压力上限，当井的压力超过该上限时，则切换为 `BHP` 控制。如果井的控制方式是 `BHP`，则 6 为压力控制的值，5 为流速上限，当流速超过上限时，则切换为 `RATE` 控制。

如果井的控制方式是 `BHP`，且 5 是默认值，则始终保持为 `BHP` 控制。

示例

```
WCONINJE
INJE1   GAS   OPEN   RATE   100000.0      10000
INJE2   GAS   OPEN   BHP    1*            10000
/
```





## <span id=_WCONPROD>WCONPROD</span> (/)

WCONPROD  用于对生产井的控制，这是与**时间相关**的控制，见 [TSTEP](#_TSTEP)。请注意，在OpenCAEPoro中，一口井是生产井还是注入井，这只是一个暂时的状态，而不是一开始就确定的结果，这意味着井的类型是可以随时切换的。因此，WCONPROD 不仅赋予了井的一般控制，也赋予了井的类型。

在 WCONPROD 中，需要匹配井的名字，方式同 [COMPDAT](#_COMPDAT)。WCONPROD 中包含如下的控制

1. 所需匹配的井的信息
2. 井的状态：开 `OPEN`，关 `CLOSE`
3. 井的控制方式：产油速率控制 `ORAT`，产气速率控制`GRAT`，产水速率控制`WRAT`，井的恒压控制(井的参考位置处) `BHP`
4. 流速控制的值或流速上限：产出气的单位为 Mscf/day，产出水或油的单位为 stb/day
5. 压力控制的值或压力上限：单位为 psia

如果井的控制方式是 `*RAT`，则 4 为流速控制的值，5 为压力下限，当井的压力低于该下限时，则切换为 `BHP` 控制。如果井的控制方式是 `BHP`，则 5 为压力控制的值，4 为流速上限，当流速超过上限时，则切换为 `RATE` 控制。

如果井的控制方式是 `BHP`，且 4 是默认值，则始终保持为 `BHP` 控制。

示例

```
WCONPROD
PROD*   OPEN    ORAT   20000.0     1000
/
```



## <span id=_WELTARG>WELTARG</span>(WELLTARG)  (e)(/)

WELTARG 用于改变井的控制方式，这是与**时间相关**的控制，见 [TSTEP](#_TSTEP)。在 WELTARG 中，需要匹配井的名字，方式同 [COMPDAT](#_COMPDAT)。

示例

```
WELTARG
PROD*  ORAT  100.0 
/
```

上述示例中，所有名字以 `PROD` 开头的井改为产油控制，速率为 100 stb/day。



## <span id=_WELLSTRE>WELLSTRE</span> (e)(/)

WELLSTRE 用于组分模型，定义了注入流体的各组分摩尔占比。如果在 [WCONINJE](#_WCONINJE) 中，给出了某一非基本注入物(即非水)的名称，则需要在 WELLSTRE 中给出各组分的摩尔占比。示例

```
WCONINJE
INJ    Solvent    OPEN    RATE     4700     4000 
/

WELLSTRE
Solvent  0.6  0.3  0.05  0.02  0.03
/
```

各组分的摩尔占比和应为 1，末尾的 0 可以省略。



## <span id=_NCNP>NCNP</span>

NCNP 关键字用于组分模型，它给出了组分的数量(不包括水组分)和烃相的最大数量。示例

```
NCNP
6 2
```



## <span id=_CNAMES>CNAMES</span> (e)
CNAMES 用于输入烃组分的名字，需要先输入 NCNP
```
CNAMES
Meth Ethane C3-C6 C7+
```



## <span id=_ZI>ZI</span> (e)

ZI 用组分模型，它给出了油藏(每个网格块)初始的组分摩尔占比(不包括水组分)

目前仍是简单的初始化，没有涉及关于深度的插值变化 **(Todo)**

```
ZI
0.5 0.03 0.07 0.2 0.15 0.05
```



## <span id=_TCRIT>TCRIT</span> (e)(/)
TCRIT 定义了烃组分的临界温度，单位 $^{\circ}$ R，如有多区域则以 / 区分，区域数量应与 NTPVT 一致，示例
```
140 270 450 670  /
100 200 300 400
/
```



## <span id=_PCRIT>PCRIT</span> (e)(/)
PCRIT 定义了烃组分的临界压力，单位 psia，格式同 [TCRIT](#_TCRIT)



## <span id=_VCRIT>VCRIT</span> (e)(/)
VCRIT 定义了烃组分的临界摩尔体积，单位 $\mathrm{ft}^{3}/\mathrm{lb\text{-}M}$，格式同 [TCRIT](#_TCRIT)



## <span id=_ZCRIT>ZCRIT</span> (e)(/)
ZCRIT 定义了烃组分的临界Z-factor，格式同 [TCRIT](#_TCRIT)



## <span id=_MW>MW</span> (e)(/)
MW 定义了烃组分的分子质量，单位 $\mathrm{lb/lb\text{-}M}$。格式同 [TCRIT](#_TCRIT)



## <span id=_ACF>ACF</span> (e)(/)
ACF 定义了烃组分的偏心因子，格式同 [TCRIT](#_TCRIT)



## <span id=_OMEGAA>OMEGAA</span>，<span id=_OMEGAB>OMEGAB</span> (e)(/)
OMEGAA，OMEGAB 分别定义了用于状态方程计算的系数 $\Omega_{A},\Omega_{B}$，格式同 TCRIT。默认时 $\Omega_{A}=0.457235529,\Omega_{B}=0.077796074$，当前仅限于 PR 方程
，格式同 [TCRIT](#_TCRIT)


## <span id=_SSHIFT>SSHIFT</span> (e)(/)
SSHIFT 定义了烃组分的体积偏移参数，默认时值为 0。格式同 [TCRIT](#_TCRIT)



## <span id=_PARACHOR>PARACHOR</span> (e)(/)
PARACHOR 用于混溶模型的表面张量计算，需要打开 MISCIBLE 选项，单位 $\mathrm{(dynes/cm)}^{1/4}\mathrm{cc}/\mathrm{gm\text{-}M}$ ，格式同 [TCRIT](#_TCRIT)



## <span id=_VCRITVIS>VCRITVIS</span> (e)(/)
VCRITVIS 定义了仅用于粘性计算的临界摩尔体积，单位 $\mathrm{ft}^{3}/\mathrm{lb\text{-}M}$。如果未输入，则使用 ZCRIT 进行计算，若 ZCRIT 也没有输入，则赋值为 VCRIT。格式同 [TCRIT](#_TCRIT)


## <span id=_BIC>BIC</span> (e)(/)

BIC 关键字用于组分模型，定义了组分间的二元相互系数，它是一个实对称矩阵，示例

```
BIC
#BIC matrix
0	   -.02	   .1	     .13    .135	0.1277	  .1	   .1	    .1	
-.02   0       .036      .05    .08     .1002     .1       .1       .1
.1     .036     0        0      0       .092810   .130663  .130663  .130663
.13    .05      0        0      0       0         .006     .006     .006
.135   .08      0        0      0       0         .006     .006     .006
.1277  .1002   .092810   0      0       0         0        0        0
.1     .1      .130663   .006   .006    0         0        0        0
.1     .1      .130663   .006   .006    0         0        0        0
.1     .1      .130663   .006   .006    0         0        0        0
/
```
也可以只输入下半部分(不含对角线)



## <span id=_LBCCOEF>LBCCOEF</span>
LBCCOEF 定义了使用 Lorentz-Bray-Clark 粘性计算公式时的参数，默认时为 0.1023, 0.023364, 0.058533, -0.040758, 0.0093324，格式同 [TCRIT](#_TCRIT)


## <span id=_RR>RR</span>

RR 关键字用于组分模型，给出了在相分裂计算中用 Newton 法求解 Rachford-Rice 方程的求解参数，包括

1. 最大迭代步数
2. 残差控制



## <span id=_SSMSTA>SSMSTA</span>

SSMSTA 关键字用于组分模型，给出了在相稳定性分析 SSM 的求解参数，包括

1. 最大迭代步数
2. 残差控制



## <span id=_NRSTA>NRSTA</span>

NRSTA 关键字用于组分模型，给出了在相稳定性分析中 Newton 法的求解参数，包括

1. 最大迭代步数
2. 残差控制



## <span id=_SSMSP>SSMSP</span>

SSMSP 关键字用于组分模型，给出了在相分裂计算中 SSM 的求解参数，包括

1. 最大迭代步数
2. 残差控制



## <span id=_NRSP>NRSP</span>

NRSP 关键字用于组分模型，给出了在相分裂计算中 Newton 法的求解参数，包括

1. 最大迭代步数
2. 残差控制



-----

## <span id=_SUMMARY>SUMMARY</span> (/)

SUMMARY 关键字用于控制输出每个时间步的各指标信息，其结果将保存在 `SUMMARY.out` 文件中。

对于一些基本指标会默认输出，包括：
* Time：当前时间点，day
* NRiter：到目前为止花费的总牛顿迭代，dimensionless
* LSiter：到目前为止花费的总线性迭代，dimensionless

可控制输出的变量如下，这些都是 SUMMARY 关键字下的**子关键字**，对于子关键字一共有三种格式，此处分别记为 A，B，C

* A：只要出现即激活该选项，例如

  ```
  SUMMARY
  FPR
  /
  ```

* B：需要给出具体井的名字，如果用 `/`，则表示输出所有井的信息。该类子关键字需要在下一行以 `/` 结束，例如

  ```
  SUMMARY
  WOPR
  INJE1  PROD2
  /
  WWPR
  /
  /
  ```

* C：需要给出具体的网格坐标。该类子关键字需要在下一行以 `/` 结束，例如

  ```
  SUMMARY
  BPR
  1 1 1
  2 3 4
  /
  /
  ```

属于 **A** 类的子关键字包括

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

属于 **B** 类的子关键字包括

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

属于 **C** 类的子关键字包括

* DG：指定井的各射孔与井处压力的压差，psia

* BPR： 指定网格块的压力，psia
* SOIL： 指定网格块的油相饱和度
* SGAS： 指定网格块的气相饱和度
* SWAT： 指定网格块的水相饱和度

-----



## <span id=_RPTSCHED>RPTSCHED</span> (/)

RPTSCHED 用来控制输出各个关键时间节点的详细信息，见 [TSTEP](#_TSTEP)。其结果将打印在`RPT.out`

可控制的输出信息包括

* PRE / PRESSURE : 各网格块的压力，psia
* PGAS：各网格块中气相压力，psia
* PWAT：各网格块中水相压力，psia
* DENO：各网格块中油相密度，lb/ft3
* DENG： 各网格块中气相密度，lb/ft3
* DENW：各网格块中水相密度，lb/ft3
* SOIL：各网格块中油相饱和度
* SGAS：各网格块中气相饱和度
* SWAT：各网格块中水相饱和度

-----

示例

```
RPTSCHED
PRES   SOIL   SWAT
/
```

