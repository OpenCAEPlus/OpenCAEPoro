### SPE1A

![compare](D:\Lsz\OpenCAEPoro\data\figure\spe1a.png)

```
---------------------------------------------------
PS_IMPES_1

0.1  1  0.1

Final time:          3655.5 Days
Total time steps:    5533
Total nonlin steps:  5533
Total linear steps:  8828
Linear solver time:  0.932773s
Simulation time:     2.44659s

---------------------------------------------------
OCP_IMPEC_1

.\OpenCAEPoro.exe ..\..\data\spe1a\spe1a.data  IMPEC  0.1  1  0.1

Final time:          3655.500 Days
Total time steps:    5271
Total Newton steps:  5271 (+39 wasted steps)
Total linear steps:  8327 (+73 wasted steps)
Linear solve time:   0.816s (57.926%)
Simulation time:     1.408s

---------------------------------------------------
PS_FIM_10

1  10  0.1

Final time:          3655.5 Days
Total time steps:    420
Total Newton steps:  434
Wasted Newton steps: 1
Total linear steps:  1235
Linear solver time:  3.62965s
Simulation time:     4.08421s


---------------------------------------------------
OCP_FIM_10

.\OpenCAEPoro.exe ..\..\data\spe1a\spe1a.data  FIM  1  10  0.1

precond_type             = 69

Final time:          3655.500 Days
Total time steps:    420
Total Newton steps:  481 (+1 wasted steps)
Total linear steps:  6301 (+18 wasted steps)
Linear solve time:   0.460s (62.080%)
Simulation time:     0.741s


```



### SPE1B

![spe1b](D:\Lsz\OpenCAEPoro\data\figure\spe1b.png)

```
---------------------------------------------------
PS_IMPES_1

0.1  1  0.1

Final time:          3655.5 Days
Total time steps:    4168
Total nonlin steps:  4168
Total linear steps:  7601
Linear solver time:  0.824397s
Simulation time:     1.99241s

---------------------------------------------------
OCP_IMPEC_1

 .\OpenCAEPoro.exe ..\..\data\spe1b\spe1b.data  IMPEC  0.1  1  0.1

Final time:          3655.500 Days
Total time steps:    3849
Total Newton steps:  3849 (+9 wasted steps)
Total linear steps:  7198 (+15 wasted steps)
Linear solve time:   0.654s (60.250%)
Simulation time:     1.086s

---------------------------------------------------
PS_FIM_10

1  10  0.1

Final time:          3655.5 Days
Total time steps:    420
Total Newton steps:  437
Wasted Newton steps: 1
Total linear steps:  1212
Linear solver time:  2.6363s
Simulation time:     3.08908s


---------------------------------------------------
OCP_FIM_10

.\OpenCAEPoro.exe ..\..\data\spe1b\spe1b.data  FIM  1  10  0.1 

precond_type             = 69

Final time:          3655.500 Days
Total time steps:    421
Total Newton steps:  507 (+1 wasted steps)
Total linear steps:  6361 (+19 wasted steps)
Linear solve time:   0.467s (61.599%)
Simulation time:     0.759s
```



### SPE9

![spe9](D:\Lsz\OpenCAEPoro\data\figure\spe9.png)

```
---------------------------------------------------
PS_IMEPS_1

0.1  1  0.1

Final time:          900 Days
Total time steps:    7767
Total nonlin steps:  7767
Total linear steps:  12465
Linear solver time:  57.9099s
Simulation time:     276.62s

---------------------------------------------------
OCP_IMPEC_1

.\OpenCAEPoro.exe ..\..\data\spe9\spe9_IMPEC.data

Final time:          900.000 Days
Total time steps:    7486
Total Newton steps:  7486 (+43 wasted steps)
Total linear steps:  9897 (+96 wasted steps)
Linear solve time:   55.544s (61.784%)
Simulation time:     89.900s
---------------------------------------------------
PS_FIM_10

1  10  0.1

Final time:          900 Days
Total time steps:    157
Total Newton steps:  216
Wasted Newton steps: 53
Total linear steps:  1061
Linear solver time:  12.3944s
Simulation time:     27.515s


---------------------------------------------------
OCP_FIM_10

.\OpenCAEPoro.exe ..\..\data\spe9\spe9_FIM.data

precond_type             = 69

Final time:          900.000 Days
Total time steps:    162
Total Newton steps:  227 (+60 wasted steps)
Total linear steps:  3284 (+824 wasted steps)
Linear solve time:   8.752s (55.793%)
Simulation time:     15.686s
```



### SPE10

![spe10](D:\Lsz\OpenCAEPoro\data\figure\spe10.png)

```
PS_FIM

Final time:          2000 Days
Total time steps:    53
Total Newton steps:  215
Wasted Newton steps: 5
Total linear steps:  1305
Linear solver time:  1413.86s
Simulation time:     2905.24s


OCP_FIM

.\OpenCAEPoro.exe ..\..\data\spe10\spe10.data

Final time:          2000.000 Days
Total time steps:    53
Total Newton steps:  215 (+5 wasted steps)
Total linear steps:  1278 (+40 wasted steps)
Linear solve time:   1403.826s (74.590%)
Simulation time:     1882.050s

```



### CP

![CP](D:\Lsz\OpenCAEPoro\data\figure\CP.png)



```
---------------------------------------------------
PS_IMPEC_1

0.1  1  0.1

Final time:          1000.000 Days
Total time steps:    4541
Total nonlin steps:  4541
Total linear steps:  4543
Linear solver time:  2.269s
Simulation time:     28.627s

---------------------------------------------------
OCP_IMPEC_1

 .\OpenCAEPoro.exe ..\..\data\cornerpoint\CP.data  IMPEC  0.1  1  0.1

Final time:          1000.000 Days
Total time steps:    1999
Total Newton steps:  1999 (+3 wasted steps)
Total linear steps:  13516 (+35 wasted steps)
Linear solve time:   0.045s (38.852%)
Simulation time:     0.115s

---------------------------------------------------
PS_FIM_10

1  10  0.1

Final time:          1000.000 Days
Total time steps:    145
Total Newton steps:  256
Wasted Newton steps: 4
Total linear steps:  263
Linear solver time:  0.745s
Simulation time:     1.654s

---------------------------------------------------
OCP_FIM_10

.\OpenCAEPoro.exe ..\..\data\cornerpoint\CP.data  FIM 1  10  0.1
 
precond_type             = 69

Final time:          1000.000 Days
Total time steps:    143
Total Newton steps:  245 (+7 wasted steps)
Total linear steps:  2293 (+68 wasted steps)
Linear solve time:   0.037s (43.154%)
Simulation time:     0.086s

```



### SPE5

![spe5](D:\Lsz\OpenCAEPoro\data\figure\spe5.png)

饱和度差异对比 (OCP_FIM)

![FIM10_0.1](D:\Lsz\Matlab_practice\PennSim\ReadPRT\SPE5\FIM10_0.1.png)

![FIM10_0.05](D:\Lsz\Matlab_practice\PennSim\ReadPRT\SPE5\FIM10_0.05.png)



饱和度差异对比 (OCP_IMPEC)

![IMPEC1_0.1](D:\Lsz\Matlab_practice\PennSim\ReadPRT\SPE5\IMPEC1_0.1.png)

![IMPEC1_0.05](D:\Lsz\Matlab_practice\PennSim\ReadPRT\SPE5\IMPEC1_0.05.png)

```

OCP 的 Flash 部分存在效率问题
PS 数据图与下列数据有些差异，但结果几无二致。

------------------------------------------------------
PS_IMPES

TMARCH
0.1  1  0.1  5  0.3  0.3  300  0.2  0.3  0.001 

Final time:          7305 Days
Total time steps:    14228
Total nonlin steps:  14228
Total linear steps:  20914
flux time(MassConerve):  0.246282s
flux time(updateproperty):  0.133042s
MassConserve time:  0.596821s
updateProperty time:  47.4394s
Linear solver time:  1.4444s
Simulation time:     51.0602s

Method     iters     maxIt    tol
SSMSTA   49846316     300    1E-12  1E-2
SSMSP    7802131      100    1E-5
NRSP     3091611      55     1E-6

------------------------------------------------------
OCP_IMPEC

.\OpenCAEPoro.exe ..\..\data\spe5\spe5.data  IMPEC  0.1  1  0.1

SSMSTA:     6511743
NRSTA:      20659
SSMSP:      9579276
NRSP:       1741235
=========================================
Final time:          7305.000 Days
Total time steps:    10923
Total Newton steps:  10923 (+80 wasted steps)
Total linear steps:  17386 (+156 wasted steps)
Linear solve time:   1.153s (7.288%)
Simulation time:     15.818s

------------------------------------------------------
PS_FIM

1  10  0.01

Final time:          7305 Days
Total time steps:    21060
Total Newton steps:  23072
Wasted Newton steps: 13992
Total linear steps:  55425
Linear solver time:  74.4161s
Simulation time:     288.01s

Method     iters     maxIt    tol
SSMSTA   312298023     100    1E-12
SSMSP    15015547     100    1E-6
NRSP     5624310      55     1E-12

------------------------------------------------------
OCP_FIM

.\OpenCAEPoro.exe ..\..\data\spe5\spe5.data  FIM  1  10  0.1

precond_type             = 64

SSMSTA:     1156364
NRSTA:      7900
SSMSP:      1710655
NRSP:       242509
=========================================
Final time:          7305.000 Days
Total time steps:    878
Total Newton steps:  1633 (+26 wasted steps)
Total linear steps:  4970 (+58 wasted steps)
Linear solve time:   13.091s (74.471%)
Simulation time:     17.579s

```



### SPE3

![compare](D:\Lsz\Matlab_practice\PennSim\ReadSummary\20220129\spe3\compare.png)

饱和度误差均在 5% 以内

```
------------------------------------------------------
PS_IMPES

0.1  1  0.1

Final time:          3650 Days
Total time steps:    3670
Total nonlin steps:  3670
Total linear steps:  7345
flux time(MassConerve):  0.179398s
flux time(updateproperty):  0.175264s
MassConserve time:  0.485313s
updateProperty time:  164.133s
Linear solver time:  1.25312s
Simulation time:     176.526s

Method     iters     maxIt    tol
SSMSTA   102544160    300    1E-12  1E-2
SSMSP    22778402     100    1E-5
NRSP     6775525      55     1E-12

------------------------------------------------------
OCP_IMPEC

.\OpenCAEPoro.exe ..\..\data\spe5\spe5.data  IMPEC  0.1  1  0.1


SSMSTA:     23835193
NRSTA:      19507
SSMSP:      24063413
NRSP:       1468770
=========================================
Final time:          3650.000 Days
Total time steps:    3670
Total Newton steps:  3670 (+0 wasted steps)
Total linear steps:  7318 (+0 wasted steps)
Linear solve time:   1.037s (2.543%)
Simulation time:     40.777s

------------------------------------------------------
PS_FIM

1  10  0.1

Final time:          3650 Days
Total time steps:    390
Total Newton steps:  588
Total linear steps:  941
Linear solver time:  5.61642s
Simulation time:     39.5691s

Method     iters     maxIt    tol
SSMSTA   17690167   300    1E-12  1E-2
SSMSP    3674930     100    1E-5
NRSP     1222903      55     1E-12

------------------------------------------------------
OCP_FIM

.\OpenCAEPoro.exe ..\..\data\spe3\spe3.data  FIM  1  10  0.1

precond_type             = 64

SSMSTA:     3984873
NRSTA:      10756
SSMSP:      3929150
NRSP:       238198
=========================================
Final time:          3650.000 Days
Total time steps:    390
Total Newton steps:  559 (+0 wasted steps)
Total linear steps:  747 (+0 wasted steps)
Linear solve time:   4.464s (30.473%)
Simulation time:     14.649s
```

