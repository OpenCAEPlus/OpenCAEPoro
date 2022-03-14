### SPE1A

![compare](D:\Lsz\Matlab_practice\PennSim\ReadSummary\1111\spe1a\compare.png)

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
OCP_IMPES_1

0.1  1  0.1

Final time:          3655.500 Days
Total time steps:    4661
Total Newton steps:  4661 (+51 wasted steps)
Total linear steps:  7785 (+99 wasted steps)
Linear solve time:   0.794s (59.696%)
Simulation time:     1.331s

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

1  10  0.1

precond_type             = 64

Final time:          3655.500 Days
Total time steps:    420
Total Newton steps:  474 (+1 wasted steps)
Total linear steps:  1414 (+4 wasted steps)
Linear solve time:   4.950s (94.546%)
Simulation time:     5.235s

---------------------------------------------------
OCP_FIM_10

1  10  0.1

precond_type             = 69

Final time:          3655.500 Days
Total time steps:    420
Total Newton steps:  480 (+1 wasted steps)
Total linear steps:  6288 (+18 wasted steps)
Linear solve time:   0.456s (60.876%)
Simulation time:     0.749s


```



### SPE1B

![compare](D:\Lsz\Matlab_practice\PennSim\ReadSummary\1111\spe1b\compare.png)

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
OCP_IMPES_1

0.1  1  0.1

Final time:          3655.500 Days
Total time steps:    3858
Total Newton steps:  3858 (+19 wasted steps)
Total linear steps:  7215 (+35 wasted steps)
Linear solve time:   0.678s (60.894%)
Simulation time:     1.113s

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

1  10  0.1

precond_type             = 64

Final time:          3655.500 Days
Total time steps:    420
Total Newton steps:  502 (+1 wasted steps)
Total linear steps:  1495 (+4 wasted steps)
Linear solve time:   3.356s (91.629%)
Simulation time:     3.663s

---------------------------------------------------
OCP_FIM_10

1  10  0.1  

precond_type             = 69

Final time:          3655.500 Days
Total time steps:    421
Total Newton steps:  507 (+1 wasted steps)
Total linear steps:  6362 (+19 wasted steps)
Linear solve time:   0.469s (61.180%)
Simulation time:     0.767s
```



### SPE9

![compare](D:\Lsz\Matlab_practice\PennSim\ReadSummary\1111\spe9\compare.png)

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
OCP_IMPES_1

0.1  1  0.1
dPlim  dSlim   dNlim   dVerrlim
300     0.0005    0.3    0.001

Final time:          900.000 Days
Total time steps:    7492
Total Newton steps:  7492 (+42 wasted steps)
Total linear steps:  9903 (+95 wasted steps)
Linear solve time:   54.309s (62.283%)
Simulation time:     87.197s

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

1  10  0.1

precond_type             = 64

Final time:          900.000 Days
Total time steps:    161
Total Newton steps:  225 (+58 wasted steps)
Total linear steps:  854 (+247 wasted steps)
Linear solve time:   12.950s (67.510%)
Simulation time:     19.182s

---------------------------------------------------
OCP_FIM_10

1  10  0.1

precond_type             = 69

Final time:          900.000 Days
Total time steps:    162
Total Newton steps:  227 (+60 wasted steps)
Total linear steps:  3260 (+808 wasted steps)
Linear solve time:   8.466s (56.374%)
Simulation time:     15.017s
```



### SPE10

![compare](D:\Lsz\Matlab_practice\PennSim\ReadSummary\1111\spe10\compare.png)

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

Final time:          2000.000 Days
Total time steps:    53
Total Newton steps:  217 (+5 wasted steps)
Total linear steps:  1294 (+40 wasted steps)
Linear solve time:   1424.577s (75.943%)
Simulation time:     1875.844s


```



### CP

![compare](D:\Lsz\Matlab_practice\PennSim\ReadSummary\1111\CP\compare.png)



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

0.1  1  0.1

Final time:          1000.000 Days
Total time steps:    2003
Total Newton steps:  2003 (+3 wasted steps)
Total linear steps:  13544 (+35 wasted steps)
Linear solve time:   0.043s (41.612%)
Simulation time:     0.103s

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

1  10  0.1
precond_type             = 69

Final time:          1000.000 Days
Total time steps:    143
Total Newton steps:  245 (+6 wasted steps)
Total linear steps:  2290 (+57 wasted steps)
Linear solve time:   0.034s (42.420%)
Simulation time:     0.079s

```



### SPE5

![compare](D:\Lsz\Matlab_practice\PennSim\ReadSummary\20220129\spe5\compare.png)

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

0.1  1  0.1

SSMSTA:     6718352
NRSTA:      19695
SSMSP:      19786221
NRSP:       1463313
=========================================
Final time:          7305.000 Days
Total time steps:    11424
Total Newton steps:  11424 (+72 wasted steps)
Total linear steps:  17859 (+139 wasted steps)
Linear solve time:   1.234s (5.886%)
Simulation time:     20.969s

SSMSTA:     6511743
NRSTA:      20659
SSMSP:      9579276
NRSP:       1741235
=========================================
Final time:          7305.000 Days
Total time steps:    10923
Total Newton steps:  10923 (+80 wasted steps)
Total linear steps:  17386 (+156 wasted steps)
Linear solve time:   1.184s (7.258%)
Simulation time:     16.309s

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

1  10  0.1
precond_type             = 64

SSMSTA:     1187658
NRSTA:      8035
SSMSP:      3085834
NRSP:       218804
=========================================
Final time:          7305.000 Days
Total time steps:    875
Total Newton steps:  1703 (+30 wasted steps)
Total linear steps:  5191 (+70 wasted steps)
Linear solve time:   12.436s (71.856%)
Simulation time:     17.306s


SSMSTA:     1154156
NRSTA:      7649
SSMSP:      1708720
NRSP:       242712
=========================================
Final time:          7305.000 Days
Total time steps:    876
Total Newton steps:  1634 (+26 wasted steps)
Total linear steps:  4966 (+58 wasted steps)
Linear solve time:   11.929s (74.470%)
Simulation time:     16.018s

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

0.1  1  0.1

SSMSTA:     23868228
NRSTA:      19537
SSMSP:      92506242
NRSP:       1829508
=========================================
Final time:          3650.000 Days
Total time steps:    3673
Total Newton steps:  3673 (+0 wasted steps)
Total linear steps:  7321 (+0 wasted steps)
Linear solve time:   1.104s (1.121%)
Simulation time:     98.496s

SSMSTA:     23835193
NRSTA:      19507
SSMSP:      24063413
NRSP:       1468770
=========================================
Final time:          3650.000 Days
Total time steps:    3670
Total Newton steps:  3670 (+0 wasted steps)
Total linear steps:  7318 (+0 wasted steps)
Linear solve time:   1.066s (2.565%)
Simulation time:     41.568s

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

1  10  0.1
precond_type             = 64

SSMSTA:     3984257
NRSTA:      10759
SSMSP:      14585718
NRSP:       288417
=========================================
Final time:          3650.000 Days
Total time steps:    390
Total Newton steps:  558 (+0 wasted steps)
Total linear steps:  715 (+0 wasted steps)
Linear solve time:   3.696s (16.837%)
Simulation time:     21.954s


SSMSTA:     3984948
NRSTA:      10757
SSMSP:      3930901
NRSP:       238102
=========================================
Final time:          3650.000 Days
Total time steps:    390
Total Newton steps:  559 (+0 wasted steps)
Total linear steps:  747 (+0 wasted steps)
Linear solve time:   3.898s (29.123%)
Simulation time:     13.384s
```

