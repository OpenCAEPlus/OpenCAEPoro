### SPE1A

![compare](D:\Lsz\Matlab_practice\PennSim\ReadSummary\1111\spe1a\compare.png)

```
PS_IMPES_1

Final time:          3655.5 Days
Total time steps:    5533
Total nonlin steps:  5533
Total linear steps:  8828
Linear solver time:  0.932773s
Simulation time:     2.44659s


OCP_IMPES_1

Final time:          3655.500 Days
Total time steps:    4661
Total Newton steps:  4661 (+51 wasted steps)
Total linear steps:  7785 (+99 wasted steps)
Linear solve time:   0.786s (58.841%)
Simulation time:     1.337s



PS_FIM_10

Final time:          3655.5 Days
Total time steps:    420
Total Newton steps:  434
Wasted Newton steps: 1
Total linear steps:  1235
Linear solver time:  3.62965s
Simulation time:     4.08421s


OCP_FIM_10

Final time:          3655.500 Days
Total time steps:    420
Total Newton steps:  474 (+1 wasted steps)
Total linear steps:  1414 (+4 wasted steps)
Linear solve time:   5.004s (94.434%)
Simulation time:     5.299s

```



### SPE1B

![compare](D:\Lsz\Matlab_practice\PennSim\ReadSummary\1111\spe1b\compare.png)

```

PS_IMPES_1


Final time:          3655.5 Days
Total time steps:    4168
Total nonlin steps:  4168
Total linear steps:  7601
Linear solver time:  0.824397s
Simulation time:     1.99241s


OCP_IMPES_1

Final time:          3655.500 Days
Total time steps:    3858
Total Newton steps:  3858 (+19 wasted steps)
Total linear steps:  7215 (+35 wasted steps)
Linear solve time:   0.666s (59.798%)
Simulation time:     1.114s


PS_FIM_10

Final time:          3655.5 Days
Total time steps:    420
Total Newton steps:  437
Wasted Newton steps: 1
Total linear steps:  1212
Linear solver time:  2.6363s
Simulation time:     3.08908s


OCP_FIM_10

Final time:          3655.500 Days
Total time steps:    420
Total Newton steps:  502 (+1 wasted steps)
Total linear steps:  1495 (+4 wasted steps)
Linear solve time:   3.356s (91.629%)
Simulation time:     3.663s
```



### SPE9

![compare](D:\Lsz\Matlab_practice\PennSim\ReadSummary\1111\spe9\compare.png)

```
PS_IMEPS_1

Final time:          900 Days
Total time steps:    7767
Total nonlin steps:  7767
Total linear steps:  12465
Linear solver time:  57.9099s
Simulation time:     276.62s


OCP_IMPES_1

Final time:          900.000 Days
Total time steps:    7492
Total Newton steps:  7492
Wasted Newton steps: 42
Wasted linear steps: 95
Simulation time:     93.851s
Total linear steps:  9903
Linear solve time:   57.371s (61.130%)


PS_FIM_10

Final time:          900 Days
Total time steps:    157
Total Newton steps:  216
Wasted Newton steps: 53
Total linear steps:  1061
Linear solver time:  12.3944s
Simulation time:     27.515s


OCP_FIM_10

Final time:          900.000 Days
Total time steps:    161
Total Newton steps:  225 (+58 wasted steps)
Total linear steps:  854 (+247 wasted steps)
Linear solve time:   13.387s (67.102%)
Simulation time:     19.949s
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
Total Newton steps:  217
Wasted Newton steps: 5
Simulation time:     1937.163s
Total linear steps:  1334
Linear solve time:   1463.019s (75.524%)


Final time:          2000.000 Days
Total time steps:    53
Total Newton steps:  214
Wasted Newton steps: 5
Simulation time:     2015.620s
Total linear steps:  1314
Linear solve time:   1531.146s (75.964%)
```



### CP

![compare](D:\Lsz\Matlab_practice\PennSim\ReadSummary\1111\CP\compare.png)



```
网格规模较少，使用 Pardiso 直接法进行求解


PS_IMPEC_1

Final time:          1000.000 Days
Total time steps:    4541
Total nonlin steps:  4541
Total linear steps:  4543
Linear solver time:  2.269s
Simulation time:     28.627s


OCP_IMPEC_1

Final time:          1000.000 Days
Total time steps:    2003
Total Newton steps:  2003 (+3 wasted steps)
Total linear steps:  13544 (+35 wasted steps)
Linear solve time:   0.043s (39.982%)
Simulation time:     0.107s


PS_FIM_10

Final time:          1000.000 Days
Total time steps:    145
Total Newton steps:  256
Wasted Newton steps: 4
Total linear steps:  263
Linear solver time:  0.745s
Simulation time:     1.654s


OCP_FIM_10

Final time:          1000.000 Days
Total time steps:    143
Total Newton steps:  245 (+6 wasted steps)
Total linear steps:  2290 (+57 wasted steps)
Linear solve time:   0.037s (43.063%)
Simulation time:     0.086s



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

------------------------------------------------------
OCP_IMPEC

0.1  1  0.1

Final time:          7305.000 Days
Total time steps:    12460
Total Newton steps:  12460 (+51 wasted steps)
Total linear steps:  18938 (+100 wasted steps)
Linear solve time:   1.354s (1.849%)
Simulation time:     73.251s

------------------------------------------------------
PS_FIM

1  10  0.01

Final time:          7305 Days
Total time steps:    21060
Total Newton steps:  23072
Wasted Newton steps: 13992
Total linear steps:  55425
Linear solver time:  76.2139s
Simulation time:     313.632s

------------------------------------------------------
OCP_FIM

1  10  0.01

Final time:          7305.000 Days
Total time steps:    2688
Total Newton steps:  3740 (+1908 wasted steps)
Total linear steps:  7796 (+4190 wasted steps)
Linear solve time:   21.751s (33.347%)
Simulation time:     65.226s

// modify dXsdXp
Final time:          7305.000 Days
Total time steps:    2738
Total Newton steps:  3741 (+1758 wasted steps)
Total linear steps:  7674 (+3945 wasted steps)
Linear solve time:   21.719s (33.430%)
Simulation time:     64.969s

// modify dXsdXp and muP,mux
// time-consuming, no improvement
Final time:          7305.000 Days
Total time steps:    3314
Total Newton steps:  4466 (+1972 wasted steps)
Total linear steps:  8788 (+4545 wasted steps)
Linear solve time:   23.122s (29.526%)
Simulation time:     78.309s
```

