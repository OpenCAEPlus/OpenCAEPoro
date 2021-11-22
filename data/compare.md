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
Total Newton steps:  4661
Wasted Newton steps: 51
Simulation time:     1.344s
Total linear steps:  7884
Linear solve time:   0.789s (58.670%)



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
Total Newton steps:  479
Wasted Newton steps: 1
Simulation time:     5.362s
Total linear steps:  1428
Linear solve time:   5.058s (94.319%)

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
Total Newton steps:  3858
Wasted Newton steps: 19
Simulation time:     1.140s
Total linear steps:  7250
Linear solve time:   0.683s (59.931%)


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
Total Newton steps:  507
Wasted Newton steps: 1
Simulation time:     3.801s
Total linear steps:  1507
Linear solve time:   3.485s (91.696%)
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
Total time steps:    6317
Total Newton steps:  6317
Wasted Newton steps: 42
Simulation time:     79.238s
Total linear steps:  8875
Linear solve time:   49.847s (62.909%)


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
Total time steps:    160
Total Newton steps:  222
Wasted Newton steps: 63
Simulation time:     19.938s
Total linear steps:  1112
Linear solve time:   13.269s (66.553%)
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
Total time steps:    1992
Total Newton steps:  1992
Wasted Newton steps: 3
Simulation time:     0.509s
Total linear steps:  1995
Linear solve time:   0.441s (86.586%)


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
Total Newton steps:  261
Wasted Newton steps: 11
Simulation time:     0.441s
Total linear steps:  272
Linear solve time:   0.385s (87.120%)



```

