### SPE1A

![compare](D:\Lsz\Matlab_practice\PennSim\ReadSummary\0928\spe1a\compare.png)

```
-----------------------------------------
PS_IMPES_1

TMARCH
0.1 1 0.1 5 0.3 0.3 300 0.2 0.3 0.001 


Final time:          3655.5 Days
Total time steps:    5533
Total nonlin steps:  5533
Total linear steps:  8828
Linear solver time:  0.932773s
Simulation time:     2.44659s

-----------------------------------------
OC_IMPES_1

0.1  1  0.1

Final time:          3655.5 Days
Total linear steps:  19177
Linear solve time:   1.22111s
Total time steps:    4970
Simulation time:     1.68812s

-----------------------------------------
OC_IMPES_0.1

0.1  0.1  0.1

Final time:          3655.5 Days
Total linear steps:  109630
Linear solve time:   6.66679s
Total time steps:    36563
Simulation time:     10.0991s
```



### SPE1B

![compare](D:\Lsz\Matlab_practice\PennSim\ReadSummary\0928\spe1b\compare.png)

```
-----------------------------------------
IMPES

TMARCH
0.1 1 0.1 5 0.3 0.3 300 0.2 0.3 0.001 

Final time:          3655.5 Days
Total time steps:    4168
Total nonlin steps:  4168
Total linear steps:  7601
Linear solver time:  0.824397s
Simulation time:     1.99241s

-----------------------------------------
OC_IMPES_1

0.1  1  0.1

Final time:          3655.5 Days
Total linear steps:  15679
Linear solve time:   1.03353s
Total time steps:    3964
Simulation time:     1.40656s

-----------------------------------------
OC_IMPES_0.1

0.1  10.1  0.1

Final time:          3655.5 Days
Total linear steps:  109441
Linear solve time:   5.84241s
Total time steps:    36563
Simulation time:     9.20546s
```



### SPE9

![compare](D:\Lsz\Matlab_practice\PennSim\ReadSummary\0928\spe9\well_control\compare.png)

![compare](D:\Lsz\Matlab_practice\PennSim\ReadSummary\0928\spe9\no_well_control\compare.png)





```
-----------------------------------------
PS_IMPES_1

TMARCH
0.1 1 0.1 5 0.3 0.3 300 0.3 0.3 0.001

-------------------
Summary
-------------------
Final time:          900 Days
Total time steps:    7767
Total nonlin steps:  7767
Total linear steps:  12465
Linear solver time:  57.9099s
Simulation time:     276.62s


-----------------------------------------
OC_IMPES_1_wc


Final time:          900 Days
Total linear steps:  16808
Linear solve time:   48.808s
Total time steps:    4106
Simulation time:     63.845s


-----------------------------------------
OC_IMPES_0.1_wc      with well control

Final time:          900 Days
Total linear steps:  32056
Linear solve time:   83.3344s
Total time steps:    9000
Simulation time:     115.682s


-----------------------------------------
OC_IMPES_1_nwc    no well control
Final time:          900 Days
Total linear steps:  16455
Linear solve time:   48.4109s
Total time steps:    4033
Simulation time:     63.0615s


-----------------------------------------
OC_IMPES_0.1_nwc

Final time:          900 Days
Total linear steps:  31955
Linear solve time:   84.0852s
Total time steps:    9000
Simulation time:     116.621s
```

