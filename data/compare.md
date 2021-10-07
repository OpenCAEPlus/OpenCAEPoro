### SPE1A

![compare](D:\Lsz\Matlab_practice\PennSim\ReadSummary\0928\spe1a\compare02.png)

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
OCP_IMPES_1

0.1  1  0.1

Final time:          3655.5 Days
Total linear steps:  20877
Linear solve time:   1.29674s
Total time steps:    5533
Simulation time:     1.80172s

-----------------------------------------
OCP_IMPES_0.1

0.1  0.1  0.1

Final time:          3655.5 Days
Total linear steps:  109617
Linear solve time:   6.53836s
Total time steps:    36563
Simulation time:     9.98476s
```



### SPE1B

![compare](D:\Lsz\Matlab_practice\PennSim\ReadSummary\0928\spe1b\compare02.png)

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
OCP_IMPES_1

0.1  1  0.1

Final time:          3655.5 Days
Total linear steps:  15348
Linear solve time:   0.977343s
Total time steps:    3847
Simulation time:     1.32726s

-----------------------------------------
OCP_IMPES_0.1

0.1  10.1  0.1

Final time:          3655.5 Days
Total linear steps:  109444
Linear solve time:   5.81064s
Total time steps:    36563
Simulation time:     9.15255s
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

![compare](D:\Lsz\Matlab_practice\PennSim\ReadSummary\0928\spe9\well_control_p3_Ve\compare.png)

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
OCP_IMPES_0.1       with well control and Volume error control

Final time:          900 Days
Total linear steps:  32151
Linear solve time:   82.9865s
Total time steps:    9005
Simulation time:     115.601s


-----------------------------------------
OCP_IMPES_1    with well control and Volume error control
Final time:          900 Days
Total linear steps:  16617
Linear solve time:   47.8777s
Total time steps:    3979
Simulation time:     62.0819s



```

