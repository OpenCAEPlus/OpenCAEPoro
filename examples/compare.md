# Numerical Results for Comparison

## SPE1A

![spe1a](figure\spe1a\spe1a.png)

```
---------------------------------------------------
OCP_IMPEC_1

.\OpenCAEPoro.exe ..\..\data\spe1a\spe1a.data  IMPEC  0.1  1  0.1

Final time:          3655.500 Days
Total time steps:    5271
Total Newton steps:  5271 (+39 wasted steps)
Total linear steps:  8327 (+73 wasted steps)
Linear solve time:   0.833s (56.367%)
Simulation time:     1.478s


---------------------------------------------------
OCP_FIM_10

.\OpenCAEPoro.exe ..\..\data\spe1a\spe1a.data  FIM  1  10  0.1

precond_type             = 69

Final time:          3655.500 Days
Total time steps:    420
Total Newton steps:  477 (+1 wasted steps)
Total linear steps:  6222 (+18 wasted steps)
Linear solve time:   0.458s (57.173%)
Simulation time:     0.800s
```



## SPE1B

![spe1b](figure\spe1b\spe1b.png)

```
---------------------------------------------------
OCP_IMPEC_1

 .\OpenCAEPoro.exe ..\..\data\spe1b\spe1b.data  IMPEC  0.1  1  0.1

Final time:          3655.500 Days
Total time steps:    3849
Total Newton steps:  3849 (+9 wasted steps)
Total linear steps:  7198 (+15 wasted steps)
Linear solve time:   0.667s (58.381%)
Simulation time:     1.143s


---------------------------------------------------
OCP_FIM_10

.\OpenCAEPoro.exe ..\..\data\spe1b\spe1b.data  FIM  1  10  0.1 

precond_type             = 69

Final time:          3655.500 Days
Total time steps:    422
Total Newton steps:  498 (+1 wasted steps)
Total linear steps:  6255 (+19 wasted steps)
Linear solve time:   0.466s (56.441%)
Simulation time:     0.826s
```


## SPE9

![spe9](figure\spe9\spe9.png)

```
---------------------------------------------------
OCP_IMPEC_1

.\OpenCAEPoro.exe ..\..\data\spe9\spe9_IMPEC.data

Final time:          900.000 Days
Total time steps:    7486
Total Newton steps:  7486 (+43 wasted steps)
Total linear steps:  9897 (+96 wasted steps)
Linear solve time:   57.486s (61.776%)
Simulation time:     93.056s


---------------------------------------------------
OCP_FIM_10

.\OpenCAEPoro.exe ..\..\data\spe9\spe9_FIM.data

precond_type             = 69

Final time:          900.000 Days
Total time steps:    165
Total Newton steps:  213 (+67 wasted steps)
Total linear steps:  3094 (+849 wasted steps)
Linear solve time:   8.896s (53.585%)
Simulation time:     16.601s
```



## SPE10

![spe10](figure\spe10\spe10.png)

```
---------------------------------------------------
PS_FIM

Final time:          2000 Days
Total time steps:    53
Total Newton steps:  215
Wasted Newton steps: 5
Total linear steps:  1305
Linear solver time:  1413.86s
Simulation time:     2905.24s

---------------------------------------------------
OCP_FIM

.\OpenCAEPoro.exe ..\..\data\spe10\spe10.data

Final time:          2000.000 Days
Total time steps:    53
Total Newton steps:  215 (+5 wasted steps)
Total linear steps:  1278 (+40 wasted steps)
Linear solve time:   1343.855s (73.547%)
Simulation time:     1827.200s
```



## CP

![CP](figure\Cornerpoint\CP.png)



```
---------------------------------------------------
OCP_IMPEC_1

 .\OpenCAEPoro.exe ..\..\data\cornerpoint\CP.data  IMPEC  0.1  1  0.1

Final time:          1000.000 Days
Total time steps:    1999
Total Newton steps:  1999 (+3 wasted steps)
Total linear steps:  13516 (+35 wasted steps)
Linear solve time:   0.043s (39.826%)
Simulation time:     0.109s


---------------------------------------------------
OCP_FIM_10

.\OpenCAEPoro.exe ..\..\data\cornerpoint\CP.data  FIM 1  10  0.1
 
precond_type             = 69

Final time:          1000.000 Days
Total time steps:    143
Total Newton steps:  240 (+6 wasted steps)
Total linear steps:  2234 (+57 wasted steps)
Linear solve time:   0.038s (30.290%)
Simulation time:     0.083s

```



## SPE5

![spe5](figure\spe5\spe5.png)


```
---------------------------------------------------
OCP_IMPEC

.\OpenCAEPoro.exe ..\..\data\spe5\spe5.data  IMPEC  0.1  1  0.1

SSMSTA:     6829509
NRSTA:      20699
SSMSP:      9985287
NRSP:       1798777
=========================================
Final time:          7305.000 Days
Total time steps:    11608
Total Newton steps:  11608 (+65 wasted steps)
Total linear steps:  18085 (+129 wasted steps)
Linear solve time:   1.239s (7.837%)
Simulation time:     15.815s


---------------------------------------------------
OCP_FIM

.\OpenCAEPoro.exe ..\..\data\spe5\spe5.data  FIM  1  10  0.1

precond_type             = 64

SSMSTA:     979435
NRSTA:      6790
SSMSP:      1491102
NRSP:       201120
=========================================
Final time:          7305.000 Days
Total time steps:    859
Total Newton steps:  1359 (+10 wasted steps)
Total linear steps:  3624 (+22 wasted steps)
Linear solve time:   8.888s (75.365%)
Simulation time:     11.793s

OCP_FIM

.\OpenCAEPoro.exe ..\..\data\spe5\spe5.data  FIM  1  50  0.1

precond_type             = 64

SSMSTA:     595065
NRSTA:      3226
SSMSP:      789150
NRSP:       119160
=========================================
Final time:          7305.000 Days
Total time steps:    336
Total Newton steps:  797 (+8 wasted steps)
Total linear steps:  2551 (+20 wasted steps)
Linear solve time:   6.803s (79.795%)
Simulation time:     8.526s

OCP_FIM

.\OpenCAEPoro.exe ..\..\data\spe5\spe5.data  FIM  1  100  0.1

precond_type             = 64

TUNING
-- Init     max    min   incre   chop    cut
   1       100     0.1      5    0.3    0.3                    /
--  dPlim  dSlim   dNlim   dVerrlim
     300     0.2    0.3    0.001                               /
-- itNRmax  NRtol  dPmax  dSmax  dPmin   dSmin   dVerrmax
       10    1E-6   200    0.2    0.1      1E-3    0.01          /
/

SSMSTA:     923073
NRSTA:      6930
SSMSP:      1061139
NRSP:       173833
=========================================
Final time:          7305.000 Days
Total time steps:    357
Total Newton steps:  1279 (+27 wasted steps)
Total linear steps:  4074 (+106 wasted steps)
Linear solve time:   10.307s (80.368%)
Simulation time:     12.825s

---------------------------------------------------
OCP_NEW_FIMn

.\OpenCAEPoro.exe ..\..\data\spe5\spe5.data  FIMn  1  10  0.1

precond_type             = 64

SSMSTA:     1005243
NRSTA:      5634
SSMSP:      1497622
NRSP:       218647
=========================================
Final time:          7305.000 Days
Total time steps:    853
Total Newton steps:  1464 (+5 wasted steps)
Total linear steps:  3862 (+7 wasted steps)
Linear solve time:   9.421s (73.430%)
Simulation time:     12.830s

.\OpenCAEPoro.exe ..\..\data\spe5\spe5.data  FIMn  1  50  0.1

precond_type             = 64

SSMSTA:     641364
NRSTA:      3560
SSMSP:      846803
NRSP:       130213
=========================================
Final time:          7305.000 Days
Total time steps:    351
Total Newton steps:  860 (+19 wasted steps)
Total linear steps:  2651 (+39 wasted steps)
Linear solve time:   6.927s (76.865%)
Simulation time:     9.011s

.\OpenCAEPoro.exe ..\..\data\spe5\spe5.data  FIMn  1  100  0.1

precond_type             = 64

TUNING
-- Init     max    min   incre   chop    cut
   1       100     0.1      5    0.3    0.3                    /
--  dPlim  dSlim   dNlim   dVerrlim
     300     0.2    0.3    0.001                               /
-- itNRmax  NRtol  dPmax  dSmax  dPmin   dSmin   dVerrmax
       10    1E-6   200    0.2    0.1      1E-3    0.01          /
/


SSMSTA:     1042349
NRSTA:      11796
SSMSP:      1226341
NRSP:       185670
=========================================
Final time:          7305.000 Days
Total time steps:    365
Total Newton steps:  1323 (+51 wasted steps)
Total linear steps:  4076 (+214 wasted steps)
Linear solve time:   10.651s (77.696%)
Simulation time:     13.709s

```



## SPE3

![spe3](figure/spe3/spe3.png)

```
---------------------------------------------------
OCP_IMPEC

.\OpenCAEPoro.exe ..\..\data\spe5\spe5.data  IMPEC  0.1  1  0.1

SSMSTA:     23835204
NRSTA:      19506
SSMSP:      24063413
NRSP:       1468737
=========================================
Final time:          3650.000 Days
Total time steps:    3670
Total Newton steps:  3670 (+0 wasted steps)
Total linear steps:  7318 (+0 wasted steps)
Linear solve time:   1.078s (2.714%)
Simulation time:     39.726s

---------------------------------------------------
OCP_FIM

.\OpenCAEPoro.exe ..\..\data\spe3\spe3.data  FIM  1  10  0.1

precond_type             = 64

SSMSTA:     3891265
NRSTA:      10763
SSMSP:      3784500
NRSP:       229583
=========================================
Final time:          3650.000 Days
Total time steps:    390
Total Newton steps:  542 (+0 wasted steps)
Total linear steps:  719 (+0 wasted steps)
Linear solve time:   4.025s (33.470%)
Simulation time:     12.026s

.\OpenCAEPoro.exe ..\..\data\spe3\spe3.data  FIM  1  50  0.1

precond_type             = 64

SSMSTA:     1867933
NRSTA:      9174
SSMSP:      1665461
NRSP:       96855
=========================================
Final time:          3650.000 Days
Total time steps:    131
Total Newton steps:  233 (+2 wasted steps)
Total linear steps:  351 (+3 wasted steps)
Linear solve time:   2.286s (38.824%)
Simulation time:     5.888s


.\OpenCAEPoro.exe ..\..\data\spe3\spe3.data  FIM  1  100  0.1

precond_type             = 64

TUNING
-- Init     max    min   incre   chop    cut
   1       100    0.1      5    0.3    0.3                    /
--  dPlim  dSlim   dNlim   dVerrlim
     300     1    0.3    0.001                               /
-- itNRmax  NRtol  dPmax  dSmax  dPmin   dSmin   dVerrmax
       10    1E-6   200    0.2    0.1    1E-3    0.01          /
/

SSMSTA:     2422183
NRSTA:      9080
SSMSP:      1666641
NRSP:       124831
=========================================
Final time:          3650.000 Days
Total time steps:    107
Total Newton steps:  336 (+26 wasted steps)
Total linear steps:  551 (+32 wasted steps)
Linear solve time:   3.258s (42.629%)
Simulation time:     7.643s

---------------------------------------------------

OCP_NEW_FIM

.\OpenCAEPoro.exe ..\..\data\spe3\spe3.data  FIMn  1  10  0.1

precond_type             = 64

SSMSTA:     3879409
NRSTA:      11023
SSMSP:      3736899
NRSP:       224124
=========================================
Final time:          3650.000 Days
Total time steps:    390
Total Newton steps:  529 (+0 wasted steps)
Total linear steps:  732 (+0 wasted steps)
Linear solve time:   4.012s (31.939%)
Simulation time:     12.560s

.\OpenCAEPoro.exe ..\..\data\spe3\spe3.data  FIMn  1  50  0.1

precond_type             = 64

SSMSTA:     1642488
NRSTA:      9156
SSMSP:      1609287
NRSP:       85851
=========================================
Final time:          3650.000 Days
Total time steps:    131
Total Newton steps:  199 (+1 wasted steps)
Total linear steps:  315 (+2 wasted steps)
Linear solve time:   2.099s (37.360%)
Simulation time:     5.618s


.\OpenCAEPoro.exe ..\..\data\spe3\spe3.data  FIMn  1  100  0.1

precond_type             = 64

TUNING
-- Init     max    min   incre   chop    cut
   1       100    0.1      5    0.3    0.3                    /
--  dPlim  dSlim   dNlim   dVerrlim
     300     1    0.3    0.001                               /
-- itNRmax  NRtol  dPmax  dSmax  dPmin   dSmin   dVerrmax
       10    1E-6   200    0.2    0.1    1E-3    0.01          /
/

SSMSTA:     1805051
NRSTA:      9007
SSMSP:      1489732
NRSP:       91106
=========================================
Final time:          3650.000 Days
Total time steps:    104
Total Newton steps:  241 (+4 wasted steps)
Total linear steps:  426 (+6 wasted steps)
Linear solve time:   2.528s (40.577%)
Simulation time:     6.229s
```

## SPE5r
![spe5r](figure/spe5refine/spe5refine.png)
```
---------------------------------------------------
OCP_FIM

.\OpenCAEPoro.exe ..\..\data\spe5refine\spe5-70x70x30-2y.data  FIM  1  10  0.1

precond_type   =  64

SSMSTA:     204184086
NRSTA:      83
SSMSP:      140490334
NRSP:       50132302
=========================================
Final time:          730.000 Days
Total time steps:    84
Total Newton steps:  263 (+29 wasted steps)
Total linear steps:  897 (+75 wasted steps)
Linear solve time:   602.064s (46.405%)
Simulation time:     1297.405s

.\OpenCAEPoro.exe ..\..\data\spe5refine\spe5-70x70x30-2y.data  FIM  1  20  0.1

precond_type   =  64

SSMSTA:     201297580
NRSTA:      72
SSMSP:      116490470
NRSP:       41681486
=========================================
Final time:          730.000 Days
Total time steps:    51
Total Newton steps:  223 (+48 wasted steps)
Total linear steps:  807 (+127 wasted steps)
Linear solve time:   584.719s (48.015%)
Simulation time:     1217.790s


.\OpenCAEPoro.exe ..\..\data\spe5refine\spe5-70x70x30-2y.data  FIM  1  50  0.1

precond_type   =  64

SSMSTA:     190604401
NRSTA:      117
SSMSP:      149990315
NRSP:       56039151
=========================================
Final time:          730.000 Days
Total time steps:    47
Total Newton steps:  213 (+85 wasted steps)
Total linear steps:  781 (+442 wasted steps)
Linear solve time:   687.898s (48.709%)
Simulation time:     1412.269s

---------------------------------------------------
OCP_NEW_FIM

.\OpenCAEPoro.exe ..\..\data\spe5refine\spe5-70x70x30-2y.data  FIMn  1  10  0.1

precond_type   =  64

SSMSTA:     180930118
NRSTA:      80
SSMSP:      140805797
NRSP:       50478353
=========================================
Final time:          730.000 Days
Total time steps:    80
Total Newton steps:  252 (+14 wasted steps)
Total linear steps:  864 (+24 wasted steps)
Linear solve time:   552.229s (42.005%)
Simulation time:     1314.664s

.\OpenCAEPoro.exe ..\..\data\spe5refine\spe5-70x70x30-2y.data  FIMn  1  20  0.1

precond_type   =  64

SSMSTA:     178769093
NRSTA:      72
SSMSP:      124291924
NRSP:       45607746
=========================================
Final time:          730.000 Days
Total time steps:    50
Total Newton steps:  218 (+34 wasted steps)
Total linear steps:  807 (+112 wasted steps)
Linear solve time:   545.601s (43.609%)
Simulation time:     1251.113s

.\OpenCAEPoro.exe ..\..\data\spe5refine\spe5-70x70x30-2y.data  FIMn  1  50  0.1

precond_type   =  64

SSMSTA:     193803346
NRSTA:      86
SSMSP:      160981671
NRSP:       60714439
=========================================
Final time:          730.000 Days
Total time steps:    50
Total Newton steps:  210 (+93 wasted steps)
Total linear steps:  757 (+397 wasted steps)
Linear solve time:   665.262s (43.112%)
Simulation time:     1543.095s

```