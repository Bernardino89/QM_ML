%CHK=dimer_0.chk 
%MEM=100GB 
%Nprocs=32
#P PBEQIDH/gen Extralinks=L608 INT=(ultrafine) SCF=tight integral=NoXCTest

0_dimer DH-SVPD --> H C N  Def2-SVPD --> Br

0 1
H 0.00000000 0.00000000 2.48162600
Br 0.00000000 0.00000000 1.06084200
N 0.00000000 0.00000000 -2.17495900
C 0.00000000 0.00000000 -3.33131000
H 0.00000000 0.00000000 -4.39852000

Br     0
S    6   1.00
113286.3877600              0.14283037779D-02
17009.6263030              0.10950417496D-01
3870.1842567              0.54421006604D-01
1093.0357227              0.19047907695
356.39721797             0.39024642737
123.12539643             0.30814432514
S    3   1.00
236.74084007            -0.11228065671
28.468661070            0.64775962312
11.883443722            0.44235575986
S    3   1.00
21.269633312           -0.22642576323
3.6129226841           0.73823712008
1.6626648969           0.42683868694
S    1   1.00
0.34823793232          1.0000000
S    1   1.00
0.13019031394          1.0000000
S    1   1.00
0.50966709763D-01      1.0000000
P    5   1.00
1560.2801881              0.87166669072D-02
368.47859205             0.66243637420D-01
117.22978849             0.26495610385
42.648909248            0.53839160587
16.087225096            0.36579387888
P    3   1.00
8.6352810058           0.34248787366
3.5613665502           0.57500678213
1.5292626609           0.24330394172
P    1   1.00
0.53064294848          1.0000000
P    1   1.00
0.15702758965          1.0000000
P    1   1.00
0.09253910997                  1.0000000
D    4   1.00
104.85518642             0.22650147581D-01
30.281143688            0.13455483230
10.651394267            0.36474454537
3.8699456233           0.49044587056
D    1   1.00
1.3240876762            .27137289040
D    1   1.00
0.3890000              1.0000000
D    1   1.00
0.16579741459                   1.0000000
****
H    0
s  3   1.00
13.0107010              0.19682158D-01
1.9622572              0.13796524
0.44453796             0.47831935
S   1   1.00
0.4617867850             1.0000000
P   1   1.00
0.8000000              1.0000000
P   1   1.00
0.0791340242          1.0000000
****
N     0
S   5   1.00
1712.8415853             -0.53934125305D-02
257.64812677            -0.40221581118D-01
58.458245853           -0.17931144990
16.198367905           -0.46376317823
5.0052600809          -0.44171422662
S   1   1.00
0.58731856571          1.0000000
S   1   1.00
0.18764592253          1.0000000
S   1   1.00
0.96171241529D-01      1.0000000
P   3   1.00
13.571470233           -0.40072398852D-01
2.9257372874          -0.21807045028
0.79927750754         -0.51294466049
P   1   1.00
0.1918861474           1.0000000
D   1   1.00
1.0000000              1.0000000
D   1   1.00
0.3095443259           1.0000000
****
C     0
S    5   1.00
1238.4016938              0.54568832082D-02
186.29004992             0.40638409211D-01
42.251176346            0.18025593888
11.676557932            0.46315121755
3.5930506482           0.44087173314
S    1   1.00
0.40245147363          1.0000000
S    1   1.00
0.13090182668          1.0000000
S    1   1.00
0.67053540256D-01      1.0000000
P    3   1.00
9.4680970621           0.38387871728D-01
2.0103545142           0.21117025112
0.54771004707          0.51328172114
P    1   1.00
0.1508036550           1.0000000
D    1   1.00
0.8000000              1.0000000
D    1   1.00
0.3229294790           1.0000000
****

-79

