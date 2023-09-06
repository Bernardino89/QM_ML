%chk=NH3dimer1Re.chk
%mem=200GB
%nproc=32
#P PBEQIDH/gen INT=ultrafine SCF=tight Extralinks=L608

NH3dimer1Re

0 1
n    -1.578718 -0.046611 0.000000
h    -2.158621 0.136396 -0.809565
h    -2.158621 0.136396 0.809565
h    -0.849471 0.658193 0.000000
n     1.578718 0.046611 0.000000
h     2.158621 -0.136396 -0.809565
h     0.849471 -0.658193 0.000000
h     2.158621 -0.136396 0.809565
   
H    0
s  3   1.00
     13.0107010              0.19682158D-01
      1.9622572              0.13796524
      0.44453796             0.47831935
S   1   1.00
      0.4617867850           1.0000000
P   1   1.00
      0.8000000              1.0000000
P   1   1.00
      0.0791340242           1.0000000
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
      ALPHA         1.0000000
D   1   1.00
      1.0000000              1.0000000
D   1   1.00
      BETA          1.0000000
****

-79
