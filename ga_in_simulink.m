clc 
clear all
f = 0.02;
t1 = 60.71;t1 = (t1/360)*f;
t2 = 37.88;t2 = (t2/360)*f;
t3 = 20.35;t3 = (t3/360)*f;
t4 = 10.34;t4 = (t4/360)*f;
tetha11 = [1,t1];tetha12 = [1,t2];tetha13 = [1,t3];tetha14 = [1,t4];
tetha11_120 = [1,(t1+0.006666666666667)];
tetha12_120 = [1,(t2+0.006666666666667)];
tetha11_240 = [1,(t3+2*0.006666666666667)];
tetha12_240 = [1,(t4+2*0.006666666666667)];
fs = [1,f];

sim ('untitled');

g = v_phase.data;
ts = t.data;
A1 = A1.data;
B1 = B1.data;
V1 = abs(0.01/2*(sum(A1)))
k2 = 0.01/2*(sum(B1));
A5 = A5.data;
B5 = B5.data;
V5 = abs(0.01/2*(sum(A5)))
k3 = 0.01/2*(sum(B5));
A7 = A7.data;
B7 = B7.data;
V7 = abs(0.01/2*(sum(A7)))
k4 = 0.01/2*(sum(B7));
A11 = A11.data;
B11 = B11.data;
V11 = abs(0.01/2*(sum(A11)))
k5 = 0.01/2*(sum(B11));