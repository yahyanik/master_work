function fun_xim = fitness_simulink(x)
global M_ma
V_dc = 100;
global v1_ma
global f
global v1
global v5
global v7
global v11
t1 = x(1);
t1 = (t1/360)*f;
t2 = x(2);
t2 = (t2/360)*f;
t3 = x(3);
t3 = (t3/360)*f;
t4 = x(4);
t4 = (t4/360)*f;
assignin('base', 'tetha11', [1,t1]);
assignin('base', 'tetha12', [1,t2]);
assignin('base', 'tetha13', [1,t3]);
assignin('base', 'tetha14', [1,t4]);
assignin('base', 'fs', [1,f]);
sim ('multilevel');
% g = v_phase.data;
A1 = A1.data;
v1 = abs(0.01/2*(sum(A1)));
A5 = A5.data;
v5 = abs(0.01/2*(sum(A5)));
A7 = A7.data;
v7 = abs(0.01/2*(sum(A7)));
A11 = A11.data;
v11 = abs(0.01/2*(sum(A11)));
v1_ma = M_ma*(4*V_dc*4/pi);

fun_xim = (100*((v1_ma-v1)/v1_ma))^4 + ((1/2)*(50*v5/v1))^2 +...
    ((1/3)*(50*v7/v1))^2 + ((1/4)*(50*v11/v1))^2;

end






