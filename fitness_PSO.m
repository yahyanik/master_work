function fun = fitness_PSO(x)

global V_dc
global v11_dc
global v22_dc
global v33_dc
global M_ma
global v1_ma

M = v11_dc*(cos(x(1))+v22_dc*cos(x(2))+v33_dc*cos(x(3)))/3;
v1 = M*(V_dc*4/pi);
v2 = v11_dc*cos(5*x(1))+v33_dc*cos(5*x(2))+v22_dc*cos(5*x(3));
v3 = v11_dc*cos(7*x(1))+v33_dc*cos(7*x(2))+v22_dc*cos(7*x(3));

v1_ma = M_ma*(V_dc*4/pi);

fun = (100*((v1_ma-v1)/v1_ma))^4 + ((1/2)*(50*v2/v1))^2 +...
    ((1/3)*(50*v3/v1))^2;

end

