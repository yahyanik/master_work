function fun = meghdar_fun_hybrid(x)

global V_dc
global v11_dc
global v22_dc
global v33_dc
global M_ma
global step1
global step2
global v1_ma
fun = (100*((v1_ma-x(1))/v1_ma))^4 + (1/5)*((50*x(2)/x(1)))^2 +...
    (1/7)*((50*x(3)/x(1)))^2 + (1/11)*((50*x(4)/x(1)))^2+(1/13)*...
    ((50*x(5)/x(1)))^2+(1/17)*((50*x(6)/x(1)))^2 +(1/19)*((50*x(7)/x(1)))^2;
end