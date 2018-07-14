function y = fitness_newtonraphson(x)

global V_dc
global M_ma
V_dc = 100;
global v1_ma
v1_ma = M_ma*(4*V_dc*4/pi);
M = (cos(x(1))+cos(x(2))+cos(x(3))+cos(x(4)))/4;
v1 = M*(4*V_dc*4/pi);
y(1) = v1_ma-v1;
y(2) = cos(5*x(1))+cos(5*x(2))+cos(5*x(3))+cos(5*x(4));
y(3) = cos(7*x(1))+cos(7*x(2))+cos(7*x(3))+cos(7*x(4));
y(4) = cos(11*x(1))+cos(11*x(2))+cos(11*x(3))+cos(11*x(4));




end



