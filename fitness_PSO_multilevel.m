function fun = fitness_PSO_multilevel(x)
global V_dc
global v11_dc
global v22_dc
global v33_dc
global M_ma
global v1_ma
M = 0;v2 = 0;v3 = 0;v4 = 0;v5 = 0;v6 = 0;v7 = 0;
V_dc = 300;
% int x(8)
% int x(9)
% int x(10)
% n1 = round(x(12));
% n2 = round(x(13));
% n3 = round(x(14));
for ii = 1:1:1
M = M + ((-1)^(ii-1))*cos(x(ii));
v2 = v2+ ((-1)^(ii-1))*cos(5*x(ii));
v3 = v3+ ((-1)^(ii-1))*cos(7*x(ii));
v4 = v4+((-1)^(ii-1))*cos(11*x(ii));
v5 = v5+((-1)^(ii-1))*cos(13*x(ii));
v6 = v6+((-1)^(ii-1))*cos(17*x(ii));
v7 = v7+((-1)^(ii-1))*cos(19*x(ii));
end
for ii = 2:1:2
    M = M+ (-1)^(ii-4)*cos(x(ii));
    v2 = v2+ (-1)^(ii-4)*cos(5*x(ii));
    v3 = v3+ (-1)^(ii-4)*cos(7*x(ii));
    v4 = v4+ (-1)^(ii-4)*cos(11*x(ii));
    v5 = v5+ (-1)^(ii-4)*cos(13*x(ii));
    v6 = v6+ (-1)^(ii-4)*cos(17*x(ii));
    v7 = v7+ (-1)^(ii-4)*cos(19*x(ii));
end
for ii =  3:1:7
    M = M+ (-1)^(ii-6)*cos(x(ii));
    v2 = v2+ (-1)^(ii-6)*cos(5*x(ii));
    v3 = v3+ (-1)^(ii-6)*cos(7*x(ii));
    v4 = v4+ (-1)^(ii-6)*cos(11*x(ii));
    v5 = v5+ (-1)^(ii-6)*cos(13*x(ii));
    v6 = v6+ (-1)^(ii-6)*cos(17*x(ii));
    v7 = v7+ (-1)^(ii-6)*cos(19*x(ii));
end

v1 = M*(V_dc*4/pi);
v1_ma = M_ma*(V_dc*4/pi);
fun = (100*((v1_ma-v1)/v1_ma))^4 + (1/5)*((50*v2/v1))^2 +...
    (1/7)*((50*v3/v1))^2 + (1/11)*((50*v4/v1))^2+(1/13)*...
    ((50*v5/v1))^2+(1/17)*((50*v6/v1))^2 +(1/19)*((50*v7/v1))^2;
% fun = 100*abs(v1_ma-v1)+50*abs(v2)+50*abs(v3)+50*abs(v4)+50*abs(v5)+50*abs(v6)+50*abs(v7);
% fun = abs(v1_ma-v1)+v2+v3+v4+v5+v6+v7;
end

