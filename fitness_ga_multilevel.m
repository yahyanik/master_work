function fun = fitness_ga_multilevel(x)
global V_dc
global v11_dc
global v22_dc
global v33_dc
global M_ma
global v_dc_optim
global step1
global step2
global v1_ma
% n1 = round(x(8));
% n2 = round(x(9));
% n3 = round(x(10));
M = 0;v2 = 0;v3 = 0;v4 = 0;v5 = 0;v6 = 0;v7 = 0;
for ii = 1:1:step1
M = M +v11_dc* ((-1)^(ii-1))*cos(x(ii));
v2 = v2+v11_dc* ((-1)^(ii-1))*cos(5*x(ii));
v3 = v3+ v11_dc*((-1)^(ii-1))*cos(7*x(ii));
v4 = v4+v11_dc*((-1)^(ii-1))*cos(11*x(ii));
v5 = v5+v11_dc*((-1)^(ii-1))*cos(13*x(ii));
v6 = v6+v11_dc*((-1)^(ii-1))*cos(17*x(ii));
v7 = v7+v11_dc*((-1)^(ii-1))*cos(19*x(ii));
end
for ii = (step1+1):1:step2
    M = M+ v22_dc*(-1)^(ii)*cos(x(ii));
    v2 = v2+ v22_dc*(-1)^(ii)*cos(5*x(ii));
    v3 = v3+ v22_dc*(-1)^(ii)*cos(7*x(ii));
    v4 = v4+ v22_dc*(-1)^(ii)*cos(11*x(ii));
    v5 = v5+v22_dc*(-1)^(ii)*cos(13*x(ii));
    v6 = v6+ v22_dc*(-1)^(ii)*cos(17*x(ii));
    v7 = v7+v22_dc* (-1)^(ii)*cos(19*x(ii));
end
for ii =  (step2+1):1:7
    M = M+ v33_dc*(-1)^(ii-1)*cos(x(ii));
    v2 = v2+ v33_dc*(-1)^(ii-1)*cos(5*x(ii));
    v3 = v3+v33_dc* (-1)^(ii-1)*cos(7*x(ii));
    v4 = v4+ v33_dc*(-1)^(ii-1)*cos(11*x(ii));
    v5 = v5+ v33_dc*(-1)^(ii-1)*cos(13*x(ii));
    v6 = v6+ v33_dc*(-1)^(ii-1)*cos(17*x(ii));
    v7 = v7+v33_dc* (-1)^(ii-1)*cos(19*x(ii));
end
    M = M/3;
    v2 = (4.*v2).*(v_dc_optim).*(1/(5*pi));
    v3 = (4.*v3).*(v_dc_optim).*(1/(7*pi));
    v4 = (4.*v4).*(v_dc_optim).*(1/(11*pi));
    v5 = (4.*v5).*(v_dc_optim).*(1/(13*pi));
    v6 = (4.*v6).*(v_dc_optim).*(1/(17*pi));
    v7 = (4.*v7).*(v_dc_optim).*(1/(19*pi));
v1 = M*(V_dc*4/pi);
v1_ma = M_ma*(V_dc*4/pi);

KK = (100*((v1_ma-v1)/v1_ma))^4 + (1/5)*((50*v2/v1))^2 +...
    (1/7)*((50*v3/v1))^2 + (1/11)*((50*v4/v1))^2+(1/13)*...
    ((50*v5/v1))^2+(1/17)*((50*v6/v1))^2 +(1/19)*((50*v7/v1))^2;
if isnan(KK)
    fun = 100;
else
    fun = KK;
end
% fun = 100*abs((v1_ma-v1))+10*abs((v2))+10*abs((v3))++10*abs((v4))+...
%     10*abs((v5))+10*abs((v6))+10*abs((v7));
end

