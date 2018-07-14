function y = fitness_hybrid_newtonraphson(x)

global V_dc
global v11_dc
global v22_dc
global v33_dc
global M_ma
global step1
global step2
global v1_ma

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
    v2 = (4.*v2).*(V_dc/3).*(1/(5*pi));
    v3 = (4.*v3).*(V_dc/3).*(1/(7*pi));
    v4 = (4.*v4).*(V_dc/3).*(1/(11*pi));
    v5 = (4.*v5).*(V_dc/3).*(1/(13*pi));
    v6 = (4.*v6).*(V_dc/3).*(1/(17*pi));
    v7 = (4.*v7).*(V_dc/3).*(1/(19*pi));
v1 = M*(V_dc*4/pi);
v1_ma = M_ma*(V_dc*4/pi);

y(1) = (v1_ma-v1);
y(2) = v2;
y(3) = v3;
y(4) = v4;
y(5) = v5;
y(6) = v6;
y(7) = v7;
end




