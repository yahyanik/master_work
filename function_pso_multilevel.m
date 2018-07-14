clc;
clear all;
step1 = 1;
step2 = 2;
step = 0.05;   %mezrabe 10 bashad
% m_return = input()
global V_dc
global v11_dc
global v22_dc
global v33_dc
global M_ma
global v1_ma
s = 3;v_dc_optim = 100;k = (1-0.05+rand(1,s)./10);
v11_dc = k(1,1);v22_dc = k(1,2);v33_dc = k(1,3);V_dc = v_dc_optim*sum(k);
k = 1;
LB = [0,0,0,0,0,0,0];
UB = [pi/2,pi/2,pi/2,pi/2,pi/2,pi/2,pi/2];
% LB = [0,0,0,0,0,0,0,0,0,0,0,1,1,1];
% UB = [pi/2,pi/2,pi/2,pi/2,pi/2,pi/2,pi/2,pi/2,pi/2,pi/2,pi/2,7,7,7];
for i = 0.8:step:0.8
    M = 0;v2 = 0;v3 = 0;v4 = 0;v5 = 0;v6 = 0;v7 = 0;
    M_ma = i;
    options = optimoptions('particleswarm');
    [tetha(:,(k)),f(:,(k))] = particleswarm(@fitness_PSO_multilevel,7,LB,UB,options);
%     n = round([tetha(8,(k)),tetha(9,(k)),tetha(10,(k))]);
%     tetha(:,(k)) = sort(tetha(:,(k)),'descend');
for jj = 1 : 100
       v22(jj) = 0;v33(jj) = 0; 
    for ii = 1:1:step1
        v22(jj) =v22(jj)+ v11_dc*(-1)^(ii-1)*cos((6*jj-1)*tetha(ii,k));
        v33(jj) =v33(jj)+ v11_dc*(-1)^(ii-1)*cos((6*jj+1)*tetha(ii,k));
    end
    for ii = (step1+1):1:step2
        v22(jj) = v22(jj) + v22_dc*(-1)^(ii)*cos((6*jj-1)*tetha(ii,k));
        v33(jj) = v33(jj) + v22_dc*(-1)^(ii)*cos((6*jj+1)*tetha(ii,k));
    end
    for ii = (step2+1):1:7
        v22(jj) = v22(jj) + v33_dc*(-1)^(ii-1)*cos((6*jj-1)*tetha(ii,k));
        v33(jj) = v33(jj) + v33_dc*(-1)^(ii-1)*cos((6*jj+1)*tetha(ii,k));
    end
        V22(jj) = (4.*v22(jj)).*(V_dc/3).*(1/((6*jj-1)*pi));
        V33(jj) = (4.*v33(jj)).*(V_dc/3).*(1/((6*jj+1)*pi));
   end
    for ii = 1:1:1
    M = M +(-1)^(ii-1)*cos(tetha(ii,k));
    v2 = v2+ (-1)^(ii-1)*cos(5*tetha(ii,k));
    v3 = v3+ (-1)^(ii-1)*cos(7*tetha(ii,k));
    v4 = v4+(-1)^(ii-1)*cos(11*tetha(ii,k));
    v5 = v5+(-1)^(ii-1)*cos(13*tetha(ii,k));
    v6 = v6+(-1)^(ii-1)*cos(17*tetha(ii,k));
    v7 = v7+(-1)^(ii-1)*cos(19*tetha(ii,k));
    end
    for ii = 2:1:2
        M = M+(-1)^(ii-4)*cos(tetha(ii,k));
        v2 = v2+ (-1)^(ii-4)*cos(5*tetha(ii,k));
        v3 = v3+ (-1)^(ii-4)*cos(7*tetha(ii,k));
        v4 = v4+ (-1)^(ii-4)*cos(11*tetha(ii,k));
        v5 = v5+ (-1)^(ii-4)*cos(13*tetha(ii,k));
        v6 = v6+ (-1)^(ii-4)*cos(17*tetha(ii,k));
        v7 = v7+ (-1)^(ii-4)*cos(19*tetha(ii,k));
    end
    for ii = 3:1:7
        M = M+ (-1)^(ii-6)*cos(tetha(ii,k));
        v2 = v2+ (-1)^(ii-6)*cos(5*tetha(ii,k));
        v3 = v3+ (-1)^(ii-6)*cos(7*tetha(ii,k));
        v4 = v4+ (-1)^(ii-6)*cos(11*tetha(ii,k));
        v5 = v5+ (-1)^(ii-6)*cos(13*tetha(ii,k));
        v6 = v6+ (-1)^(ii-6)*cos(17*tetha(ii,k));
        v7 = v7+(-1)^(ii-6)*cos(19*tetha(ii,k));
    end
    v1 = M*(V_dc*4/pi);
    v1_persent(k) = (abs(v1)./v1_ma)*100;
    v2_persent(k) = (abs(v2.*V_dc./(5*pi))./V_dc)*100;
    v3_persent(k) = (abs(v3.*V_dc./(7*pi))./V_dc)*100;
    v4_persent(k) = (abs(v4.*V_dc./(11*pi))./V_dc)*100;
    v5_persent(k) = (abs(v5.*V_dc./(13*pi))./V_dc)*100;
    v6_persent(k) = (abs(v6.*V_dc./(17*pi))./V_dc)*100;
    v7_persent(k) = (abs(v7.*V_dc./(19*pi))./V_dc)*100;
THD(k) = 100*((sqrt(v2^2+v3^2+v4^2+v5^2+v6^2+v7^2))/v1);
THD1(k) = 100*((sqrt(sum(V22.^2)+sum(V33.^2)))/v1);
    xt(k) = i;
    k = k+1;
end

%%
tetha(1,:) = (tetha(1,:)./pi).*180;
tetha(2,:) = (tetha(2,:)./pi).*180;
tetha(3,:) = (tetha(3,:)./pi).*180;
tetha(4,:) = (tetha(4,:)./pi).*180;
tetha(5,:) = (tetha(5,:)./pi).*180;
tetha(6,:) = (tetha(6,:)./pi).*180;
tetha(7,:) = (tetha(7,:)./pi).*180;

%%
figure(1)
semilogy(xt,f)
figure (2)
subplot(3,3,1)
plot(xt,tetha(1,:));
legend('teta1')
subplot(3,3,2)
plot(xt,tetha(2,:),'r');
legend('teta2')
subplot(3,3,3)
plot(xt,tetha(3,:),'g');
legend('teta3')
subplot(3,3,4)
plot(xt,tetha(4,:),'b');
legend('teta4')
subplot(3,3,5)
plot(xt,tetha(5,:),'y');
legend('teta5')
subplot(3,3,6)
plot(xt,tetha(6,:),'m');
legend('teta6')
subplot(3,3,7)
plot(xt,tetha(7,:),'c');
legend('teta7')
figure(3)

plot(xt,v1_persent,xt,v2_persent,xt,v3_persent,xt,v4_persent,xt,v5_persent,xt,v6_persent,xt,v7_persent);
grid on
legend('V1','Vh5','Vh7','Vh11','Vh13','Vh17','Vh19')


%%
tetha = sort(tetha,'descend')
THD
THD1