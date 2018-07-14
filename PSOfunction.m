clc;
clear all;
step = 0.01;   %mezrabe 10 bashad
% m_return = input()
global V_dc;
global v11_dc
global v22_dc
global v33_dc
global M_ma
global v1_ma
s = 3;v_dc_optim = 100;k = [1,1,1];
v11_dc = k(1,1);v22_dc = k(1,2);v33_dc = k(1,3);V_dc = v_dc_optim*sum(k);
k = 1;
LB = [0,0,0];
UB = [pi/2,pi/2,pi/2];
% tetha = zeros(4,100);
% f = zeros(1,100);
% v1 = zeros(1,100);
% v2 = zeros(1,100);
% v3 = zeros(1,100);
% v4 = zeros(1,100);
% xt = zeros(1,100);
for i = 0.8:step:0.8
%     k = i*20;
M_ma = i;
% temp = 0;
% temp1 = 0;
% for j = 1 : 5
% [x,fval,exitflag,output] = particleswarm(fun,2,lb,ub,options)
    options = optimoptions(@fmincon);
    [tetha(:,(k)),f(:,(k))] = particleswarm(@fitness_PSO,3,LB,UB,options);
%     [f(:,(k)), tetha(:,(k)), Swarm, history] = PSO(psoOptions);
%     tetha(:,(k)) = sort(tetha(:,(k)),'descend');
%     if j == 1
%         temp1 = tetha(:,(k));
%         temp = f(:,(k));
%     end
%     if f(:,(k)) < temp
%         temp = f(:,(k));
%         temp1 = tetha(:,(k));
%     end
% end

% f(:,(k)) = temp;
% tetha(:,(k)) = temp1;
if f(:,(k)) <=0.01
    f(:,(k)) =0.01;
end
for ii=1 : 20
    v22(ii) = v11_dc*cos((6*ii-1)*tetha(1,k))+v22_dc*cos((6*ii-1)*tetha(2,k))+v33_dc*cos((6*ii-1)*tetha(3,k));
    V22(ii) = (4.*v22(ii)).*(V_dc).*(1/((6*ii-1)*pi));
    v33(ii) = v11_dc*cos((6*ii+1)*tetha(1,k))+v22_dc*cos((6*ii+1)*tetha(2,k))+v33_dc*cos((6*ii+1)*tetha(3,k));
    V33(ii) = (4.*v33(ii)).*(V_dc).*(1/((6*ii+1)*pi));
end
M = (cos(tetha(1,k))+cos(tetha(2,k))+cos(tetha(3,k)))/3;
v1 = M*(4*100*4/pi);
v1_persent(k) = (abs(v1)./v1_ma)*100;
v2 = cos(5*tetha(1,k))+cos(5*tetha(2,k))+cos(5*tetha(3,k));
v2_persent(k) = (abs(v2.*V_dc./(5*pi))./v1_ma)*100;
v3 = cos(7*tetha(1,k))+cos(7*tetha(2,k))+cos(7*tetha(3,k));
v3_persent(k) = (abs(v3.*V_dc./(7*pi))./v1_ma)*100;
THD1(k) = 100*((sqrt(sum(V22.^2)+sum(V33.^2)))/v1);
% v4 = cos(11*tetha(1,k))+cos(11*tetha(2,k))+cos(11*tetha(3,k));
% v4_persent(k) = (abs(v4.*V_dc./(11*pi))./v1_ma)*100;
xt(k) = i;
k = k+1;
end

% if tetha <= 10^-4
%     tetha = 0;
% end
% if f <= 10^-4
%     f = 0;
% end
% f
% tetha
%%

tetha(1,:) = (tetha(1,:)./pi).*180;
tetha(2,:) = (tetha(2,:)./pi).*180;
tetha(3,:) = (tetha(3,:)./pi).*180;
%%
% x(1)=1.36959618152425
% V_dc =100
% M_ma = 0.2
% x(2) = 1.56749586575157
% x(3)=1.54573607774858
% x(4)=0.961918381942918
% M = (cos(x(1))+cos(x(2))+cos(x(3))+cos(x(4)))/4;
% v1 = M*(4*V_dc*9/pi);
% v2 = cos(5*x(1))+cos(5*x(2))+cos(5*x(3))+cos(5*x(4));
% v3 = cos(7*x(1))+cos(7*x(2))+cos(7*x(3))+cos(7*x(4));
% v4 = cos(11*x(1))+cos(11*x(2))+cos(11*x(3))+cos(11*x(4));
% v1_ma = M_ma*(4*V_dc*9/pi);
% 
% fun = (100*((v1_ma-v1)/220))^4 + ((1/2)*(50*v2/v1))^2 +...
%     ((1/3)*(50*v3/v1))^2 + ((1/4)*(50*v4/v1))^2



%%

figure(1)
semilogy(xt,f)
xlabel ('Modulation Index');
ylabel('Objective function');
figure (2)
plot(xt,tetha(1,:),xt,tetha(2,:),'r',xt,tetha(3,:),'g');
legend('teta1','teta2','teta3');
xlabel ('Modulation Index');
ylabel('Optimized Swiching Angeles');
% subplot(2,2,1)
% plot(xt,tetha(1,:));
% legend('teta1')
% subplot(2,2,2)
% plot(xt,tetha(2,:),'r');
% legend('teta2')
% subplot(2,2,3)
% plot(xt,tetha(3,:),'g');
% legend('teta3')
% subplot(2,2,4)
% % plot(xt,tetha(4,:),'k');
% % legend('teta4')
figure(3)
plot(xt,v1_persent,xt,v2_persent,xt,v3_persent);
grid off
legend('V1','Vh5','Vh7','Vh11','Vh13','Vh17','Vh19')
xlabel ('Modulation Index');
ylabel('Harmonic Amplitude Percent');
figure(4)
plot(xt,THD1,'r',xt,DF21,'g');
ylim([0 100])
legend('THD','DF2');
xlabel ('Modulation Index');
%DF2
%f
% tetha
%v1_persent
%v2_persent
%v3_persent
%v4_persent
%v5_persent
%v6_persent
%v7_persent
THD1
