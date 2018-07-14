clc;
clear all;
step = 0.1;   %mezrabe 10 bashad
% m_return = input()
global V_dc;
global M_ma
global v1_ma
k = 1;
% LB = [0,0,0,0];
% UB = [pi/2,pi/2,pi/2,pi/2];
% tetha = zeros(4,100);
% f = zeros(1,100);
% v1 = zeros(1,100);
% v2 = zeros(1,100);
% v3 = zeros(1,100);
% v4 = zeros(1,100);
% xt = zeros(1,100);
for i = step:step:1
%     k = i*20;
M_ma = i;
v1_ma = M_ma*(4*V_dc*4/pi);

% temp = 0;
% temp1 = 0;
% for j = 1 : 5
    x0 = (pi/4)*(rand(4,1))
    [te,f] = newtonraphson(@fitness_newtonraphson,[],x0,[],[]);
    [a,b] = size (te);
    tetha(:,(k)) = te(:,b);
    if tetha (1,(k)) >= (pi/4)
        tetha(1,(k)) = pi/4;
    end
    if tetha (2,(k)) >= (pi/4)
        tetha(2,(k)) = pi/4;
    end
    if tetha (3,(k)) >= (pi/4)
        tetha(3,(k)) = pi/4;
    end
    if tetha (4,(k)) >= (pi/4)
        tetha(4,(k)) = pi/4;
    end
    if tetha (1,(k)) <= 0
        tetha(1,(k)) = 0;
    end
    if tetha (2,(k)) <= 0
        tetha(2,(k)) = 0;
    end
    if tetha (3,(k)) <= 0
        tetha(3,(k)) = 0;
    end
    if tetha (4,(k)) <= 0
        tetha(4,(k)) = 0;
    end
    tetha(:,(k)) = sort(tetha(:,(k)),'descend');


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
% if f(:,(k)) <=0.01
%     f(:,(k)) =0.01;
% end
M = (cos(tetha(1,k))+cos(tetha(2,k))+cos(tetha(3,k))+cos(tetha(4,k)))/4;
v1 = M*(4*V_dc*4/pi);
v1_persent(k) = (abs(v1)./v1_ma)*100;
v2 = cos(5*tetha(1,k))+cos(5*tetha(2,k))+cos(5*tetha(3,k))+cos(5*tetha(4,k));
v2_persent(k) = (abs(v2.*4.*V_dc./(5*pi))./v1_ma)*100;
v3 = cos(7*tetha(1,k))+cos(7*tetha(2,k))+cos(7*tetha(3,k))+cos(7*tetha(4,k));
v3_persent(k) = (abs(v3.*4.*V_dc./(7*pi))./v1_ma)*100;
v4 = cos(11*tetha(1,k))+cos(11*tetha(2,k))+cos(11*tetha(3,k))+cos(11*tetha(4,k));
v4_persent(k) = (abs(v4.*4.*V_dc./(11*pi))./v1_ma)*100;
fun = (100*((v1_ma-v1)/v1_ma))^4 + ((1/2)*(50*v2/v1))^2 +...
    ((1/3)*(50*v3/v1))^2 + ((1/4)*(50*v4/v1))^2;
fitt(:,(k)) = fun;
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
tetha(4,:) = (tetha(4,:)./pi).*180;
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
semilogy(xt,fitt)
figure (2)
subplot(2,2,1)
plot(xt,tetha(1,:));
legend('teta1')
subplot(2,2,2)
plot(xt,tetha(2,:),'r');
legend('teta2')
subplot(2,2,3)
plot(xt,tetha(3,:),'g');
legend('teta3')
subplot(2,2,4)
plot(xt,tetha(4,:),'k');
legend('teta4')
figure(3)
hold on
plot(xt,v1_persent);
plot(xt,v2_persent,'r');
plot(xt,v3_persent,'g');
plot(xt,v4_persent,'k');
legend('V1','Vh5','Vh7','Vh11')

