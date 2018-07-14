clc;
clear all;
warning('off')
global M_ma
global v1_ma
global f
global v1
global v5
global v7
global v11
step = 0.1;
k = 1;
f = 0.02;
LB = [0,0,0,0];
UB = [90,90,90,90];
for i = 0.1:step:1

M_ma = i;
temp = 0;
temp1 = 0;

% for j = 1 : 2
    [tetha(:,(k)),fit(:,(k))] = ga(@fitness_simulink,4,[],[],[],[],LB,UB);
    tetha(:,(k)) = sort(tetha(:,(k)),'descend');
%     if j == 1
%         temp1 = tetha(:,(k));
%         temp = fit(:,(k));
%     end
%     if fit(:,(k)) < temp
%         temp = fit(:,(k));
%         temp1 = tetha(:,(k));
%     end
% end

% fit(:,(k)) = temp;
% tetha(:,(k)) = temp1;
v1_persent(k) = (v1./v1_ma)*100;
v5_persent(k) = (v5./v1_ma)*100;
v7_persent(k) = (v7./v1_ma)*100;
v11_persent(k) = (v11./v1_ma)*100;
if fit(:,(k)) <=0.01
    fit(:,(k)) =0.01;
end
    xt(k) = i;

    k = k+1;
    
 end

%%

figure(1)
semilogy(xt,fit)
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
plot(xt,v5_persent,'r');
plot(xt,v7_persent,'g');
plot(xt,v11_persent,'k');
legend('V1','Vh5','Vh7','Vh11')

