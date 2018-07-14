
clear all;
format long
core = 1;
global step1
global step2
step1 = 1;
step2 = 2;
% for core = 1:6
% %     clear all;
%     if core == 1
%         step1 = 5;
%         step2 = 6;
%     elseif core == 2
%         step1 = 1;
%         step2 = 6;
%     elseif core == 3
%         step1 = 1;
%         step2 = 2;
%     elseif core == 4
%         step1 = 3;
%         step2 = 6;
%     elseif core == 5
%         step1 = 3;
%         step2 = 4;
%     elseif core == 6
%         step1 = 1;
%         step2 = 4;
%     end
    %%

% parpool('local')
step = 0.05;   %mezrabe 10 bashad
% m_return = input()
global V_dc
global v11_dc
global v22_dc
global v33_dc
global M_ma
global v1_ma
% V_dc = 300;

v_dc_optim = 100;
V_dc = 300;
v11_dc = 1;
v22_dc = 1;
v33_dc = 1;k = 1;
LB = [0,0,0,0,0,0,0];
UB = [pi/2,pi/2,pi/2,pi/2,pi/2,pi/2,pi/2];
% LB = [0,0,0,0,0,0,0,0,0,0,0,1,1,1];
% UB = [pi/2,pi/2,pi/2,pi/2,pi/2,pi/2,pi/2,pi/2,pi/2,pi/2,pi/2,7,7,7];
for i = 0.8:step:0.8
     temp = 0;
     temp1 = 0;
%      for j =1:1:10
        M=0;v2=0;v3=0;v4=0;v5=0;v6=0;v7=0;v8=0;v9=0;v10=0;
        M_ma = i;
%         options = gaoptimset('Generations',2200,...
%             'PopInitRange',[-1;2],'PopulationSize',2000,...
%             'mutationfcn',@mutationadaptfeasible,'TolCon',1e-12);

        options = gaoptimset('Generations',10000,...
            'PopInitRange',[0;1.6],'PopulationSize',200,...
            'mutationfcn',@mutationadaptfeasible,'TolCon',1e-6,...
             'SelectionFcn',@garollet);
        [tetha(:,(k)),f(k)] = ga(@fitness_ga_multilevel...
            ,7,[],[],[],[],LB,UB,@simple_constraint_multilevel,options)
%         if j == 1
%             temp1 = tetha(:,(k));
%             temp = f(:,(k));
%         end 
%         if f(:,(k)) < temp
%             temp = f(:,(k));
%             temp1 = tetha(:,(k));
%         end
%     end
%     f(:,(k)) = temp;
%     tetha(:,(k)) = temp1;

%%
    Df = 0;
   for jj = 1 : 10
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
        Df = Df+(V22(jj)/((6*jj-1)^2))^2+(V33(jj)/((6*jj+1)^2))^2;
   end
for ii = 1:1:step1
    M = M +v11_dc*(-1)^(ii-1)*cos(tetha(ii,k));
    v2 = v2+ v11_dc*(-1)^(ii-1)*cos(5*tetha(ii,k));
    v3 = v3+v11_dc* (-1)^(ii-1)*cos(7*tetha(ii,k));
    v4 = v4+v11_dc*(-1)^(ii-1)*cos(11*tetha(ii,k));
    v5 = v5+v11_dc*(-1)^(ii-1)*cos(13*tetha(ii,k));
    v6 = v6+v11_dc*(-1)^(ii-1)*cos(17*tetha(ii,k));
    v7 = v7+v11_dc*(-1)^(ii-1)*cos(19*tetha(ii,k));
    v8 = v8+v11_dc*(-1)^(ii-1)*cos(23*tetha(ii,k));
    v9 = v9+v11_dc*(-1)^(ii-1)*cos(25*tetha(ii,k));
    v10 = v10+v11_dc*(-1)^(ii-1)*cos(29*tetha(ii,k));
end
for ii = (step1+1):1:step2
        M = M+ v22_dc*(-1)^(ii)*cos(tetha(ii,k));
        v2 = v2+  v22_dc*(-1)^(ii)*cos(5*tetha(ii,k));
        v3 = v3+  v22_dc*(-1)^(ii)*cos(7*tetha(ii,k));
        v4 = v4+  v22_dc*(-1)^(ii)*cos(11*tetha(ii,k));
        v5 = v5+ v22_dc* (-1)^(ii)*cos(13*tetha(ii,k));
        v6 = v6+ v22_dc* (-1)^(ii)*cos(17*tetha(ii,k));
        v7 = v7+  v22_dc*(-1)^(ii)*cos(19*tetha(ii,k));
        v8 = v8+  v22_dc*(-1)^(ii)*cos(23*tetha(ii,k));
        v9 = v9+ v22_dc* (-1)^(ii)*cos(25*tetha(ii,k));
        v10 = v10+  v22_dc*(-1)^(ii)*cos(29*tetha(ii,k));
end
  for ii = (step2+1):1:7
        M = M+v33_dc* (-1)^(ii-1)*cos(tetha(ii,k));
        v2 = v2+v33_dc* (-1)^(ii-1)*cos(5*tetha(ii,k));
        v3 = v3+ v33_dc*(-1)^(ii-1)*cos(7*tetha(ii,k));
        v4 = v4+ v33_dc*(-1)^(ii-1)*cos(11*tetha(ii,k));
        v5 = v5+ v33_dc*(-1)^(ii-1)*cos(13*tetha(ii,k));
        v6 = v6+v33_dc* (-1)^(ii-1)*cos(17*tetha(ii,k));
        v7 = v7+v33_dc*(-1)^(ii-1)*cos(19*tetha(ii,k));
        v8 = v8+v33_dc*(-1)^(ii-1)*cos(23*tetha(ii,k));
        v9 = v9+v33_dc*(-1)^(ii-9)*cos(25*tetha(ii,k));
        v10 = v10+v33_dc*(-1)^(ii-9)*cos(29*tetha(ii,k));
  end
   M = M/3;
    v2 = (4.*v2).*(V_dc/3).*(1/(5*pi));
    v3 = (4.*v3).*(V_dc/3).*(1/(7*pi));
    v4 = (4.*v4).*(V_dc/3).*(1/(11*pi));
    v5 = (4.*v5).*(V_dc/3).*(1/(13*pi));
    v6 = (4.*v6).*(V_dc/3).*(1/(17*pi));
    v7 = (4.*v7).*(V_dc/3).*(1/(19*pi));
    v8 = (4.*v8).*(V_dc/3).*(1/(23*pi));
    v9 = (4.*v9).*(V_dc/3).*(1/(25*pi));
    v10 = (4.*v10).*(V_dc/3).*(1/(29*pi));
    v1 = M*(V_dc*4/pi);
    v1_persent(k) = (abs(v1)./v1_ma)*100;
    v2_persent(k) = (abs(v2)./v1_ma)*100;
    v3_persent(k) = (abs(v3)./v1_ma)*100;
    v4_persent(k) = (abs(v4)./v1_ma)*100;
    v5_persent(k) = (abs(v5)./v1_ma)*100;
    v6_persent(k) = (abs(v6)./v1_ma)*100;
    v7_persent(k) = (abs(v7)./v1_ma)*100;
    THD(k) = 100*((sqrt(v2^2+v3^2+v4^2+v5^2+v6^2+v7^2+v8^2+v9^2+v10^2))/v1);
    THD1(k) = 100*((sqrt(sum(V22.^2)+sum(V33.^2)))/v1);
    DF2(k) = (100/v1)*(sqrt((v2/25)^2+(v3/49)^2+(v4/121)^2+(v5/169)^2+...
        (v6/289)^2+(v7/361)^2+(v8/529)^2+(v9/625)^2+(v10/841)^2));
    DF21(k) = (100/v1)*(sqrt(Df));
    xt(k) = i;
    k = k+1;
end

%
tetha(1,:) = mod(((tetha(1,:)./pi).*180),360);
tetha(2,:) = mod(((tetha(2,:)./pi).*180),360);
tetha(3,:) = mod(((tetha(3,:)./pi).*180),360);
tetha(4,:) = mod(((tetha(4,:)./pi).*180),360);
tetha(5,:) = mod(((tetha(5,:)./pi).*180),360);
tetha(6,:) = mod(((tetha(6,:)./pi).*180),360);
tetha(7,:) = mod(((tetha(7,:)./pi).*180),360);
% tetha(8,:) = (tetha(8,:)./pi).*180;
% tetha(9,:) = (tetha(9,:)./pi).*180;
% tetha(10,:) = (tetha(10,:)./pi).*180;
% tetha(11,:) = (tetha(11,:)./pi).*180;

%
figure(1*core)
semilogy(xt,f)
figure (2*core)
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
figure(3*core)
plot(xt,v1_persent,xt,v2_persent,xt,v3_persent,xt,v4_persent,xt,v5_persent,xt,v6_persent,xt,v7_persent);
grid on
legend('V1','Vh5','Vh7','Vh11','Vh13','Vh17','Vh19')

%save('211.txt','tetha','v1_persent','v2_persent','v3_persent',...
%    'v4_persent','v5_persent','v6_persent','v7_persent','f','THD','DF2');
%%
THD
THD1
DF2
DF21
f
tetha
v1_persent
v2_persent
v3_persent
v4_persent
v5_persent
v6_persent
v7_persent
% if core == 1
%             temp1 = tetha(:,(k));
%             temp = f(:,(k));
%         end 
%         if f(:,(k)) < temp
%             temp = f(:,(k));
%             temp1 = tetha(:,(k));
%             temp2 = THD1;
%             temp3 = DF21;
%             temp4 = DF2;
%             temp5 = THD;
%             temp6 = v1_persent;
%             temp7 = v2_persent;
%             temp8 = v3_persent;
%             temp9 = v4_persent;
%             temp10 = v5_persent;
%             temp11 = v6_persent;
%             temp12 = v7_persent;
%         end
%     end
%     f(:,(k)) = temp;
%     tetha(:,(k)) = temp1;
   
% end