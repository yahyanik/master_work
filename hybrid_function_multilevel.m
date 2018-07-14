clc;

% clear all;
format long
warning off
% core = 1;
% fig = 46;

% Etelaat = zeros (10,30,6);
% TETHA = zeros (7,30,6);
global step1
global step2
global V_dc
global v11_dc
global v22_dc
global v33_dc
global M_ma
global v1_ma
global v_dc_optim
v_dc_optim = 100;
% V_dc = 300;
for vol = 1: 1 :14
switch vol
    case 1
    v11_dc = 0.95;v22_dc = 0.95;v33_dc = 0.95;
    case 2
    v11_dc = 0.95;v22_dc = 0.95;v33_dc = 1;
    case 3
    v11_dc = 0.95;v22_dc = 0.95;v33_dc = 1.05;
    case 4
    v11_dc = 0.95;v22_dc = 1;v33_dc = 0.95;
    case 5
    v11_dc = 0.95;v22_dc = 1;v33_dc = 1;
    case 6
    v11_dc = 0.95;v22_dc = 1;v33_dc = 1.05;
    case 7
    v11_dc = 0.95;v22_dc = 1.05;v33_dc = 0.95;
    case 8
    v11_dc = 0.95;v22_dc = 1.05;v33_dc = 1;
    case 9
    v11_dc = 0.95;v22_dc = 1.05;v33_dc = 1.05;
     case 10
    v11_dc = 1;v22_dc = 0.95;v33_dc = 0.95;
    case 11
    v11_dc = 1;v22_dc = 0.95;v33_dc = 1;
    case 12
    v11_dc = 1;v22_dc = 0.95;v33_dc = 1.05;
    case 13
    v11_dc = 1;v22_dc = 1;v33_dc = 0.95;
    case 14
    v11_dc = 1;v22_dc = 1;v33_dc = 1;
    case 15
    v11_dc = 1;v22_dc = 1;v33_dc = 1.05;
    case 16
    v11_dc = 1;v22_dc = 1.05;v33_dc = 0.95;
    case 17
    v11_dc = 1;v22_dc = 1.05;v33_dc = 1;
    case 18
    v11_dc = 1;v22_dc = 1.05;v33_dc = 1.05;
     case 19
    v11_dc = 1.05;v22_dc = 0.95;v33_dc = 0.95;
    case 20
    v11_dc = 1.05;v22_dc = 0.95;v33_dc = 1;
    case 21
    v11_dc = 1.05;v22_dc = 0.95;v33_dc = 1.05;
    case 22
    v11_dc = 1.05;v22_dc = 1;v33_dc = 0.95;
    case 23
    v11_dc = 1.05;v22_dc = 1;v33_dc = 1;
    case 24
    v11_dc = 1.05;v22_dc = 1;v33_dc = 1.05;
    case 25
    v11_dc = 1.05;v22_dc = 1.05;v33_dc = 0.95;
    case 26
    v11_dc = 1.05;v22_dc = 1.05;v33_dc = 1;
    case 27
    v11_dc = 1.05;v22_dc = 1.05;v33_dc = 1.05;
end
    
% step1 = 1;
% step2 = 2;
for core = 1:1:6
%     clear all;
    if core == 1
        step1 = 5;
        step2 = 6;
    elseif core == 2
        step1 = 1;
        step2 = 6;
    elseif core == 3
        step1 = 1;
        step2 = 2;
    elseif core == 4
        step1 = 3;
        step2 = 6;
    elseif core == 5
        step1 = 3;
        step2 = 4;
    elseif core == 6
        step1 = 1;
        step2 = 4;
    end
    %%

% parpool('local')
step = 0.04;   %mezrabe 10 bashad
% m_return = input()

    

% s = 3;
% v_dc_optim = 100;
% % k = (1-0.05+rand(1,s)./10);
% v11_dc = 1;
% v22_dc = 1;
% v33_dc = 1;
V_dc = v_dc_optim*sum(v11_dc+v22_dc+v33_dc);

k = 1;
LB = [0,0,0,0,0,0,0];
UB = [pi/2,pi/2,pi/2,pi/2,pi/2,pi/2,pi/2];
% LB = [0,0,0,0,0,0,0,0,0,0,0,1,1,1];
% UB = [pi/2,pi/2,pi/2,pi/2,pi/2,pi/2,pi/2,pi/2,pi/2,pi/2,pi/2,7,7,7];
%%
% tetha = zeros(7,30);
% f = zeros(1,30);
% v1_persent = zeros(1,30);
% v2_persent = zeros(1,30);
% v3_persent = zeros(1,30);
% v4_persent = zeros(1,30);
% v5_persent =zeros(1,30);
% v6_persent = zeros(1,30);
% v7_persent = zeros(1,30);
% xt = zeros(1,30);
% THD1 = zeros(1,30);
% DF2 = zeros(1,30);
% DF21 = zeros(1,30);
%%
for i = 0.04:step:1
     temp = 0;
     temp1 = 0;
     M_ma = i;
%      Jacubian()
     
     %%
      for j =1:1:10
        
        
    
%         options = gaoptimset('Generations',2200,...
%             'PopInitRange',[-1;2],'PopulationSize',2000,...
%             'mutationfcn',@mutationadaptfeasible,'TolCon',1e-12);

        options = gaoptimset('Generations',1000,...
                    'PopInitRange',[0;1.6],...
                    'PopulationSize',70,...
                    'mutationfcn',@mutationadaptfeasible,...
                    'Vectorized','on',...
                    'TolCon',1e-4,...
                    'CrossoverFcn',@crossoverfunc,'SelectionFcn',@garollet,...
                    'TolFun',1e-4);
%             'SelectionFcn',@garollet);
if M_ma<=1
        [tetha(:,(k)),f((k))] = ga(@fitness_ga_multilevel_Hybrid...
            ,7,[],[],[],[],LB,UB,@simple_constraint_multilevel_Hybrid,options);
else
   [tetha(:,(k)),f((k))] = ga(@fitness_ga_multilevel_Hybrid1...
         ,7,[],[],[],[],LB,UB,@simple_constraint_multilevel_Hybrid,options);
end
        
        %%
        x0 = tetha(:,(k));
%             x0 = [0;0;0;0;0;0;0]
%         [te,fw] = newtonraphson(@fitness_hybrid_newtonraphson,@cal_Jacubian,x0,1e-5,1e4);
%         [a,b] = size (te);
%         tetha(:,(k)) = te(:,b);
% %         f1 = struct2cell(fw);
% %         f2 = f1{3,1};
% %         len = f1{5,1};
%         f((k)) = fitness_ga_multilevel(tetha(:,(k)));
% %         x0 = (pi/4)*(rand(4,1));
% %         [te,f] = newtonraphson(@fitness_newtonraphson,[],x0,[],[]);
        A = [1 -1 0 0 0 0 0;0 1 -1 0 0 0 0;0 0 1 -1 0 0 0
            0 0 0 1 -1 0 0;0 0 0 0 1 -1 0;0 0 0 0 0 1 -1];
        B = [0.0088;0.0088;0.0088;0.0088;0.0088;0.0088];
        
        options = optimoptions('fmincon','Algorithm','sqp',...
            'TolFun',1e-20,'TolX',1e-20,'TolCon',1e-20,'MaxFunEvals',...
            1e4,'OutputFcn',@outfun);
        setappdata(0,'FMINCONsStopFlag',false);
        T = timer('startdelay',10,'timerfcn',@(src,evt)setappdata(0,'FMINCONsStopFlag',true));
        start(T);
        if M_ma<=1
        x = fmincon(@fitness_ga_multilevel,x0,[],[],[],[],LB,UB,@simple_constraint_multilevel,options);
        else
            x = fmincon(@fitness_ga_multilevel1,x0,[],[],[],[],LB,UB,@simple_constraint_multilevel,options);
        end
        [~,b] = size (x);
        tetha(:,(k)) = x(:,b);
        f((k)) = fitness_ga_multilevel(tetha(:,(k)));
%%
            if j == 1
            temp1 = tetha(:,(k));
            temp = f(k);
           
            end 
    if f(k) < temp && ~isnan(f(k))
            temp = f(k);
            temp1 = tetha(:,(k));
         
    end
     end
    
    f(k) = temp;
   
    tetha(:,(k)) = temp1;

%%
    Df = 0;
%     v22 = zeros(1,10);
%     v33 = zeros(1,10);
%     V22 = zeros(1,10);
%     V33 = zeros(1,10);
   for jj = 1 : 20
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
        V22(jj) = (4.*v22(jj)).*(v_dc_optim).*(1/((6*jj-1)*pi));
        V33(jj) = (4.*v33(jj)).*(v_dc_optim).*(1/((6*jj+1)*pi));
        Df = Df+(V22(jj)/((6*jj-1)^2))^2+(V33(jj)/((6*jj+1)^2))^2;
   end

%    for ii = 1:1:step1
%     M = M +v11_dc*(-1)^(ii-1)*cos(tetha(ii,k));
%    end
%     for ii = (step1+1):1:step2
%         M = M+ v22_dc*(-1)^(ii)*cos(tetha(ii,k));
%     end
%   for ii = (step2+1):1:7
%         M = M+v33_dc* (-1)^(ii-1)*cos(tetha(ii,k));
%   end
%     v1 = (M/3)*(V_dc*4/pi);
    
    
%     THD(k) = 100*((sqrt(v2^2+v3^2+v4^2+v5^2+v6^2+v7^2+v8^2+v9^2+v10^2))/v1);
   
%     V = V22+V33;
%     x = thd(V)
%     DF2(k) = (100/v1)*(sqrt((v2/25)^2+(v3/49)^2+(v4/121)^2+(v5/169)^2+...
%         (v6/289)^2+(v7/361)^2+(v8/529)^2+(v9/625)^2+(v10/841)^2));
   
%     
%     %%
    M=0;v2=0;v3=0;v4=0;v5=0;v6=0;v7=0;v8=0;v9=0;v10=0;
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
    v2 = (4.*v2).*(v_dc_optim).*(1/(5*pi));
    v3 = (4.*v3).*(v_dc_optim).*(1/(7*pi));
    v4 = (4.*v4).*(v_dc_optim).*(1/(11*pi));
    v5 = (4.*v5).*(v_dc_optim).*(1/(13*pi));
    v6 = (4.*v6).*(v_dc_optim).*(1/(17*pi));
    v7 = (4.*v7).*(v_dc_optim).*(1/(19*pi));
    v8 = (4.*v8).*(v_dc_optim).*(1/(23*pi));
    v9 = (4.*v9).*(v_dc_optim).*(1/(25*pi));
    v10 = (4.*v10).*(v_dc_optim).*(1/(29*pi));
    v1 = M*(V_dc*4/pi);
     THD1(k) = 100*((sqrt(sum(V22.^2)+sum(V33.^2)))/v1);
      DF21(k) = (100/v1)*(sqrt(Df));

    v1_persent(k) = (abs(v1)./v1_ma)*100;
    v2_persent(k) = (abs(v2)./v1_ma)*100;
    v3_persent(k) = (abs(v3)./v1_ma)*100;
    v4_persent(k) = (abs(v4)./v1_ma)*100;
    v5_persent(k) = (abs(v5)./v1_ma)*100;
    v6_persent(k) = (abs(v6)./v1_ma)*100;
    v7_persent(k) = (abs(v7)./v1_ma)*100;
% %     
    mio(k) = i;
    k = k+1;
    
    
end
%  mio(k) = 1.1;
% 
% %%
tetha(1,:) = mod(((tetha(1,:)./pi).*180),360);
tetha(2,:) = mod(((tetha(2,:)./pi).*180),360);
tetha(3,:) = mod(((tetha(3,:)./pi).*180),360);
tetha(4,:) = mod(((tetha(4,:)./pi).*180),360);
tetha(5,:) = mod(((tetha(5,:)./pi).*180),360);
tetha(6,:) = mod(((tetha(6,:)./pi).*180),360);
tetha(7,:) = mod(((tetha(7,:)./pi).*180),360);
% % tetha(8,:) = (tetha(8,:)./pi).*180;
% % tetha(9,:) = (tetha(9,:)./pi).*180;
% % tetha(10,:) = (tetha(10,:)./pi).*180;
% % tetha(11,:) = (tetha(11,:)./pi).*180;
% 
% %%
%  THD1
%  f
%  tetha

% d1 = num2str(core);
% d0 = 'D:\E Book\lib\Uni books\ARSHAD\the projec of the Master\variables\';
% d3 = 'core';
% d4 = '_';
% d5 = 'vol';
% d2 = num2str(vol);
% d6 = '.mat';
% string esm;
% esm = [d0 d3 d1 d4 d5 d2 d6];
% d7 = '.png';
% d8 = num2str(fig);
% esm1 = [d0 d3 d1 d4 d5 d2 d4 d8 d7];
%%
figure(fig)
semilogy(mio,f)
xlabel ('Modulation Index');
ylabel('Objective function');
% saveas(gcf,esm1)
fig = fig+1;
figure (fig)
plot(mio,tetha(1,:),mio,tetha(2,:),'r',mio,tetha(3,:),'g',mio,tetha(4,:),'b',...
   mio,tetha(5,:),'y',mio,tetha(6,:),'m',mio,tetha(7,:),'c' );
legend('teta1','teta2','teta3','teta4','teta5','teta6','teta7');
xlabel ('Modulation Index');
ylabel('Optimized Swiching Angeles');
% d8 = num2str(fig);
%  esm1 = [d0 d3 d1 d4 d5 d2 d4 d8 d7];
% saveas(gcf,esm1)
% subplot(3,3,1)
% plot(mio,tetha(1,:));
% legend('teta1')
% subplot(3,3,2)
% plot(mio,tetha(2,:),'r');
% legend('teta2')
% subplot(3,3,3)
% plot(mio,tetha(3,:),'g');
% legend('teta3')
% subplot(3,3,4)
% plot(mio,tetha(4,:),'b');
% legend('teta4')
% subplot(3,3,5)
% plot(mio,tetha(5,:),'y');
% legend('teta5')
% subplot(3,3,6)
% plot(mio,tetha(6,:),'m');
% legend('teta6')
% subplot(3,3,7)
% plot(mio,tetha(7,:),'c');
% legend('teta7')
fig = fig+1;
figure(fig)
plot(mio,v1_persent,mio,v2_persent,mio,v3_persent,mio,v4_persent,mio,v5_persent,mio,v6_persent,mio,v7_persent);
grid on
legend('V1','Vh5','Vh7','Vh11','Vh13','Vh17','Vh19')
xlabel ('Modulation Index');
ylabel('Harmonic Amplitude Percent');
% d8 = num2str(fig);
% esm1 = [d0 d3 d1 d4 d5 d2 d4 d8 d7];
% saveas(gcf,esm1)
fig = fig+1;

%save('211.txt','tetha','v1_persent','v2_persent','v3_persent',...
%    'v4_persent','v5_persent','v6_persent','v7_persent','f','THD','DF2');
%%

THD_koli(core,:) = (THD1');
DF21_koli(core,:) = (DF21');
v1_persent_koli(core,:) = (v1_persent');
v2_persent_koli(core,:) = (v2_persent');
v3_persent_koli(core,:) = (v3_persent');
v4_persent_koli(core,:) = (v4_persent');
v5_persent_koli(core,:) = (v5_persent');
v6_persent_koli(core,:) = (v6_persent');
v7_persent_koli(core,:) = (v7_persent');
f_koli(core,:) = (f');
mio_koli = (mio');
% Etelaat(2,:,core) = DF21;
% Etelaat(3,:,core) = v1_persent;
% Etelaat(4,:,core) = v2_persent;
% Etelaat(5,:,core) = v3_persent;
% Etelaat(6,:,core) = v4_persent;
% Etelaat(7,:,core) = v5_persent;
% Etelaat(8,:,core) = v6_persent;
% Etelaat(9,:,core) = v7_persent;
% Etelaat(10,:,core) = f;

TETHA(:,:,core) = tetha;


%    save(esm)
 end

   %%
[q,w] = min(THD_koli,[],1);
[a,b] = size(THD_koli);

for i = 1:b
    
Etelaat_koli(1,i) = THD_koli(w(i),i);
Etelaat_koli(2,i) = DF21_koli(w(i),i);
Etelaat_koli(3,i) = v1_persent_koli(w(i),i);
Etelaat_koli(4,i) = v2_persent_koli(w(i),i);
Etelaat_koli(5,i) = v3_persent_koli(w(i),i);
Etelaat_koli(6,i) = v4_persent_koli(w(i),i);
Etelaat_koli(7,i) = v5_persent_koli(w(i),i);
Etelaat_koli(8,i) = v6_persent_koli(w(i),i);
Etelaat_koli(9,i) = v7_persent_koli(w(i),i);
c = f_koli(w(i),i);
noe_SWICING1(i) = w(i);
% TETHA_koli = TETHA(:,i,w(i));
TETHA_koli(:,:) = TETHA(:,i,w(i));

end
beep()
Jadval(:,:,vol) = Etelaat_koli;
TETHA_KAMEL(:,:,vol) = TETHA_koli;
NOE_SWICING(vol) = w(i);
NOE_SWICING1(vol,:) = noe_SWICING1;
VORUDI(1,:) = mio_koli;
VORUDI(2,vol) = v11_dc ;
VORUDI(3,vol) = v22_dc ;
VORUDI(4,vol) = v33_dc ;
end
 

%%

% save('D:\E Book\lib\Uni books\\ARSHAD\the projec of the Master\variables\KOLI1.mat')
% dom(:,:) = Etelaat(1,:,:)
% [q,w] = min((dom'));
% Etelaat_koli = zeros(10,30);
% TETHA_koli = zeros(7,30);
% for i = 1:30
%     
% Etelaat_koli(:,:) = Etelaat(:,i,w(i));
% TETHA_koli(:,:) = TETHA(:,i,w(i));
% 
% end
% Halat_switch = w;
%%
% VORUDI = VORUDI';
% JVORUDI = 0;
be = 1;
for beta = 1: 27
for i = 1:25
    JVORUDI(be,5) = VORUDI(beta,4);
    be= be+1;
end
end
AANFIS_1 = JVORUDI;

