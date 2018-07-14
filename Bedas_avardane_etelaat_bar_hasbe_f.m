clc
clear all;
for vol = 1: 1 :27
    d1 = '6';
d0 = 'D:\E Book\lib\Uni books\ARSHAD\the projec of the Master\variables\';
d3 = 'core';
d4 = '_';
d5 = 'vol';
d2 = num2str(vol);
d6 = '.mat';
string esm;
esm = [d0 d3 d1 d4 d5 d2 d6];
load(esm);

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
Etelaat_koli(10,i) = f_koli(w(i),i);
noe_SWICING1(i) = w(i);
% TETHA_koli = TETHA(:,i,w(i));
TETHA_koli(:,i) = TETHA(:,i,w(i));

end
beep()
Jadval(:,:,vol) = Etelaat_koli;
TETHA_KAMEL_d(:,:,vol) = TETHA_koli;
NOE_SWICING(vol) = w(i);
NOE_SWICING1(vol,:) = noe_SWICING1;
VORUDI(1,:) = mio_koli;
VORUDI(2,vol) = v11_dc ;
VORUDIz(3,vol) = v22_dc ;
VORUDI(4,vol) = v33_dc ;
end
%%

save('D:\E Book\lib\Uni books\ARSHAD\the projec of the Master\variables\VORUDI_ANFIS_f.mat')
%%
k = 1;
for i = 1:27
    for j = 1:25
        ss1(k) = TETHA_KAMEL_f(1,j,i);
        ss2(k) = TETHA_KAMEL_f(2,j,i);
        ss3(k) = TETHA_KAMEL_f(3,j,i);
        ss4(k) = TETHA_KAMEL_f(4,j,i);
        ss5(k) = TETHA_KAMEL_f(5,j,i);
        ss6(k) = TETHA_KAMEL_f(6,j,i);
        ss7(k) = TETHA_KAMEL_f(7,j,i);
        ss8(k) = NOE_SWICING1_f(i,j);
        k = k+1;
    end
end
ss1 = ss1';ss2 = ss2';ss3 = ss3';ss4 = ss4';ss5 = ss5';
ss6 = ss6';ss7 = ss7';ss8 = ss8';
ANFIS_1 = [JVORUDI ss1];
ANFIS_2 = [JVORUDI ss2];
ANFIS_3 = [JVORUDI ss3];
ANFIS_4 = [JVORUDI ss4];
ANFIS_5 = [JVORUDI ss5];
ANFIS_6 = [JVORUDI ss6];
ANFIS_7 = [JVORUDI ss7];
ANFIS_8 = [JVORUDI ss8];


%%
clc

% for i = 1:26
%     nemune1(i,:) = ANFIS_1(i,:);
% end

    g = VORUDI(1,:);

g = g'
w = nemune1(:,5);

kmos = [g w]
plot(g,w)

%%

l = 1;
for i = 1:25
    for k = 0:26
        ANFIS_1_j(l,:) = ANFIS_1(25*k+i,:);
        l = l+1;
    end
end

%%
o = 1;
for i = 1:675
    if ANFIS_1(i,1) ~= 0
            kh_2(o,:) = ANFIS_1(i,:);
            o = o+1;
    end
end
%%
o = 1;
for i = 1 : 675
    if kh(i,1) == 0 
       tes(o,:) = ANFIS_1(i,:);
       o = o+1;
    end
end

%%
plot(g',kk50(1,:),g',kk50(2,:),'r',g',kk50(3,:),'g',g',kk50(4,:),'b',...
   g',kk50(5,:),'y',g',kk50(6,:),'m',g',kk50(7,:),'c' );
legend('teta1','teta2','teta3','teta4','teta5','teta6','teta7');
xlabel ('Modulation Index');
ylabel('Optimized Swiching Angeles');
%%

g = VORUDI(1,:);
g = g';
g = g (1:25);
plot(g',Jadval_f(1,1:25,14),g',(Jadval_f(2,1:25,14)*100),g',Jadval_f(3,1:25,14),g',Jadval_f(4,1:25,14),...
   g',Jadval_f(5,1:25,14),g',Jadval_f(6,1:25,14),g',Jadval_f(7,1:25,14),g',Jadval_f(8,1:25,14),g',Jadval_f(9,1:25,14) );
legend('THD','DF2','V1','Vh5','Vh7','Vh11','Vh13','Vh17','Vh19');
xlabel ('Modulation Index');
ylabel('Harmonics Amplitude Percent');
ylim ([0 110])
xlim([0 1.2])

%%
g = VORUDI(1,:);
g = g';
g = g (1:25);
for i = 1 : 2
    if i ==2
       gama(:,i) = (Jadval_f(2,1:25,14)*100);
    else
       gama(:,i) = Jadval_f(i,1:25,14);
    end
end
gama
bar(g',gama );
legend('THD','DF2','V1','Vh5','Vh7','Vh11','Vh13','Vh17','Vh19');
xlabel ('Modulation Index');
ylabel('THD & DF2 Percent');
ylim ([0 110])
xlim([0 1.2])

%%

g = VORUDI(1,:);
g = g';
g = g (1:25);
plot(g',Jadval_f(3,1:25,14),g',Jadval_f(4,1:25,14),'--',...
   g',Jadval_f(5,1:25,14),':',g',Jadval_f(6,1:25,14),'b--o',g',Jadval_f(7,1:25,14),'c*',g',Jadval_f(8,1:25,14),'--gs',g',Jadval_f(9,1:25,14),'*' );
legend('V1','Vh5','Vh7','Vh11','Vh13','Vh17','Vh19');
xlabel ('Modulation Index');
ylabel('Harmonics Amplitude Percent');
ylim ([0 110])
xlim([0 1.2])

%%

for i = 1: 9
    kk_ss (i,:) = Jadval_f(i,1:25,14);
end
kk_ss(2,:) = kk_ss(2,:).*100;

%%
plot(g',TETHA_KAMEL_f(1,1:25,14),g',TETHA_KAMEL_f(2,1:25,14),'--',...
    g',TETHA_KAMEL_f(3,1:25,14),':',g',TETHA_KAMEL_f(4,1:25,14),'*',...
    g',TETHA_KAMEL_f(5,1:25,14),'b--o',g',TETHA_KAMEL_f(6,1:25,14),'c*',...
    g',TETHA_KAMEL_f(7,1:25,14),'--gs')
xlabel('Modulation Index');
ylabel('Optimized Switching Angles');
xlim([0 1.2])
legend('Tetha1','Tetha2','Tetha3','Tetha4','Tetha5','Tetha6','Tetha7');


%%
plot(g',v1_persent_koli(1,1:25),g',v2_persent_koli(1,1:25),...
    g',v3_persent_koli(1,1:25),g',v4_persent_koli(1,1:25),...
    g',v5_persent_koli(1,1:25),g',v6_persent_koli(1,1:25),...
    g',v7_persent_koli(1,1:25));
xlim([0 1.2])
legend('V1','Vh5','Vh7','Vh11','Vh13','Vh17','Vh19');
xlabel('Modulation Index');
ylabel('Harmonics Amplitude Percent');