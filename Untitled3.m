% 
% clc
% clear all
% f = 0.02;
% t1 = 0;
% t1 = (t1/360)*f;
% t2 = 0;
% t2 = (t2/360)*f;
% t3 = 0;
% t3 = (t3/360)*f;
% t4 = 0;
% t4 = (t4/360)*f;
% tetha11 = [1,t1];
% assignin('base', 'tetha11', [1,t1])
% tetha12 = [1,t2];
% tetha13 = [1,t3];
% tetha14 = [1,t4];
% fs = [1,f];
% sim ('multilevel');
% g = v_phase.data;
% A1 = A1.data;
% v1 = abs(0.01/2*(sum(A1)));
% A5 = A5.data;
% v5 = abs(0.01/2*(sum(A5)));
% A7 = A7.data;
% v7 = abs(0.01/2*(sum(A7)));
% A11 = A11.data;
% v11 = abs(0.01/2*(sum(A11)));
% y = A1;
%%
clc
x = 0;
for i = 1 : 1 :2
    x = x + cos(i)
end
round(6.7)
clear all
x = 0;
for i = 1 : 1 :2
    x = x + cos(i)
end
round(6.7)
figure(1)
plot(x,sin(x))
figure(12)
plot(x,cos(x))

%%
A = A10.data
x = rms(A)


%%
clc
B = sym('w_%d',[4,1])
% syms ('w',7)
M = 0
M = M+ v22_dc*(-1)^(1)*cos(B(3))
M1 = diff(M,B(3))
f = subs(M1,4)

%%
clc
clear all
z = zeros(4,2);
z(:,:,1) = [1,2;3,4;5,6;7,8];



[q,w] = min(z,[],2)

%%
r = 6.2;
round(r)