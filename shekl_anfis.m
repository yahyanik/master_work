clc;
x = (1:135)';
y = ANFIS_test(:,5);
kk = Unti;
zzzz = ANFIS_test(:,1:4);
plot(x,y,'o',x,evalfis(zzzz,Unti),'*')
legend('test data','FIS output')