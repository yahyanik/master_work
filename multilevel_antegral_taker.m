clc
v1_ma = 305.5775;
t = time.data;
a1 = A1.data;
w = sum(A1);
a5 = A5.data;
w2 = sum(A5);
a7 = A7.data;
a11 = A11.data;
a13 = A13.data;
a17 = A17.data;
a19 = A19.data;
a23 = A23.data;
a25 = A25.data;
a29 = A29.data;
a31 = A31.data;
a35 = A35.data;
a37 = A37.data;
a41 = A41.data;
a43 = A43.data;
a47 = A47.data;
a49 = A49.data;

b1 = B1.data;
t1 = t.*(2*50*pi);

for i = 2:1:size(t)
    c1(i) = (-t1(i-1)+t1(i))*a1(i);
    c2(i) = (-t1(i-1)+t1(i))*a5(i);
    c3(i) = (-t1(i-1)+t1(i))*a7(i);
    c4(i) = (-t1(i-1)+t1(i))*a11(i);
    c5(i) = (-t1(i-1)+t1(i))*a13(i);
    c6(i) = (-t1(i-1)+t1(i))*a17(i);
    c7(i) = (-t1(i-1)+t1(i))*a19(i);
    c8(i) = (-t1(i-1)+t1(i))*a23(i);
    c9(i) = (-t1(i-1)+t1(i))*a25(i);
    c10(i) = (-t1(i-1)+t1(i))*a29(i);
    c11(i) = (-t1(i-1)+t1(i))*a31(i);
    c12(i) = (-t1(i-1)+t1(i))*a35(i);
    c13(i) = (-t1(i-1)+t1(i))*a37(i);
    c14(i) = (-t1(i-1)+t1(i))*a41(i);
    c15(i) = (-t1(i-1)+t1(i))*a43(i);
    c16(i) = (-t1(i-1)+t1(i))*a47(i);
    c17(i) = (-t1(i-1)+t1(i))*a49(i);
    h1(i) = (-t1(i-1)+t1(i))*b1(i);
    
end
v111 = (1/pi*sum(c1))
v222 = (1/pi*sum(c2))
v333 = (1/pi*sum(c3))
v444 = (1/pi*sum(c4))
v555 = (1/pi*sum(c5))
v666 = (1/pi*sum(c6))
v777 = (1/pi*sum(c7));
v888 = (1/pi*sum(c8));
v999 = (1/pi*sum(c9));
v100 = (1/pi*sum(c10));
v110 = (1/pi*sum(c11));
v120 = (1/pi*sum(c12));
v130 = (1/pi*sum(c13));
v140 = (1/pi*sum(c14));
v150 = (1/pi*sum(c15));
v160 = (1/pi*sum(c16));
v170 = (1/pi*sum(c17));
% V111 = (1/pi*sum(h1));
THD12 = 100*((sqrt(v222^2+v333^2+v444^2+v555^2+v666^2+...
    v777^2+v888^2+v999^2+v100^2+v110^2+...
v120^2+v130^2+v140^2+v150^2+v160^2+v170^2))/v111)
%     DF233 = (100/v111)*(sqrt((v222/25)^2+(v333/49)^2+(v444/121)^2+(v555/169)^2+...
%         (v666/289)^2+(v777/361)^2))
% V1 = abs(2/0.02*sum(c7))
% volt = [v111,v222,v333,v444,v555,v666,v777]
% x = [50,250,350,550,650,850,950,]
% figure()
% bar(x,volt)
% figure()
% plot(x,volt)