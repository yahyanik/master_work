function [c, ceq] = simple_constraint_multilevel(x)
   c = [x(1)+0.0088-x(2);
       x(2)+0.0088-x(3);
       x(3)+0.0088-x(4);
       x(4)+0.0088-x(5);
       x(5)+0.0088-x(6);
       x(6)+0.0088-x(7)];
   ceq = [];