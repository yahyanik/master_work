function [c, ceq] = simple_constraint(x)
   c = [x(1)+0.0088-x(2);
   x(2)+0.0088-x(3)];
   ceq = [];