function [data y] = obj_dfunc(prob, data, u)

x1 = u(1);
x2 = u(2);

y = [2*(x1-2) 4*(x2-1)];


end