function [data y] = obj_func(prob, data, u)

x1 = u(1);
x2 = u(2);

y = (x1-2)^2+2*(x2-1)^2;

end