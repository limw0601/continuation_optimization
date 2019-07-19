function [data, y] = cr3bp_bc(prob, data, u) %#ok<INUSL>

x0   = u(1:6);
x1 = u(7:12);

y = x0-x1;

end