function [data, y] = T_bc(prob, data, u) %#ok<INUSL>


T0 = u(1);
T  = u(2);

y = [T0+1; T-2];

end