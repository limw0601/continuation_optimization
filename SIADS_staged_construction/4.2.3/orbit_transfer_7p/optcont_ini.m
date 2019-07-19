function [data, y] = optcont_ini(prob, data, u) %#ok<INUSL>

x0 = u(1:6);
T0 = u(7);
T  = u(8); 

y = [x0-data.xa; T0+1; T-2];% periodicity in R2 x S1 on Poincare section

end