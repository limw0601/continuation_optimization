function [data, y] = optcont_ini(prob, data, u) %#ok<INUSL>

x0 = u(1:6);
x1 = u(7:12);
T0 = u(13);
T  = u(14); 

y = [x0-data.xa; x1-data.xb; T0+1; T-2];% periodicity in R2 x S1 on Poincare section

end