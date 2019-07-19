function J = cr3bp_dp(x, p)
%BISTABLE_DP   'coll'-compatible encoding of Jacobian with respect to parameters.

mu   = 3.04036e-6;
gama = 1.00782e-2;
x1 = x(1,:);
x2 = x(2,:);
x3 = x(3,:);
x4 = x(4,:);
x5 = x(5,:);
x6 = x(6,:);

d1 = sqrt((x1+1+1/gama).^2+x2.^2+x3.^2);
d2 = sqrt((x1+1).^2+x2.^2+x3.^2);
d1dx1 = (x1+1+1/gama)./d1;
d1dx2 = x2./d1;
d1dx3 = x3./d1;
d2dx1 = (x1+1)./d2;
d2dx2 = x2./d2;
d2dx3 = x3./d2;

dEdx1 = -(x1+1+(1-mu)/gama)+(1-mu)/gama^3./d1.^2.*d1dx1+mu/gama^3./d2.^2.*d2dx1;
dEdx2 = -x2+(1-mu)/gama^3./d1.^2.*d1dx2+mu/gama^3./d2.^2.*d2dx2;
dEdx3 = (1-mu)/gama^3./d1.^2.*d1dx3+mu/gama^3./d2.^2.*d2dx3;
dEdx4 = x4;
dEdx5 = x5;
dEdx6 = x6;

J = zeros(6,1,numel(x1));
J(1,1,:) = dEdx1;
J(2,1,:) = dEdx2;
J(3,1,:) = dEdx3;
J(4,1,:) = dEdx4;
J(5,1,:) = dEdx5;
J(6,1,:) = dEdx6;

end