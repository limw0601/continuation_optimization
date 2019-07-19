function y = cr3bp(x, p)
%CATN   'coll'-compatible encoding of catenary vector field

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
% U = 0.5*((x1+1+(1-mu)/gama).^2+x2.^2)+(1-mu)/gama^3./d1+mu/gama^3./d2;
% E = 0.5*(v1.^2+v2.^2+v3.^2)-U-mu*(1-mu)/2/gama^2;

dEdx1 = -(x1+1+(1-mu)/gama)+(1-mu)/gama^3./d1.^2.*d1dx1+mu/gama^3./d2.^2.*d2dx1;
dEdx2 = x2+(1-mu)/gama^3./d1.^2.*d1dx2+mu/gama^3./d2.^2.*d2dx2;
dEdx3 = (1-mu)/gama^3./d1.^2.*d1dx3+mu/gama^3./d2.^2.*d2dx3;
dEdx4 = x4;
dEdx5 = x5;
dEdx6 = x6;

dx = zeros(6,numel(x1));
dx(1,:) = x4+gama*p.*dEdx1;
dx(2,:) = x5+gama*p.*dEdx2;
dx(3,:) = x6+gama*p.*dEdx3;
dx(4,:) = x1+2*x5-(1-mu)*(x1+1+1/gama)./gama^3./d1.^3-mu*(x1+1)/gama^3./d2.^3+(1-mu+gama)/gama+gama*p.*dEdx4;
dx(5,:) = x2-2*x4-(1-mu)/gama^3./d1.^3.*x2-mu./gama^3./d2.^3.*x2+gama*p.*dEdx5;
dx(6,:) = -(1-mu)/gama^3./d1.^3.*x3-mu./gama^3./d2.^3.*x3+gama*p.*dEdx6;

y = dx;

end
