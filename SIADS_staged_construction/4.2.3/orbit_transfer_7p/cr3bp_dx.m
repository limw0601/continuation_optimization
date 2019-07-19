function J = cr3bp_dx(x, p) %#ok<INUSL>
%NONLINODE_DFDX   'coll'-compatible encoding of Jacobian with respect to state.


mu   = 3.04036e-6;
gama = 1.00782e-2;
x1 = x(1,:);
x2 = x(2,:);
x3 = x(3,:);
d1 = sqrt((x1+1+1/gama).^2+x2.^2+x3.^2);
d2 = sqrt((x1+1).^2+x2.^2+x3.^2);
d1dx1 = (x1+1+1/gama)./d1;
d1dx2 = x2./d1;
d1dx3 = x3./d1;
d2dx1 = (x1+1)./d2;
d2dx2 = x2./d2;
d2dx3 = x3./d2;

d1dx1dx1 = 1./d1-(x1+1+1/gama)./d1.^2.*d1dx1;
d1dx1dx2 = -(x1+1+1/gama)./d1.^2.*d1dx2;
d1dx1dx3 = -(x1+1+1/gama)./d1.^2.*d1dx3;
d1dx2dx2 = 1./d1-x2./d1.^2.*d1dx2;
d1dx2dx3 = -x2./d1.^2.*d1dx3;
d1dx3dx3 = 1./d1-x3./d1.^2.*d1dx3;
d2dx1dx1 = 1./d2-(x1+1)./d2.^2.*d2dx1;
d2dx1dx2 = -(x1+1)./d2.^2.*d2dx2;
d2dx1dx3 = -(x1+1)./d2.^2.*d2dx3;
d2dx2dx2 = 1./d2-x2./d2.^2.*d2dx2;
d2dx2dx3 = -x2./d2.^2.*d2dx3;
d2dx3dx3 = 1./d2-x3./d2.^2.*d2dx3;

dEdx1dx1 = -1-2*(1-mu)/gama^3./d1.^3.*d1dx1.*d1dx1+(1-mu)/gama^3./d1.^2.*d1dx1dx1-2*mu/gama^3./d2.^3.*d2dx1.*d2dx1+mu/gama^3./d2.^2.*d2dx1dx1;
dEdx1dx2 = -2*(1-mu)/gama^3./d1.^3.*d1dx2.*d1dx1+(1-mu)/gama^3./d1.^2.*d1dx1dx2-2*mu/gama^3./d2.^3.*d2dx2.*d2dx1+mu/gama^3./d2.^2.*d2dx1dx2;
dEdx1dx3 = -2*(1-mu)/gama^3./d1.^3.*d1dx3.*d1dx1+(1-mu)/gama^3./d1.^2.*d1dx1dx3-2*mu/gama^3./d2.^3.*d2dx3.*d2dx1+mu/gama^3./d2.^2.*d2dx1dx3;
dEdx2dx1 = dEdx1dx2;
dEdx2dx2 = -1-2*(1-mu)/gama^3./d1.^3.*d1dx2.*d1dx2+(1-mu)/gama^3./d1.^2.*d1dx2dx2-2*mu/gama^3./d2.^3.*d2dx2.*d2dx2+mu/gama^3./d2.^2.*d2dx2dx2;
dEdx2dx3 = -2*(1-mu)/gama^3./d1.^3.*d1dx3.*d1dx2+(1-mu)/gama^3./d1.^2.*d1dx2dx3-2*mu/gama^3./d2.^3.*d2dx3.*d2dx2+mu/gama^3./d2.^2.*d2dx2dx3;
dEdx3dx1 = dEdx1dx3;
dEdx3dx2 = dEdx2dx3;
dEdx3dx3 = -2*(1-mu)/gama^3./d1.^3.*d1dx3.*d1dx3+(1-mu)/gama^3./d1.^2.*d1dx3dx3-2*mu/gama^3./d2.^3.*d2dx3.*d2dx3+mu/gama^3./d2.^2.*d2dx3dx3;

J = zeros(6,6,numel(x1));
J(1,1,:) = gama*p.*dEdx1dx1;
J(1,2,:) = gama*p.*dEdx1dx2;
J(1,3,:) = gama*p.*dEdx1dx3;
J(1,4,:) = 1;
J(2,1,:) = gama*p.*dEdx2dx1;
J(2,2,:) = gama*p.*dEdx2dx2;
J(2,3,:) = gama*p.*dEdx2dx3;
J(2,5,:) = 1;
J(3,1,:) = gama*p.*dEdx3dx1;
J(3,2,:) = gama*p.*dEdx3dx2;
J(3,3,:) = gama*p.*dEdx3dx3;
J(3,6,:) = 1;
J(4,1,:) = 1-(1-mu)/gama^3./d1.^3+3*(1-mu)*(x1+1+1/gama)/gama^3./d1.^4.*d1dx1...
    -mu/gama^3./d2.^3+3*mu*(x1+1)/gama^3./d2.^4.*d2dx1;
J(4,2,:) = 3*(1-mu)*(x1+1+1/gama)./gama^3./d1.^4.*d1dx2...
    +3*mu*(x1+1)/gama^3./d2.^4.*d2dx2;
J(4,3,:) = 3*(1-mu)*(x1+1+1/gama)./gama^3./d1.^4.*d1dx3...
    +3*mu*(x1+1)/gama^3./d2.^4.*d2dx3;
J(4,4,:) = gama*p;
J(4,5,:) = 2;
J(5,1,:) = 3*(1-mu)*x2./gama^3./d1.^4.*d1dx1...
    +3*mu*x2/gama^3./d2.^4.*d2dx1;
J(5,2,:) = 1-(1-mu)./gama^3./d1.^3+3*(1-mu)*x2./gama^3./d1.^4.*d1dx2...
    -mu/gama^3./d2.^3+3*mu*x2/gama^3./d2.^4.*d2dx2;
J(5,3,:) = 3*(1-mu)*x2./gama^3./d1.^4.*d1dx3...
    +3*mu*x2/gama^3./d2.^4.*d2dx3;
J(5,4,:) = -2;
J(5,5,:) = gama*p;
J(6,1,:) = 3*(1-mu)*x3./gama^3./d1.^4.*d1dx1...
    +3*mu*x3/gama^3./d2.^4.*d2dx1;
J(6,2,:) = 3*(1-mu)*x3./gama^3./d1.^4.*d1dx2...
    +3*mu*x3/gama^3./d2.^4.*d2dx2;
J(6,3,:) = -(1-mu)./gama^3./d1.^3+3*(1-mu)*x3./gama^3./d1.^4.*d1dx3...
    -mu/gama^3./d2.^3+3*mu*x3/gama^3./d2.^4.*d2dx3;
J(6,6,:) = gama*p;

end