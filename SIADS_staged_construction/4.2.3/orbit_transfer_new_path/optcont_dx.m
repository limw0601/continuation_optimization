function J = optcont_dx(t, x, p) 

global tT

J = zeros(6,6,numel(t));
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

J(1,4,:) = 1;
J(2,5,:) = 1;
J(3,6,:) = 1;
J(4,1,:) = 1-(1-mu)/gama^3./d1.^3+3*(1-mu)*(x1+1+1/gama)/gama^3./d1.^4.*d1dx1...
    -mu/gama^3./d2.^3+3*mu*(x1+1)/gama^3./d2.^4.*d2dx1;
J(4,2,:) = 3*(1-mu)*(x1+1+1/gama)./gama^3./d1.^4.*d1dx2...
    +3*mu*(x1+1)/gama^3./d2.^4.*d2dx2;
J(4,3,:) = 3*(1-mu)*(x1+1+1/gama)./gama^3./d1.^4.*d1dx3...
    +3*mu*(x1+1)/gama^3./d2.^4.*d2dx3;
J(4,5,:) = 2;
J(5,1,:) = 3*(1-mu)*x2./gama^3./d1.^4.*d1dx1...
    +3*mu*x2/gama^3./d2.^4.*d2dx1;
J(5,2,:) = 1-(1-mu)./gama^3./d1.^3+3*(1-mu)*x2./gama^3./d1.^4.*d1dx2...
    -mu/gama^3./d2.^3+3*mu*x2/gama^3./d2.^4.*d2dx2;
J(5,3,:) = 3*(1-mu)*x2./gama^3./d1.^4.*d1dx3...
    +3*mu*x2/gama^3./d2.^4.*d2dx3;
J(5,4,:) = -2;
J(6,1,:) = 3*(1-mu)*x3./gama^3./d1.^4.*d1dx1...
    +3*mu*x3/gama^3./d2.^4.*d2dx1;
J(6,2,:) = 3*(1-mu)*x3./gama^3./d1.^4.*d1dx2...
    +3*mu*x3/gama^3./d2.^4.*d2dx2;
J(6,3,:) = -(1-mu)./gama^3./d1.^3+3*(1-mu)*x3./gama^3./d1.^4.*d1dx3...
    -mu/gama^3./d2.^3+3*mu*x3/gama^3./d2.^4.*d2dx3;
                                                                                                                                     
J  = tT*J;

end