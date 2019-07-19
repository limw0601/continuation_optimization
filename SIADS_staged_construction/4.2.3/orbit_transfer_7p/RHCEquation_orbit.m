function dx=RHCEquation_orbit(t,x,p)
mu    = 3.04036e-6;                                                        % 地球质量/(地球质量+太阳质量)mu = 3.0404234945077e-6;
gama = 1.00782e-2;                                                         % gama:L2点距离地球的距离(1.507677260085612e6km)

d1 = sqrt((x(1)+1+1/gama)^2+x(2)^2+x(3)^2);
d2 = sqrt((x(1)+1)^2+x(2)^2+x(3)^2);
dx = zeros(6,1);
dx(1) = x(4);
dx(2) = x(5);
dx(3) = x(6);
dx(4) = x(1)+2*x(5)-(1-mu)*(x(1)+1+1/gama)/gama^3/d1^3-mu*(x(1)+1)/gama^3/d2^3+(1-mu+gama)/gama;
dx(5) = x(2)-2*x(4)-(1-mu)/gama^3/d1^3*x(2)-mu/gama^3/d2^3*x(2);
dx(6) = -(1-mu)/gama^3/d1^3*x(3)-mu/gama^3/d2^3*x(3);