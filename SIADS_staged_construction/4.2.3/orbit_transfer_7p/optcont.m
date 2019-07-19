function y = optcont(t, x, p)

global ncheb tT
% ncheb = 5;
mu   = 3.04036e-6;
gama = 1.00782e-2;
x1 = x(1,:);
x2 = x(2,:);
x3 = x(3,:);
x4 = x(4,:);
x5 = x(5,:);
x6 = x(6,:);


u = zeros(3,numel(t));
for i=1:3
    for j=1:ncheb
        u(i,:) = u(i,:) + p((i-1)*ncheb+j,:).*mychebyshevT(j-1,t);
    end
end

u1 = u(1,:);
u2 = u(2,:);
u3 = u(3,:);
d1 = sqrt((x1+1+1/gama).^2+x2.^2+x3.^2);
d2 = sqrt((x1+1).^2+x2.^2+x3.^2);
dx = zeros(6,numel(t));
dx(1,:) = x4;
dx(2,:) = x5;
dx(3,:) = x6;
dx(4,:) = x1+2*x5-(1-mu)*(x1+1+1/gama)./gama^3./d1.^3-mu*(x1+1)/gama^3./d2.^3+(1-mu+gama)/gama+u1;
dx(5,:) = x2-2*x4-(1-mu)/gama^3./d1.^3.*x2-mu./gama^3./d2.^3.*x2+u2;
dx(6,:) = -(1-mu)/gama^3./d1.^3.*x3-mu./gama^3./d2.^3.*x3+u3;

y = dx;
y = tT*y;

end