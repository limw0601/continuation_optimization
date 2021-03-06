function y = ddchebyshevT(j,t)

% y = (j*sin(j*acos(t)))./(1 - t.^2).^(1/2); % 0/0 at t=-1 and t=1
if j<=10
y = j*dmychebyshevU(j-1,t);
y = y/(pi/2)^0.5;

else
y = coco_ezDFDP('f(x,p)v', @dchebyshevT, j,t);
end
end