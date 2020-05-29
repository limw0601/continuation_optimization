function y = dchebyshevT(j,t)

if j<=10
% y = (j*sin(j*acos(t)))./(1 - t.^2).^(1/2); % 0/0 at t=-1 and t=1
y = j*mychebyshevU(j-1,t);
y = y/(pi/2)^0.5;

else
y = coco_ezDFDP('f(x,p)v', @mychebyshevT, j+zeros(size(t)),t);
y = reshape(y,size(t));
end

end