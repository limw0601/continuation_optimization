function J = optcont_dp(t, x, p) 

global ncheb tT

J = zeros(6,3*ncheb,numel(t));
for i=1:3
    for j=1:ncheb
        J(3+i,(i-1)*ncheb+j,:) = mychebyshevT(j-1,t);
    end
end
J  = tT*J;

end