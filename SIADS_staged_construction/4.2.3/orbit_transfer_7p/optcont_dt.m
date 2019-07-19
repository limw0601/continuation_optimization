function J = optcont_dt(t, x, p) 

global ncheb tT

J = zeros(6,numel(t));

dudt = zeros(3,numel(t));
for i=1:3
    for j=2:ncheb
        dudt(i,:) = dudt(i,:) + p((i-1)*ncheb+j,:).*dchebyshevT(j-1,t);
    end
end
du1dt = dudt(1,:);
du2dt = dudt(2,:);
du3dt = dudt(3,:);

J(4,:) = du1dt;
J(5,:) = du2dt;
J(6,:) = du3dt;

J  = tT*J;

end