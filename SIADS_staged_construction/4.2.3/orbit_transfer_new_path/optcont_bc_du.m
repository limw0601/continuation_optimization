function [data, J] = optcont_bc_du(prob, data, u) %#ok<INUSD,INUSL>

J = zeros(12,24);
J(1:6,1:6) = eye(6);
J(1:6,7:12) = -eye(6);
J(7:12,13:18) = eye(6);
J(7:12,19:24) = -eye(6);

end