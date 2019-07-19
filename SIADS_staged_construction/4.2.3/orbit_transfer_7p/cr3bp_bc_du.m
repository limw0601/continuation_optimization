function [data, J] = cr3bp_bc_du(prob, data, u) %#ok<INUSL>

% x0   = u(1:6);
% x1 = u(7:12);
% 
% y = x0-x1;
J = [eye(6) -eye(6)];

end