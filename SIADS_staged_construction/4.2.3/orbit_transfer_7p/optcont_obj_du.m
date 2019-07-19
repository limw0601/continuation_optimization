function [data, J] = optcont_obj_du(prob, data, u) %#ok<INUSD,INUSL>

global tT
J = data.W*u;
J = J';
J = tT*J;

end