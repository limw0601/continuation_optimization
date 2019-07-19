function [data, y] = optcont_obj(prob, data, u) %#ok<INUSL>

global tT
y = u'*data.W*u/2;
y = tT*y;

end