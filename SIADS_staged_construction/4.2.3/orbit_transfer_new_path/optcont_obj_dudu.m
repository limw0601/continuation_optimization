function [data, dJ] = optcont_obj_dudu(prob, data, u) %#ok<INUSD,INUSL>

global tT
dJ(1,:,:) = tT*data.W;

end 
