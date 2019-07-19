function [data, y] = optcont_bc(prob, data, u) %#ok<INUSL>

x0_halo1   = u(1:6);
x0_connect = u(7:12);
x1_connect = u(13:18);
x0_halo2   = u(19:24);

y = [x0_halo1-x0_connect; x1_connect-x0_halo2];

end