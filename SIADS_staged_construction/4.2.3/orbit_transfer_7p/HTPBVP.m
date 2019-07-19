function HTPBVP()

prob = coco_prob;
prob = coco_set(prob, 'cont', 'NAdapt', 2, 'h_max', 5);
prob1 = coco_set(prob, 'coll', 'NTST', 50);


%% Initilization by continuation methods to yield terminal states
load ini_term_states
              
la = zeros(6,1);
tf = 2;
xa  = bdata.xa;
xb  = bdata.xb;
bdata.tf = tf;
options  = odeset('reltol',1e-13,'abstol',1e-15);
[t0, x0]  = ode45(@(t,x) hamilton(t, x), [0,tf], [xa;la], options);
coll_args = {@optcont, t0, x0};
prob1 = ode_isol2coll(prob1, '', coll_args{:});
[data, uidx] = coco_get_func_data(prob1, 'coll', 'data', 'uidx');
maps = data.coll_seg.maps;
prob1 = coco_add_func(prob1, 'bc', @optcont_ini, bdata, 'zero', 'uidx', ...
  uidx([maps.x0_idx(1:6); maps.T_idx]));
prob1 = coco_add_pars(prob1, 'pars_xb', maps.x1_idx(1:6), ...
  {'xb1' 'xb2' 'xb3' 'xb4' 'xb5' 'xb6'});
prob1 = coco_add_pars(prob1, 'pars_la', maps.x0_idx(7:12), ...
  {'la1' 'la2' 'la3' 'la4' 'la5' 'la6'});
% let y1e approach to xb(1);
if x0(end,1)>bdata.xb(1)
    dx_int = [bdata.xb(1), x0(end,1)];
else
    dx_int = [x0(end,1),bdata.xb(1)];
end
bd = coco(prob1, 'collH1', [], 1, {'xb1' 'xb2' 'xb3' 'xb4' 'xb5' 'xb6' 'la1'},dx_int);


pid = [2 3 5 6 4]; % [2 3 4 5 6] fails
for k=2:6
    kk = pid(k-1);
    prunid = sprintf('collH%d',k-1);
    crunid = sprintf('collH%d',k);
    lab = coco_bd_labs(bd, 'EP');
    [lb, id] = sort(lab);
    prob1 = ode_coll2coll(prob, '', prunid, lb(end));
    [data, uidx] = coco_get_func_data(prob1, 'coll', 'data', 'uidx');
    maps = data.coll_seg.maps;
    prob1 = coco_add_func(prob1, 'bc', @optcont_ini, bdata, 'zero', 'uidx', ...
          uidx([maps.x0_idx(1:6); maps.T_idx]));
    prob1 = coco_add_pars(prob1, 'pars_xb', maps.x1_idx(1:6), ...
          {'xb1' 'xb2' 'xb3' 'xb4' 'xb5' 'xb6'});
    prob1 = coco_add_pars(prob1, 'pars_la', maps.x0_idx(7:12), ...
          {'la1' 'la2' 'la3' 'la4' 'la5' 'la6'});
    ye = sprintf('xb%d',kk);
    yy = coco_bd_col(bd,ye);
    idx = coco_bd_idxs(bd, 'EP');
    dye = yy(idx(id(end)));
    if dye>bdata.xb(kk)
      dy_int = [bdata.xb(kk) dye];
    else
      dy_int = [dye bdata.xb(kk)];
    end
    cont_args = cell(1,7);
    for kid=k:6
        cont_args{kid-k+1} = sprintf('xb%d',pid(kid-1));
    end
    for kid=1:k
        cont_args{6-k+1+kid} = sprintf('la%d',kid);
    end
    % let yke approach to xb(k)
    bd = coco(prob1, crunid, [], 1, cont_args, dy_int);
end

lab = coco_bd_labs(bd, 'EP');
sol = coll_read_solution('', crunid, max(lab));

figure(1)
plot(sol.tbp(1:5:end),-sol.xbp(1:5:end,10),'r.','LineWidth',2,'MarkerSize',8);hold on
plot(sol.tbp(1:5:end),-sol.xbp(1:5:end,11),'b.','LineWidth',2,'MarkerSize',8);
plot(sol.tbp(1:5:end),-sol.xbp(1:5:end,12),'k.','LineWidth',2,'MarkerSize',8);

figure(2)
plot(sol.tbp(1:5:end),sol.xbp(1:5:end,1),'r.','LineWidth',2,'MarkerSize',8);hold on
plot(sol.tbp(1:5:end),sol.xbp(1:5:end,2),'b.','LineWidth',2,'MarkerSize',8);
plot(sol.tbp(1:5:end),sol.xbp(1:5:end,3),'k.','LineWidth',2,'MarkerSize',8);

end

function y = optcont(x, p)

mu   = 3.04036e-6;
gama = 1.00782e-2;
x1 = x(1,:);
x2 = x(2,:);
x3 = x(3,:);
x4 = x(4,:);
x5 = x(5,:);
x6 = x(6,:);
l1 = x(7,:);
l2 = x(8,:);
l3 = x(9,:);
l4 = x(10,:);
l5 = x(11,:);
l6 = x(12,:);
u1 = -l4;
u2 = -l5;
u3 = -l6;
% original system
dx = zeros(6,numel(x1));
d1 = sqrt((x1+1+1/gama).^2+x2.^2+x3.^2);
d2 = sqrt((x1+1).^2+x2.^2+x3.^2);
dx(1,:) = x4;
dx(2,:) = x5;
dx(3,:) = x6;
dx(4,:) = x1+2*x5-(1-mu)*(x1+1+1/gama)./gama^3./d1.^3-mu*(x1+1)/gama^3./d2.^3+(1-mu+gama)/gama+u1;
dx(5,:) = x2-2*x4-(1-mu)/gama^3./d1.^3.*x2-mu./gama^3./d2.^3.*x2+u2;
dx(6,:) = -(1-mu)/gama^3./d1.^3.*x3-mu./gama^3./d2.^3.*x3+u3;
% adjoint system
dl = zeros(6,numel(x1));
d1dx1 = (x1+1+1/gama)./d1;
d1dx2 = x2./d1;
d1dx3 = x3./d1;
d2dx1 = (x1+1)./d2;
d2dx2 = x2./d2;
d2dx3 = x3./d2;
J14 = 1;
J25 = 1;
J36 = 1;
J41 = 1-(1-mu)/gama^3./d1.^3+3*(1-mu)*(x1+1+1/gama)/gama^3./d1.^4.*d1dx1...
    -mu/gama^3./d2.^3+3*mu*(x1+1)/gama^3./d2.^4.*d2dx1;
J42 = 3*(1-mu)*(x1+1+1/gama)./gama^3./d1.^4.*d1dx2...
    +3*mu*(x1+1)/gama^3./d2.^4.*d2dx2;
J43 = 3*(1-mu)*(x1+1+1/gama)./gama^3./d1.^4.*d1dx3...
    +3*mu*(x1+1)/gama^3./d2.^4.*d2dx3;
J45 = 2;
J51 = 3*(1-mu)*x2./gama^3./d1.^4.*d1dx1...
    +3*mu*x2/gama^3./d2.^4.*d2dx1;
J52 = 1-(1-mu)./gama^3./d1.^3+3*(1-mu)*x2./gama^3./d1.^4.*d1dx2...
    -mu/gama^3./d2.^3+3*mu*x2/gama^3./d2.^4.*d2dx2;
J53 = 3*(1-mu)*x2./gama^3./d1.^4.*d1dx3...
    +3*mu*x2/gama^3./d2.^4.*d2dx3;
J54 = -2;
J61 = 3*(1-mu)*x3./gama^3./d1.^4.*d1dx1...
    +3*mu*x3/gama^3./d2.^4.*d2dx1;
J62 = 3*(1-mu)*x3./gama^3./d1.^4.*d1dx2...
    +3*mu*x3/gama^3./d2.^4.*d2dx2;
J63 = -(1-mu)./gama^3./d1.^3+3*(1-mu)*x3./gama^3./d1.^4.*d1dx3...
    -mu/gama^3./d2.^3+3*mu*x3/gama^3./d2.^4.*d2dx3;

dl(1,:) = J41.*l4+J51.*l5+J61.*l6;
dl(2,:) = J42.*l4+J52.*l5+J62.*l6;
dl(3,:) = J43.*l4+J53.*l5+J63.*l6;
dl(4,:) = J14.*l1+J54.*l5;
dl(5,:) = J25.*l2+J45.*l4;
dl(6,:) = J36.*l3;

y = [dx; -dl];

end

function [data, y] = optcont_ini(prob, data, u) %#ok<INUSL>

x0 = u(1:6);
T  = u(7); 

y = [x0-data.xa; T-data.tf];% periodicity in R2 x S1 on Poincare section

end


function y = hamilton(t, x)

mu   = 3.04036e-6;
gama = 1.00782e-2;
x1 = x(1,:);
x2 = x(2,:);
x3 = x(3,:);
x4 = x(4,:);
x5 = x(5,:);
x6 = x(6,:);
l1 = x(7,:);
l2 = x(8,:);
l3 = x(9,:);
l4 = x(10,:);
l5 = x(11,:);
l6 = x(12,:);
% original system
dx = zeros(6,numel(t));
d1 = sqrt((x1+1+1/gama).^2+x2.^2+x3.^2);
d2 = sqrt((x1+1).^2+x2.^2+x3.^2);
dx(1,:) = x4;
dx(2,:) = x5;
dx(3,:) = x6;
dx(4,:) = x1+2*x5-(1-mu)*(x1+1+1/gama)./gama^3./d1.^3-mu*(x1+1)/gama^3./d2.^3+(1-mu+gama)/gama;
dx(5,:) = x2-2*x4-(1-mu)/gama^3./d1.^3.*x2-mu./gama^3./d2.^3.*x2;
dx(6,:) = -(1-mu)/gama^3./d1.^3.*x3-mu./gama^3./d2.^3.*x3;
% adjoint system
dl = zeros(6,numel(t));
d1dx1 = (x1+1+1/gama)./d1;
d1dx2 = x2./d1;
d1dx3 = x3./d1;
d2dx1 = (x1+1)./d2;
d2dx2 = x2./d2;
d2dx3 = x3./d2;
J14 = 1;
J25 = 1;
J36 = 1;
J41 = 1-(1-mu)/gama^3./d1.^3+3*(1-mu)*(x1+1+1/gama)/gama^3./d1.^4.*d1dx1...
    -mu/gama^3./d2.^3+3*mu*(x1+1)/gama^3./d2.^4.*d2dx1;
J42 = 3*(1-mu)*(x1+1+1/gama)./gama^3./d1.^4.*d1dx2...
    +3*mu*(x1+1)/gama^3./d2.^4.*d2dx2;
J43 = 3*(1-mu)*(x1+1+1/gama)./gama^3./d1.^4.*d1dx3...
    +3*mu*(x1+1)/gama^3./d2.^4.*d2dx3;
J45 = 2;
J51 = 3*(1-mu)*x2./gama^3./d1.^4.*d1dx1...
    +3*mu*x2/gama^3./d2.^4.*d2dx1;
J52 = 1-(1-mu)./gama^3./d1.^3+3*(1-mu)*x2./gama^3./d1.^4.*d1dx2...
    -mu/gama^3./d2.^3+3*mu*x2/gama^3./d2.^4.*d2dx2;
J53 = 3*(1-mu)*x2./gama^3./d1.^4.*d1dx3...
    +3*mu*x2/gama^3./d2.^4.*d2dx3;
J54 = -2;
J61 = 3*(1-mu)*x3./gama^3./d1.^4.*d1dx1...
    +3*mu*x3/gama^3./d2.^4.*d2dx1;
J62 = 3*(1-mu)*x3./gama^3./d1.^4.*d1dx2...
    +3*mu*x3/gama^3./d2.^4.*d2dx2;
J63 = -(1-mu)./gama^3./d1.^3+3*(1-mu)*x3./gama^3./d1.^4.*d1dx3...
    -mu/gama^3./d2.^3+3*mu*x3/gama^3./d2.^4.*d2dx3;

dl(1,:) = J41.*l4+J51.*l5+J61.*l6;
dl(2,:) = J42.*l4+J52.*l5+J62.*l6;
dl(3,:) = J43.*l4+J53.*l5+J63.*l6;
dl(4,:) = J14.*l1+J54.*l5;
dl(5,:) = J25.*l2+J45.*l4;
dl(6,:) = J36.*l3;

y = [dx; -dl];

end
