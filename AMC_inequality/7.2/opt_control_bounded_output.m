function y=opt_control_bounded_output()

% in this version, we minimize objective funtional combined with both
% output and input. We expect an unique optimal solutions. Here, we impose
% inequality constraints in terms of output. We first try the case that
% initial parameters are selected from infeasibel region. We release two
% probelem parameters to first find a fold point. We then drive
% d.p3,...,d.p_ncheb to zeros.


global ncheb tT
global M m g l TT
nterm = 10;
ncheb = nterm; 
tf = 2;
tT =(tf-0)/2;
M  = 2;
m  = 0.1;
l  = 0.5;
g  = 9.8;

diary bounded_output
diary on

prob = coco_prob;
prob = coco_set(prob, 'ode', 'autonomous', false);
prob = coco_set(prob, 'cont', 'h_max', 500);
prob = coco_set(prob, 'cont', 'MaxRes', 1.5);
prob = coco_set(prob, 'cont', 'al_max', 15);
prob = coco_set(prob, 'coll', 'NTST', 40);
prob = coco_set(prob, 'cont', 'ItMX', 250);
prob = coco_set(prob, 'cont', 'NAdapt', 5);

%% first run to find initial fold
% zero problems
p0 = zeros(ncheb,1)+0.0;
p0(1:4) = [0.0;0.0;0.00;0.00];
par_args = cell(1,ncheb);   
for j=1:ncheb
    par_args{j} = sprintf('p%d',j);
end
xa        = [0.1;0;0;0]; %
bdata.x0  = xa;
[t0, x0]  = ode45(@(t,x) optcont(t, x, p0), -1:1/400:1, xa);
coll_args = {@optcont, @optcont_dx, @optcont_dp, @optcont_dt, t0, x0, par_args, p0};
prob1 = ode_isol2coll(prob, '', coll_args{:});
[data, uidx] = coco_get_func_data(prob1, 'coll', 'data', 'uidx');
maps = data.coll_seg.maps;
prob1 = coco_add_func(prob1, 'bc', @optcont_bc, @optcont_bc_du, ...
  @optcont_bc_dudu, bdata, 'zero', 'uidx', ...
  uidx([maps.x0_idx; maps.T0_idx; maps.T_idx]));

% objective
syms xx
TT = zeros(ncheb);
for i=1:ncheb
    for j=1:ncheb
        TT(i,j) = int(mychebyshevT(i-1,xx)*mychebyshevT(j-1,xx),xx,[-1,1]);
    end
end

data  = obj_init_data(data);
prob1 = coco_add_func(prob1, 'obj', @optcont_obj, @optcont_obj_du, data,...
    'inactive', 'obj', 'uidx', uidx, 'remesh', @obj_remesh);

data  = out_init_data(data);
prob1 = coco_add_func(prob1, 'output', @optcont_output, @optcont_output_du,...
    data, 'inequality', 'out', 'uidx', uidx, 'remesh', @output_remesh);

% adjoints
prob1 = adjt_isol2coll(prob1, '');
[data, axidx] = coco_get_adjt_data(prob1, 'coll', 'data', 'axidx');
opt = data.coll_opt;
prob1 = coco_add_adjt(prob1, 'bc', 'aidx', ...
  axidx([opt.x0_idx; opt.T0_idx; opt.T_idx]));

data  = adj_obj_init_data(data);
prob1 = coco_add_adjt(prob1, 'obj', @adj_optcont_obj, @adj_optcont_obj_du, data, 'd.obj',...
    'aidx', axidx([opt.xcn_idx; opt.T0_idx; opt.T_idx; opt.p_idx]),...
    'remesh', @adj_obj_remesh);

data  = adj_out_init_data(data);
prob1 = coco_add_adjt(prob1, 'output', @adj_optcont_output,...
    @adj_optcont_output_du, data, 'ncp.output',...
    'aidx', axidx([opt.xcn_idx; opt.T0_idx; opt.T_idx; opt.p_idx]),...
    'remesh', @adj_output_remesh);

cont_args = cell(1,ncheb+4);
cont_args{1} = 'obj';
cont_args{2} = 'd.obj';
cont_args{3} = sprintf('p%d',1);
cont_args{4} = sprintf('p%d',2);
for kid=1:ncheb-2
    cont_args{4+kid} = sprintf('d.p%d',kid+2);
end

cont_args{ncheb+3} = 'ncp.output';
cont_args{ncheb+4} = 'out';

coco(prob1, 'optcont1', [], 1, cont_args, [0 4000]);


prob = coco_prob;
prob = coco_set(prob, 'ode', 'autonomous', false);
prob = coco_set(prob, 'cont', 'h_max', 500);
prob = coco_set(prob, 'cont', 'MaxRes', 1.5);
prob = coco_set(prob, 'cont', 'al_max', 15);
prob = coco_set(prob, 'cont', 'NAdapt', 5);


%% branch switch from fold to grow nontrivial adjoint
bd1   = coco_bd_read('optcont1');
BPlab = coco_bd_labs(bd1, 'BP');
lab   = BPlab(2);

% zero problems
prob2 = ode_BP2coll(prob, '', 'optcont1', lab); 

[data, uidx] = coco_get_func_data(prob2, 'coll', 'data', 'uidx');
maps = data.coll_seg.maps;

prob2 = coco_add_func(prob2, 'bc', @optcont_bc, @optcont_bc_du, ...
  @optcont_bc_dudu, bdata, 'zero', 'uidx', ...
  uidx([maps.x0_idx; maps.T0_idx; maps.T_idx])); 

data  = obj_init_data(data);
prob2 = coco_add_func(prob2, 'obj', @optcont_obj, @optcont_obj_du, data,...
    'inactive', 'obj', 'uidx', uidx, 'remesh', @obj_remesh);

data  = out_init_data(data);
prob2 = coco_add_func(prob2, 'output', @optcont_output, @optcont_output_du,...
    data, 'inequality', 'out', 'uidx', uidx, 'remesh', @output_remesh);

% branch switch data
chart = coco_read_solution('optcont1', lab, 'chart');
cdata = coco_get_chart_data(chart, 'lsol');

% adjoints
prob2 = adjt_BP2coll(prob2, '', 'optcont1', lab);
[data, axidx] = coco_get_adjt_data(prob2, 'coll', 'data', 'axidx');
opt = data.coll_opt;

[chart, lidx] = coco_read_adjoint('bc', 'optcont1', lab, 'chart', 'lidx');
prob2   = coco_add_adjt(prob2, 'bc', 'aidx', ...
  axidx([opt.x0_idx; opt.T0_idx; opt.T_idx]), ...
  'l0', chart.x, 'tl0', cdata.v(lidx));

[chart, lidx] = coco_read_adjoint('obj', 'optcont1', lab, 'chart', 'lidx');
data  = adj_obj_init_data(data);
prob2 = coco_add_adjt(prob2, 'obj', @adj_optcont_obj, @adj_optcont_obj_du, data, 'd.obj',...
    'aidx', axidx([opt.xcn_idx; opt.T0_idx; opt.T_idx; opt.p_idx]),...
    'l0', chart.x, 'tl0', cdata.v(lidx), 'remesh', @adj_obj_remesh);

[chart, lidx] = coco_read_adjoint('output', 'optcont1', lab, 'chart', 'lidx');
data  = adj_out_init_data(data);
prob2 = coco_add_adjt(prob2, 'output', @adj_optcont_output,...
    @adj_optcont_output_du, data, 'ncp.output',...
    'aidx', axidx([opt.xcn_idx; opt.T0_idx; opt.T_idx; opt.p_idx]),...
    'l0', chart.x, 'tl0', cdata.v(lidx), 'remesh', @adj_output_remesh);

% computational domain
dobj_int = [0 1.1];
prob2 = coco_add_event(prob2, 'opt', 'BP', 'd.obj', '>', 1);

ss = cont_args{1};
cont_args{1} = cont_args{2};
cont_args{2} = ss;
coco(prob2, 'optcont2', [], cont_args, dobj_int);


%% continue to let d.pk=0
for k=1:ncheb-2
    prun = sprintf('optcont%d',k+1); % previous run
    crun = sprintf('optcont%d',k+2); % current run
    bd2 = coco_bd_read(prun);
    lab = coco_bd_labs(bd2, 'opt');

    % zero problems
    prob3 = ode_coll2coll(prob, '', prun, lab); 

    [data, uidx] = coco_get_func_data(prob3, 'coll', 'data', 'uidx');
    maps = data.coll_seg.maps;

    prob3 = coco_add_func(prob3, 'bc', @optcont_bc, @optcont_bc_du, ...
      @optcont_bc_dudu, bdata, 'zero', 'uidx', ...
      uidx([maps.x0_idx; maps.T0_idx; maps.T_idx])); 
  
    data  = obj_init_data(data);
    prob3 = coco_add_func(prob3, 'obj', @optcont_obj, @optcont_obj_du, data,...
        'inactive', 'obj', 'uidx', uidx, 'remesh', @obj_remesh);

    data  = out_init_data(data);
    prob3 = coco_add_func(prob3, 'output', @optcont_output, @optcont_output_du,...
        data, 'inequality', 'out', 'uidx', uidx, 'remesh', @output_remesh);

    % adjoints
    prob3 = adjt_coll2coll(prob3, '', prun, lab);

    chart = coco_read_adjoint('bc', prun, lab, 'chart');
    [data, axidx] = coco_get_adjt_data(prob3, 'coll', 'data', 'axidx');
    opt = data.coll_opt;

    prob3   = coco_add_adjt(prob3, 'bc', 'aidx', ...
      axidx([opt.x0_idx; opt.T0_idx; opt.T_idx]), 'l0', chart.x);
  
    chart = coco_read_adjoint('obj', prun, lab, 'chart');
    data  = adj_obj_init_data(data);
    prob3 = coco_add_adjt(prob3, 'obj', @adj_optcont_obj, @adj_optcont_obj_du, data, 'd.obj',...
    'aidx', axidx([opt.xcn_idx; opt.T0_idx; opt.T_idx; opt.p_idx]), 'l0', chart.x,...
    'remesh', @adj_obj_remesh);

    chart = coco_read_adjoint('output', prun, lab, 'chart');
    data  = adj_out_init_data(data);
    prob3 = coco_add_adjt(prob3, 'output', @adj_optcont_output,...
        @adj_optcont_output_du, data, 'ncp.output',...
        'aidx', axidx([opt.xcn_idx; opt.T0_idx; opt.T_idx; opt.p_idx]),...
        'l0', chart.x, 'remesh', @adj_output_remesh);
    % computational domain
    chart = coco_read_adjoint('coll.pars', prun, lab, 'chart');
    dp10  = chart.x(k+2);
    dmu = sprintf('d.p%d',k+2);
    if dp10>0
      dp1_int = [-0.1 dp10];
      prob3 = coco_add_event(prob3, 'opt', 'BP', dmu, '<', 0);
    else
      dp1_int = [dp10 0.1];
      prob3 = coco_add_event(prob3, 'opt', 'BP', dmu, '>', 0);
    end

    cont_args = cell(1,ncheb+4);
    cont_args{1} = dmu;
    cont_args{2} = 'obj';
    for kid=1:k+2
        cont_args{2+kid} = sprintf('p%d',kid);
    end
    for kid=k+3:ncheb
        cont_args{2+kid} = sprintf('d.p%d',kid);
    end
    cont_args{ncheb+3} = 'ncp.output';
    cont_args{ncheb+4} = 'out';

    coco(prob3, crun, [], cont_args, dp1_int);
end


%% drive ncp.ineq to zero
k = k+1;
prun = sprintf('optcont%d',k+1); % previous run
crun = sprintf('optcont%d',k+2); % current run
bd2 = coco_bd_read(prun);
lab = coco_bd_labs(bd2, 'opt');

% zero problems
prob3 = ode_coll2coll(prob, '', prun, lab); 

[data, uidx] = coco_get_func_data(prob3, 'coll', 'data', 'uidx');
maps = data.coll_seg.maps;

prob3 = coco_add_func(prob3, 'bc', @optcont_bc, @optcont_bc_du, ...
  @optcont_bc_dudu, bdata, 'zero', 'uidx', ...
  uidx([maps.x0_idx; maps.T0_idx; maps.T_idx])); 

data  = obj_init_data(data);
prob3 = coco_add_func(prob3, 'obj', @optcont_obj, @optcont_obj_du, data,...
    'inactive', 'obj', 'uidx', uidx, 'remesh', @obj_remesh);

data  = out_init_data(data);
prob3 = coco_add_func(prob3, 'output', @optcont_output, @optcont_output_du,...
    data, 'inequality', 'out', 'uidx', uidx, 'remesh', @output_remesh);

% adjoints
prob3 = adjt_coll2coll(prob3, '', prun, lab);

chart = coco_read_adjoint('bc', prun, lab, 'chart');
[data, axidx] = coco_get_adjt_data(prob3, 'coll', 'data', 'axidx');
opt = data.coll_opt;

prob3   = coco_add_adjt(prob3, 'bc', 'aidx', ...
  axidx([opt.x0_idx; opt.T0_idx; opt.T_idx]), 'l0', chart.x);

chart = coco_read_adjoint('obj', prun, lab, 'chart');
data  = adj_obj_init_data(data);
prob3 = coco_add_adjt(prob3, 'obj', @adj_optcont_obj, @adj_optcont_obj_du, data, 'd.obj',...
'aidx', axidx([opt.xcn_idx; opt.T0_idx; opt.T_idx; opt.p_idx]), 'l0', chart.x,...
'remesh', @adj_obj_remesh);

chart = coco_read_adjoint('output', prun, lab, 'chart');
data  = adj_out_init_data(data);
prob3 = coco_add_adjt(prob3, 'output', @adj_optcont_output,...
    @adj_optcont_output_du, data, 'ncp.output',...
    'aidx', axidx([opt.xcn_idx; opt.T0_idx; opt.T_idx; opt.p_idx]),...
    'l0', chart.x, 'remesh', @adj_output_remesh);

prob3 = coco_add_event(prob3, 'opt', 'BP', 'ncp.output', '==', 0);

% computational domain

cont_args = cell(1,ncheb+3);
cont_args{1} = 'ncp.output';
cont_args{2} = 'obj';
for kid=1:k+2
    cont_args{2+kid} = sprintf('p%d',kid);
end
cont_args{ncheb+3} = 'out';


coco(prob3, crun, [], cont_args, [-0.1 65.9]);



%% plot optimal results
%% numerical result
% nterm = 10;
crun = sprintf('optcont%d',nterm+1);
bd6 = coco_bd_read(crun);
% control input
par_args = cell(1,nterm);   
for j=1:nterm
    par_args{j} = sprintf('p%d',j);
end
p = coco_bd_col(bd6, par_args);
idx = coco_bd_idxs(bd6, 'opt');
popt = p(:,idx);
tp = (0:0.01:tf)';     % physical domain [0,2]
t = (tp-tf/2)/tf*2;    % standard domain [-1,1]
u = control_expansion(t,popt);

figure(1)
plot(tp,u,'b:','LineWidth',2); hold on
set(gca,'LineWidth',1.2);
set(gca,'FontSize',14);
xlabel('t');
ylabel('u(t)');
% state variable
lab = coco_bd_labs(bd6, 'opt');
sol = coll_read_solution('', crun, lab);
x  = sol.xbp;
[t0, x0]  = ode45(@(t,x) optcont(t, x, zeros(nterm,1)), -1:0.1:1, xa);
figure(2)
plot(sol.tbp*tf/2+tf/2,x(:,1),'b:','LineWidth',2); hold on
xlabel('t');
ylabel('x(t)');
set(gca,'LineWidth',1.2);
set(gca,'FontSize',14);
xlabel('t');
ylabel('\theta(t)');
axis([0,2,-0.08,0.1])
figure(3)
plot(sol.tbp*tf/2+tf/2,x(:,3),'b:','LineWidth',2); hold on
set(gca,'LineWidth',1.2);
set(gca,'FontSize',14);
xlabel('t');
ylabel('x(t)')

fprintf('control cost:%d\n',popt'*TT*popt);

diary off
end

function y = optcont(t, x, p)

global ncheb tT
global M m g l
u = 0;
for j=1:ncheb
    u = u + p(j,:).*mychebyshevT(j-1,t);
end
x1 = x(1,:);
x2 = x(2,:);
x4 = x(4,:);
f1 = x2;
f2 = (u.*cos(x1)-(M+m)*g*sin(x1)+m*l*cos(x1).*sin(x1).*x2.^2)./(m*l*(cos(x1)).^2-(M+m)*l);
f3 = x4;
f4 = (u+m*l*x2.^2.*sin(x1)-m*g*sin(x1).*cos(x1))./(M+m-m*(cos(x1)).^2);

y = tT*[f1;f2;f3;f4];




end

function J = optcont_dx(t, x, p) 

global ncheb tT
global M m g l
x1 = x(1,:);
x2 = x(2,:);

J = zeros(4,4,numel(t));
u = 0;
for j=1:ncheb
    u = u + p(j,:).*mychebyshevT(j-1,t);
end
J(1,2,:) = 1;
J(2,1,:) = (u.*sin(x1) + g.*cos(x1).*(M + m) - l.*m.*x2.^2.*cos(x1).^2 + l.*m.*x2.^2.*sin(x1).^2)./(l.*(M + m) - l.*m.*cos(x1).^2) + (2.*l.*m.*cos(x1).*sin(x1).*(l.*m.*cos(x1).*sin(x1).*x2.^2 + u.*cos(x1) - g.*sin(x1).*(M + m)))./(l.*(M + m) - l.*m.*cos(x1).^2).^2;
J(2,2,:) = -(2.*l.*m.*x2.*cos(x1).*sin(x1))./(l.*(M + m) - l.*m.*cos(x1).^2);
J(3,4,:) = 1;
J(4,1,:) = (g.*m.*sin(x1).^2 - g.*m.*cos(x1).^2 + l.*m.*x2.^2.*cos(x1))./(M + m - m.*cos(x1).^2) - (2.*m.*cos(x1).*sin(x1).*(l.*m.*sin(x1).*x2.^2 + u - g.*m.*cos(x1).*sin(x1)))./(- m.*cos(x1).^2 + M + m).^2;
J(4,2,:) = (2.*l.*m.*x2.*sin(x1))./(M + m - m.*cos(x1).^2);

J = J*tT;

% % 
% f  = @(x,p) optcont(p(1,:), x, p(2:end,:));
% tp = [ t ; p ];
% Jx = coco_ezDFDX('f(x,p)v', f, x, tp);
% max(abs(J(:)-Jx(:)))
end

function J = optcont_dp(t, x, p) 

global ncheb tT
global M m l
x1 = x(1,:);
u  = 0;
for j=1:ncheb
    u = u + p(j,:).*mychebyshevT(j-1,t);
end
dy2du = cos(x1)./(m*l*(cos(x1)).^2-(M+m)*l);
dy4du = 1./(M+m-m*(cos(x1)).^2);
J = zeros(4,ncheb,numel(t));
for j=1:ncheb
    J(2,j,:) = dy2du.*mychebyshevT(j-1,t);
    J(4,j,:) = dy4du.*mychebyshevT(j-1,t);
end
J  = tT*J;

% %
% f  = @(x,p) optcont(x(1,:), x(2:end,:), p);
% tx = [ t ; x ];
% Jp = coco_ezDFDP('f(x,p)v', f, tx, p);
% max(abs(J(:)-Jp(:)))
end

function J = optcont_dt(t, x, p) 

global ncheb tT
global M m l
x1 = x(1,:);
u  = 0;
for j=1:ncheb
    u = u + p(j,:).*mychebyshevT(j-1,t);
end
dudt = 0;
for j=2:ncheb
    dudt = dudt + p(j,:).*dchebyshevT(j-1,t);
end
dy2du = cos(x1)./(m*l*(cos(x1)).^2-(M+m)*l);
dy4du = 1./(M+m-m*(cos(x1)).^2);

J = zeros(4,numel(t));
J(2,:) = dy2du.*dudt;
J(4,:) = dy4du.*dudt;
J  = tT*J;

% %
% [mm, n] = size(x);
% f  = @(x,p) optcont(x(1,:), p(1:mm,:), p(mm+1:end,:));
% xp = [ x ; p ];
% Jt = reshape(coco_ezDFDX('f(x,p)v', f, t, xp), [mm n]);
% max(abs(J(:)-Jt(:)))
end



function [data, y] = optcont_bc(prob, data, u) %#ok<INUSL>

x0 = u(1:4);
T0 = u(5);
T  = u(6);

y = [x0-data.x0; T0+1; T-2];

end

function [data, J] = optcont_bc_du(prob, data, u) %#ok<INUSD,INUSL>

J = eye(6);

end

function [data, dJ] = optcont_bc_dudu(prob, data, u) %#ok<INUSD,INUSL>

dJ = zeros(6,6,6);

end





function data = out_init_data(fdata)

data.coll_seg = fdata.coll_seg;
data.yhan     = @yhan;
data.yhan_dx  = @yhan_dx;

data = coco_func_data(data);

end



function [prob, status, xtr] = output_remesh(prob, data, chart, old_u, old_V) %#ok<INUSD>

[fdata, uidx] = coco_get_func_data(prob, 'coll', 'data', 'uidx');
data.coll_seg = fdata.coll_seg;

xtr    = [];
prob   = coco_change_func(prob, data, 'uidx', uidx);
status = 'success';

end



function [data, y] = optcont_output(prob, data, u)

global tT

pr   = data.pr;
maps = pr.coll_seg.maps;
mesh = pr.coll_seg.mesh;

T    = u(maps.T_idx);
x    = u(maps.xbp_idx);
p    = u(maps.p_idx);

xcn = reshape(maps.W*x, maps.x_shp);
gcn = pr.yhan(xcn);
gcn = mesh.gka.*gcn;   % mesh here is related adaptive meshing

y = (0.5*T/maps.NTST)*mesh.gwt*gcn';
y = tT*y-0.01;

end


function [data, J] = optcont_output_du(prob, data, u)

global tT ncheb


pr   = data.pr;
maps = pr.coll_seg.maps;
mesh = pr.coll_seg.mesh;

T = u(maps.T_idx);
x = u(maps.xbp_idx);
p = u(maps.p_idx);

xcn = reshape(maps.W*x, maps.x_shp);
gcn = pr.yhan(xcn);
gcn = mesh.gka.*gcn;

gdxcn = pr.yhan_dx(xcn);
gdxcn = mesh.gdxka.*gdxcn;
gdxcn = sparse(maps.gdxrows, maps.gdxcols, gdxcn(:));

J_xbp = (0.5*T/maps.NTST)*mesh.gwt*gdxcn*maps.W;
J_T0  = 0;
J_T   = (0.5/maps.NTST)*mesh.gwt*gcn';
J_p   = zeros(1,ncheb);
J     = [J_xbp J_T0 J_T J_p];
J     = tT*J;

end


function data = adj_out_init_data(fdata)

data.coll_seg  = fdata.coll_seg;
data.yhan      = @yhan;
data.yhan_dx   = @yhan_dx;
data.yhan_dxdx = @yhan_dxdx;

seg  = fdata.coll_seg;
maps = seg.maps;
int  = seg.int;

NCOL = int.NCOL;
NTST = maps.NTST;
xdim = int.dim;
pdim = maps.pdim;

rows = NCOL*NTST*kron(0:(xdim-1), ones(1,xdim));
opt.gdxdxrows1 = repmat(rows, [1 NCOL*NTST]) + ...
  kron(1:NCOL*NTST, ones(1,xdim^2));
cols = reshape(1:xdim*NCOL*NTST, [xdim NCOL*NTST]);
opt.gdxdxcols1 = repmat(cols, [xdim 1]);

% Derivative of (T/2N)*gxcn with respect to xbp:
step = 1+xdim*(0:NCOL-1); % NCOL
step = repmat(step(:), [1 xdim]) + repmat(0:xdim-1, [NCOL 1]); % xdim*NCOL
step = repmat(step(:), [1 xdim*(NCOL+1)]) + ...
  (xdim*NCOL*NTST+2+pdim)*repmat(0:xdim*(NCOL+1)-1, [xdim*NCOL 1]); % xdim^2*NCOL*(NCOL+1)
step = repmat(step(:), [1 NTST]) + ...
  (xdim*NCOL+xdim*(NCOL+1)*(xdim*NCOL*NTST+2+pdim))*...
  repmat(0:NTST-1, [xdim^2*NCOL*(NCOL+1) 1]); % xdim^2*NCOL*(NCOL+1)*NTST
opt.gdxdxcols2 = step(:);
opt.gdxdxrows2 = ones(xdim^2*NCOL*(NCOL+1)*NTST, 1);

step = 1:NCOL; % NCOL
step = repmat(step(:), [1 xdim]) + NTST*NCOL*repmat(0:xdim-1, [NCOL 1]); % xdim*NCOL
step = repmat(step(:), [1 xdim*(NCOL+1)]) + ...
  xdim*NTST*NCOL*repmat(0:xdim*(NCOL+1)-1, [xdim*NCOL 1]); % xdim^2*NCOL*(NCOL+1)
step = repmat(step(:), [1 NTST]) + (NCOL+xdim^2*NTST*NCOL*(NCOL+1))*...
  repmat(0:NTST-1, [xdim^2*NCOL*(NCOL+1) 1]); % xdim^2*NCOL*(NCOL+1)*NTST
opt.gdxdxidx = step(:);

% Derivative of (T/2N)*gxcn with respect to T:
opt.gdxdTrows = ones(xdim*NTST*NCOL, 1);
opt.gdxdTcols = (xdim*NCOL*NTST+2+pdim)*(xdim*(NCOL+1)*NTST+1) + ...
  (1:xdim*NTST*NCOL)';

% Derivative of (1/2N)*w*g' with respect to xbp:
opt.gdTdxcols = NTST*NCOL*xdim+2 + ...
  (xdim*NTST*NCOL+2+pdim)*(0:xdim*(NCOL+1)*NTST-1)';
opt.gdTdxrows = ones(xdim*(NCOL+1)*NTST, 1);

% Derivative of [p1/5, p2/5, p3/5] with respect to p:

opt.gdpdprows = ones(pdim*pdim,1);
% opt.gdpdpcols = (NTST*xdim*NCOL+2+pdim)*(xdim*(NCOL+1)*NTST+2) + ...
%   xdim*NTST*NCOL+2+[1 NTST*xdim*NCOL+2+pdim+2 2*(NTST*xdim*NCOL+2+pdim)+3]';
gdpdpcols = (NTST*xdim*NCOL+2+pdim)*kron(0:(pdim-1), ones(1,pdim));
opt.gdpdpcols = (NTST*xdim*NCOL+2+pdim)*(xdim*(NCOL+1)*NTST+2)+...
    xdim*NTST*NCOL+2+gdpdpcols+repmat((1:pdim),[1,pdim]);

opt.dJrows = 1;
opt.dJcols = (xdim*NTST*NCOL+2+pdim)*(xdim*NTST*(NCOL+1)+2+pdim);

data.coll_opt = opt;

data = coco_func_data(data);

end

function [prob, status, xtr, ftr] = adj_output_remesh(prob, data, chart, lb, Vlb)  %#ok<INUSL>

[fdata, axidx] = coco_get_adjt_data(prob, 'coll', 'data', 'axidx');
data = adj_out_init_data(fdata);
fopt = fdata.coll_opt;

aidx = axidx([fopt.xcn_idx; fopt.T0_idx; fopt.T_idx; fopt.p_idx]);

xtr    = [];
ftr    = 1;
prob   = coco_change_adjt(prob, data, 'aidx', aidx, 'l0', lb, 'vecs', Vlb);
status = 'success';

end




function [data, y] = adj_optcont_output(prob, data, u)

global tT ncheb

pr   = data.pr;
maps = pr.coll_seg.maps;
mesh = pr.coll_seg.mesh;

T = u(maps.T_idx);
x = u(maps.xbp_idx);
p = u(maps.p_idx);

xcn = reshape(maps.W*x, maps.x_shp);
gcn = pr.yhan(xcn);
gcn = mesh.gka.*gcn;

gdxcn = pr.yhan_dx(xcn);
gdxcn = mesh.gdxka.*gdxcn;

J_xbp = (0.5*T/maps.NTST)*gdxcn(:)';
J_T0  = 0;
J_T   = (0.5/maps.NTST)*mesh.gwt*gcn';
J_p   = zeros(1,ncheb);

y = [J_xbp J_T0 J_T J_p];
y = tT*y;

end


function [data, J] = adj_optcont_output_du(prob, data, u)

global tT ncheb

pr   = data.pr;
maps = pr.coll_seg.maps;
mesh = pr.coll_seg.mesh;
opt  = pr.coll_opt;

T = u(maps.T_idx);
x = u(maps.xbp_idx);

xcn = reshape(maps.W*x, maps.x_shp);

gdxdxcn = pr.yhan_dxdx(xcn);
gdxdxcn = mesh.gdxdxka.*gdxdxcn;
gdxdxcn = sparse(opt.gdxdxrows1, opt.gdxdxcols1, gdxdxcn(:))*maps.W;
J       = (0.5*T/maps.NTST)*sparse(opt.gdxdxrows2, opt.gdxdxcols2, ...
  gdxdxcn(opt.gdxdxidx), opt.dJrows, opt.dJcols);

gdxcn   = pr.yhan_dx(xcn);
gdxcn   = mesh.gdxka.*gdxcn;
J       = J + (0.5/maps.NTST)*sparse(opt.gdxdTrows, opt.gdxdTcols, ...
  gdxcn(:), opt.dJrows, opt.dJcols);

gdxcn   = mesh.gwt*sparse(maps.gdxrows, maps.gdxcols, gdxcn(:))*maps.W;
J       = J + (0.5/maps.NTST)*sparse(opt.gdTdxrows, opt.gdTdxcols, ...
  gdxcn(:), opt.dJrows, opt.dJcols);

J = tT*J;
% [data, Jd] = coco_ezDFDX('f(o,d,x)', prob, data, @adj_optcont_output, u);
% max(max(abs(J-Jd(:)')))

end



function data = obj_init_data(fdata)

data.coll_seg = fdata.coll_seg;
data.ghan     = @ghan;
data.ghan_dx  = @ghan_dx;

data = coco_func_data(data);

end

function [prob, status, xtr] = obj_remesh(prob, data, chart, old_u, old_V) %#ok<INUSD>

[fdata, uidx] = coco_get_func_data(prob, 'coll', 'data', 'uidx');
data.coll_seg = fdata.coll_seg;

xtr    = [];
prob   = coco_change_func(prob, data, 'uidx', uidx);
status = 'success';

end


function [data, y] = optcont_obj(prob, data, u)

global tT TT

pr   = data.pr;
maps = pr.coll_seg.maps;
mesh = pr.coll_seg.mesh;

T    = u(maps.T_idx);
x    = u(maps.xbp_idx);
p    = u(maps.p_idx);

xcn = reshape(maps.W*x, maps.x_shp);
gcn = pr.ghan(xcn);
gcn = mesh.gka.*gcn;   % mesh here is related adaptive meshing

y   = (0.5*T/maps.NTST)*mesh.gwt*gcn' + p'*TT*p;
y   = tT*y;



end


function [data, J] = optcont_obj_du(prob, data, u)

global tT TT

pr   = data.pr;
maps = pr.coll_seg.maps;
mesh = pr.coll_seg.mesh;

T = u(maps.T_idx);
x = u(maps.xbp_idx);
p = u(maps.p_idx);

xcn = reshape(maps.W*x, maps.x_shp);
gcn = pr.ghan(xcn);
gcn = mesh.gka.*gcn;

gdxcn = pr.ghan_dx(xcn);
gdxcn = mesh.gdxka.*gdxcn;
gdxcn = sparse(maps.gdxrows, maps.gdxcols, gdxcn(:));

J_xbp = (0.5*T/maps.NTST)*mesh.gwt*gdxcn*maps.W;
J_T0  = 0;
J_T   = (0.5/maps.NTST)*mesh.gwt*gcn';
J_p   = 2*TT*p;
J_p   = J_p';

J = [J_xbp J_T0 J_T J_p];
J = tT*J;

% [data, Jd] = coco_ezDFDX('f(o,d,x)', prob, data, @optcont_obj, u);
% max(max(abs(J-Jd)))
end



function data = adj_obj_init_data(fdata)

data.coll_seg  = fdata.coll_seg;
data.ghan      = @ghan;
data.ghan_dx   = @ghan_dx;
data.ghan_dxdx = @ghan_dxdx;

seg  = fdata.coll_seg;
maps = seg.maps;
int  = seg.int;

NCOL = int.NCOL;
NTST = maps.NTST;
xdim = int.dim;
pdim = maps.pdim;

rows = NCOL*NTST*kron(0:(xdim-1), ones(1,xdim));
opt.gdxdxrows1 = repmat(rows, [1 NCOL*NTST]) + ...
  kron(1:NCOL*NTST, ones(1,xdim^2));
cols = reshape(1:xdim*NCOL*NTST, [xdim NCOL*NTST]);
opt.gdxdxcols1 = repmat(cols, [xdim 1]);

% Derivative of (T/2N)*gxcn with respect to xbp:
step = 1+xdim*(0:NCOL-1); % NCOL
step = repmat(step(:), [1 xdim]) + repmat(0:xdim-1, [NCOL 1]); % xdim*NCOL
step = repmat(step(:), [1 xdim*(NCOL+1)]) + ...
  (xdim*NCOL*NTST+2+pdim)*repmat(0:xdim*(NCOL+1)-1, [xdim*NCOL 1]); % xdim^2*NCOL*(NCOL+1)
step = repmat(step(:), [1 NTST]) + ...
  (xdim*NCOL+xdim*(NCOL+1)*(xdim*NCOL*NTST+2+pdim))*...
  repmat(0:NTST-1, [xdim^2*NCOL*(NCOL+1) 1]); % xdim^2*NCOL*(NCOL+1)*NTST
opt.gdxdxcols2 = step(:);
opt.gdxdxrows2 = ones(xdim^2*NCOL*(NCOL+1)*NTST, 1);

step = 1:NCOL; % NCOL
step = repmat(step(:), [1 xdim]) + NTST*NCOL*repmat(0:xdim-1, [NCOL 1]); % xdim*NCOL
step = repmat(step(:), [1 xdim*(NCOL+1)]) + ...
  xdim*NTST*NCOL*repmat(0:xdim*(NCOL+1)-1, [xdim*NCOL 1]); % xdim^2*NCOL*(NCOL+1)
step = repmat(step(:), [1 NTST]) + (NCOL+xdim^2*NTST*NCOL*(NCOL+1))*...
  repmat(0:NTST-1, [xdim^2*NCOL*(NCOL+1) 1]); % xdim^2*NCOL*(NCOL+1)*NTST
opt.gdxdxidx = step(:);

% Derivative of (T/2N)*gxcn with respect to T:
opt.gdxdTrows = ones(xdim*NTST*NCOL, 1);
opt.gdxdTcols = (xdim*NCOL*NTST+2+pdim)*(xdim*(NCOL+1)*NTST+1) + ...
  (1:xdim*NTST*NCOL)';

% Derivative of (1/2N)*w*g' with respect to xbp:
opt.gdTdxcols = NTST*NCOL*xdim+2 + ...
  (xdim*NTST*NCOL+2+pdim)*(0:xdim*(NCOL+1)*NTST-1)';
opt.gdTdxrows = ones(xdim*(NCOL+1)*NTST, 1);

% Derivative of 2p^\top TT with respect to p:

opt.gdpdprows = ones(pdim*pdim,1);
% opt.gdpdpcols = (NTST*xdim*NCOL+2+pdim)*(xdim*(NCOL+1)*NTST+2) + ...
%   xdim*NTST*NCOL+2+[1 NTST*xdim*NCOL+2+pdim+2 2*(NTST*xdim*NCOL+2+pdim)+3]';
gdpdpcols = (NTST*xdim*NCOL+2+pdim)*kron(0:(pdim-1), ones(1,pdim));
opt.gdpdpcols = (NTST*xdim*NCOL+2+pdim)*(xdim*(NCOL+1)*NTST+2)+...
    xdim*NTST*NCOL+2+gdpdpcols+repmat((1:pdim),[1,pdim]);


opt.dJrows = 1;
opt.dJcols = (xdim*NTST*NCOL+2+pdim)*(xdim*NTST*(NCOL+1)+2+pdim);

data.coll_opt = opt;

data = coco_func_data(data);

end

function [prob, status, xtr, ftr] = adj_obj_remesh(prob, data, chart, lb, Vlb)  %#ok<INUSL>

[fdata, axidx] = coco_get_adjt_data(prob, 'coll', 'data', 'axidx');
data = adj_obj_init_data(fdata);
fopt = fdata.coll_opt;

aidx = axidx([fopt.xcn_idx; fopt.T0_idx; fopt.T_idx; fopt.p_idx]);

xtr    = [];
ftr    = 1;
prob   = coco_change_adjt(prob, data, 'aidx', aidx, 'l0', lb, 'vecs', Vlb);
status = 'success';

end




function [data, y] = adj_optcont_obj(prob, data, u)

global tT TT

pr   = data.pr;
maps = pr.coll_seg.maps;
mesh = pr.coll_seg.mesh;

T = u(maps.T_idx);
x = u(maps.xbp_idx);
p = u(maps.p_idx);

xcn = reshape(maps.W*x, maps.x_shp);
gcn = pr.ghan(xcn);
gcn = mesh.gka.*gcn;

gdxcn = pr.ghan_dx(xcn);
gdxcn = mesh.gdxka.*gdxcn;

J_xbp = (0.5*T/maps.NTST)*gdxcn(:)';
J_T0  = 0;
J_T   = (0.5/maps.NTST)*mesh.gwt*gcn';
J_p   = 2*TT*p;
J_p   = J_p';

y = [J_xbp J_T0 J_T J_p];
y = tT*y;

end


function [data, J] = adj_optcont_obj_du(prob, data, u)

global tT TT


pr   = data.pr;
maps = pr.coll_seg.maps;
mesh = pr.coll_seg.mesh;
opt  = pr.coll_opt;

T = u(maps.T_idx);
x = u(maps.xbp_idx);

xcn = reshape(maps.W*x, maps.x_shp);

gdxdxcn = pr.ghan_dxdx(xcn);
gdxdxcn = mesh.gdxdxka.*gdxdxcn;
gdxdxcn = sparse(opt.gdxdxrows1, opt.gdxdxcols1, gdxdxcn(:))*maps.W;
J       = (0.5*T/maps.NTST)*sparse(opt.gdxdxrows2, opt.gdxdxcols2, ...
  gdxdxcn(opt.gdxdxidx), opt.dJrows, opt.dJcols);

gdxcn   = pr.ghan_dx(xcn);
gdxcn   = mesh.gdxka.*gdxcn;
J       = J + (0.5/maps.NTST)*sparse(opt.gdxdTrows, opt.gdxdTcols, ...
  gdxcn(:), opt.dJrows, opt.dJcols);

gdxcn   = mesh.gwt*sparse(maps.gdxrows, maps.gdxcols, gdxcn(:))*maps.W;
J       = J + (0.5/maps.NTST)*sparse(opt.gdTdxrows, opt.gdTdxcols, ...
  gdxcn(:), opt.dJrows, opt.dJcols);

J       = J + sparse(opt.gdpdprows, opt.gdpdpcols, 2*TT', ...
  opt.dJrows, opt.dJcols);


J = tT*J;
% [data, Jd] = coco_ezDFDX('f(o,d,x)', prob, data, @adj_optcont_obj, u);
% max(max(abs(J-Jd(:)')))

end




function y = ghan(x)

x1 = x(1,:);
x3 = x(3,:);
y  = 100*x1.^2+40*x3.^2;

end


function J = ghan_dx(x)

x1 = x(1,:);
x3 = x(3,:);
J  = zeros(1, 4, numel(x1));
J(1,1,:) = 200*x1;
J(1,3,:) = 80*x3;

end


function J = ghan_dxdx(x)

x1 = x(1,:);
J  = zeros(1, 4, 4, numel(x1));
J(1,1,1,:) = 200;
J(1,3,3,:) = 80;

end





function y = yhan(x)

x1 = x(1,:);
x3 = x(3,:);
y  = x1.^2+x3.^2;

end


function J = yhan_dx(x)

x1 = x(1,:);
x3 = x(3,:);
J  = zeros(1, 4, numel(x1));
J(1,1,:) = 2*x1;
J(1,3,:) = 2*x3;

end


function J = yhan_dxdx(x)

x1 = x(1,:);
J  = zeros(1, 4, 4, numel(x1));
J(1,1,1,:) = 2;
J(1,3,3,:) = 2;

end