

%% Construction initial solutions for two periodic orbits and connected trajectory
global ncheb tT % tT =(tf-t0)/2
ncheb = 10;
tT = 1;
%% generate two periodic orbits with given period T=2
% get initial one from ref.15
xa = [-0.252566564360231;0;0.273165147711868;0;1.254563849725983;0];
p0 = 0;
tf = 179.5443;
n1 = 1.99099e-7;
tf = tf*n1*3600*24;
options  = odeset('reltol',1e-13,'abstol',1e-15);
[t0,x0] = ode45(@(t,x) RHCEquation_orbit(t,x,p0),[0,tf],xa,options);

prob = coco_prob();
prob = coco_set(prob,'coll','NTST',20);
prob = coco_set(prob, 'cont', 'NAdapt', 5, 'h_max', 10);

funcs = { @cr3bp, @cr3bp_dx, @cr3bp_dp};
coll_args = { funcs{:}, t0, x0, p0}; % Regarding {'alpha1'} as continuation variables
prob1 = ode_isol2coll(prob, 'orbit1', coll_args{:});
[data, uidx] = coco_get_func_data(prob1, 'orbit1.coll', 'data', 'uidx');
maps = data.coll_seg.maps;
cr3bp_bc_funcs = {@cr3bp_bc, @cr3bp_bc_du, @cr3bp_bc_dudu};
prob1 = coco_add_func(prob1, 'orbit1_bc', cr3bp_bc_funcs{:}, [], 'zero', 'uidx', [uidx(maps.x0_idx);uidx(maps.x1_idx)]);
prob1 = coco_add_pars(prob1, 'orbit1_par_x0', uidx(maps.x0_idx(2)),'x0_1');
prob1 = coco_add_pars(prob1, 'orbit1_par_T', uidx(maps.T_idx),'T_1');
prob1 = coco_add_event(prob1, 'T2', 'BP', 'T_1', '=', 2);
cont_args = { 1, {'T_1'} };
bd01 = coco(prob1, 'Halo1', [], cont_args{:});

% extract the periodic orbit
lab = coco_bd_labs(bd01, 'T2');
lab1 = lab(1);
sol = coll_read_solution('orbit1', 'Halo1', lab1);
xp1 = sol.xbp;
tp1 = sol.tbp;

%% construct composite problem to generate H3 and connected trajectory
% add two instances of periodic orbits
prob = coco_prob(); % use this to deal with coexist of autonomous and non-autonomous trajectories,
prob = coco_set(prob,'coll','NTST',50);
prob1 = coco_set(prob, 'cont', 'NAdapt', 5, 'h_max', 10);

% construct two halo orbits
% reconstruct the periodic orbit with T=2 (H2)
prob1 = ode_coll2coll(prob1, 'halo1', 'Halo1', 'orbit1', lab1);
[data_halo1, uidx_halo1] = coco_get_func_data(prob1, 'halo1.coll', 'data', 'uidx');
maps_halo1 = data_halo1.coll_seg.maps;
prob1 = coco_add_func(prob1, 'halo1_bc', cr3bp_bc_funcs{:}, [], 'zero', 'uidx', [uidx_halo1(maps_halo1.x0_idx);uidx_halo1(maps_halo1.x1_idx)]);
prob1 = coco_add_pars(prob1, 'halo1_par_x0', uidx_halo1(maps_halo1.x0_idx(2)),'x0_1');
prob1 = coco_add_pars(prob1, 'halo1_par_T', uidx_halo1(maps_halo1.T_idx),'T_1');

% generate another periodic orbit with T=2 (initial H3)
coll_args = { funcs{:}, tp1, xp1, p0}; % Regarding {'alpha2'} as continuatino variables
prob1 = ode_isol2coll(prob1, 'halo2', coll_args{:});
[data_halo2, uidx_halo2] = coco_get_func_data(prob1, 'halo2.coll', 'data', 'uidx');
maps_halo2 = data_halo2.coll_seg.maps;
prob1 = coco_add_func(prob1, 'halo2_bc', cr3bp_bc_funcs{:}, [], 'zero', 'uidx', [uidx_halo2(maps_halo2.x0_idx);uidx_halo2(maps_halo2.x1_idx)]);
prob1 = coco_add_pars(prob1, 'halo2_par_x0', uidx_halo2(maps_halo2.x0_idx(2)),'x0_2');
prob1 = coco_add_pars(prob1, 'halo2_par_T', uidx_halo2(maps_halo2.T_idx),'T_2');

% generate connected trajectory in terms of periodic solutions
p0 = ones(3*ncheb,1)*0.0000;
t0 = 2*tp1/(tp1(end)-tp1(1))-1;
x0 = xp1;
par_args = cell(1,3*ncheb);   
for i=1:3
    for j=1:ncheb
        par_args{ncheb*(i-1)+j} = sprintf('p%d',ncheb*(i-1)+j);
    end
end
coll_args = {@optcont, @optcont_dx, @optcont_dp, @optcont_dt,...
    t0, x0, par_args, p0};
prob1 = coco_set(prob1, 'connect.ode', 'autonomous', false); % set the transfer trajectory to be non-autonomous
prob1 = ode_isol2coll(prob1, 'connect', coll_args{:});
[data_connect, uidx_connect] = coco_get_func_data(prob1, 'connect.coll', 'data', 'uidx');
maps_connect = data_connect.coll_seg.maps;
% impose conditions for T0 and T
T_bc_funcs = {@T_bc, @T_bc_du, @T_bc_dudu};
prob1 = coco_add_func(prob1, 'T_bc', T_bc_funcs{:}, [], 'zero', ...
  'uidx', [uidx_connect(maps_connect.T0_idx);uidx_connect(maps_connect.T_idx)]);

% add links between these instances
[data_connect, uidx_connect] = coco_get_func_data(prob1, 'connect.coll', 'data', 'uidx');
maps_connect = data_connect.coll_seg.maps;
cont_bc_funcs = {@optcont_bc, @optcont_bc_du, @optcont_bc_dudu};
prob1 = coco_add_func(prob1, 'cont_bc', cont_bc_funcs{:}, [], 'zero', ...
  'uidx', [uidx_halo1(maps_halo1.x0_idx); uidx_connect(maps_connect.x0_idx);...
  uidx_connect(maps_connect.x1_idx); uidx_halo2(maps_halo2.x0_idx)]);

prob1 = coco_add_event(prob1, 'T3', 'T_2', '=', 3);
cont_args = { 1, {'T_2' 'p1' 'p11' 'p21' 'p2' 'p12' 'p22'} };
bd = coco(prob1, 'HALO3', [], cont_args{:}, [1.8 3.2]);

% plot continuation process
% klab = [3 6 9 12];
% for k=klab
%     sol1 = coll_read_solution('halo1','HALO3', k);
%     sol2 = coll_read_solution('halo2','HALO3', k);
%     sol3 = coll_read_solution('connect','HALO3', k);
%     plot3(sol1.xbp(:,1),sol1.xbp(:,2),sol1.xbp(:,3),'b-','LineWidth',2); hold on
%     plot3(sol2.xbp(:,1),sol2.xbp(:,2),sol2.xbp(:,3),'b-','LineWidth',2);
%     plot3(sol3.xbp(:,1),sol3.xbp(:,2),sol3.xbp(:,3),'b-','LineWidth',2);
%     pause;
% end
% set(gca,'FontSize', 14);
% set(gca,'LineWidth',1.5);
% xlabel('$x$','interpreter','latex','FontSize',18);
% ylabel('$y$','interpreter','latex','FontSize',18);
% zlabel('$z$','interpreter','latex','FontSize',18);
lab = coco_bd_labs(bd, 'T3');
lab = max(lab); 
% 2nd is in the other family, then we have North-South transfer
% sol = coll_read_solution('halo2','HALO3', lab);



%% Construct optimization problems
% ========
% original system
% add two instances of periodic orbits

prob = coco_prob(); % use this to deal with coexist of autonomous and non-autonomous trajectories,
prob1 = coco_set(prob, 'cont', 'NAdapt', 5, 'h_max', 10);

% construct two halo orbits
prob1 = ode_coll2coll(prob1, 'halo1', 'HALO3', 'halo1', lab);
prob1 = ode_coll2coll(prob1, 'halo2', 'HALO3', 'halo2', lab);
% extract data
[data_halo1, uidx_halo1] = coco_get_func_data(prob1, 'halo1.coll', 'data', 'uidx');
maps_halo1 = data_halo1.coll_seg.maps;
[data_halo2, uidx_halo2] = coco_get_func_data(prob1, 'halo2.coll', 'data', 'uidx');
maps_halo2 = data_halo2.coll_seg.maps;
% impose periodic BCs for two halo orbits
cr3bp_bc_funcs = {@cr3bp_bc, @cr3bp_bc_du, @cr3bp_bc_dudu};
prob1 = coco_add_func(prob1, 'halo1_bc', cr3bp_bc_funcs{:}, [], 'zero', 'uidx', [uidx_halo1(maps_halo1.x0_idx);uidx_halo1(maps_halo1.x1_idx)]);
prob1 = coco_add_func(prob1, 'halo2_bc', cr3bp_bc_funcs{:}, [], 'zero', 'uidx', [uidx_halo2(maps_halo2.x0_idx);uidx_halo2(maps_halo2.x1_idx)]);
% Define inactive prameters
prob1 = coco_add_pars(prob1, 'halo1_par_x0', uidx_halo1(maps_halo1.x0_idx(2)),'x0_1');
prob1 = coco_add_pars(prob1, 'halo1_par_T', uidx_halo1(maps_halo1.T_idx),'T_1');
prob1 = coco_add_pars(prob1, 'halo2_par_x0', uidx_halo2(maps_halo2.x0_idx(2)),'x0_2');
prob1 = coco_add_pars(prob1, 'halo2_par_T', uidx_halo2(maps_halo2.T_idx),'T_2');

% add one instance of connected trajectory
prob1 = ode_coll2coll(prob1, 'connect', 'HALO3', 'connect', lab);
[data_connect, uidx_connect] = coco_get_func_data(prob1, 'connect.coll', 'data', 'uidx');
maps_connect = data_connect.coll_seg.maps;
% impose conditions for T0 and T
T_bc_funcs = {@T_bc, @T_bc_du, @T_bc_dudu};
prob1 = coco_add_func(prob1, 'T_bc', T_bc_funcs{:}, [], 'zero', ...
  'uidx', [uidx_connect(maps_connect.T0_idx);uidx_connect(maps_connect.T_idx)]);
% add links between these instances
[data_connect, uidx_connect] = coco_get_func_data(prob1, 'connect.coll', 'data', 'uidx');
maps_connect = data_connect.coll_seg.maps;
cont_bc_funcs = {@optcont_bc, @optcont_bc_du, @optcont_bc_dudu};
prob1 = coco_add_func(prob1, 'cont_bc', cont_bc_funcs{:}, [], 'zero', ...
  'uidx', [uidx_halo1(maps_halo1.x0_idx); uidx_connect(maps_connect.x0_idx);...
  uidx_connect(maps_connect.x1_idx); uidx_halo2(maps_halo2.x0_idx)]);
% add objetive
syms xx
T=zeros(ncheb);
for i=1:ncheb
    for j=1:ncheb
        T(i,j) = int(mychebyshevT(i-1,xx)*mychebyshevT(j-1,xx),xx,[-1,1]);
    end
end
odata.W = blkdiag(T,T,T);
obj_funcs = {@optcont_obj, @optcont_obj_du, @optcont_obj_dudu};
prob1 = coco_add_func(prob1, 'obj', obj_funcs{:}, odata, 'inactive', 'obj', 'uidx', uidx_connect(maps_connect.p_idx));


% ==============
% adjoint systems
prob1 = adjt_isol2coll(prob1, 'halo1');
prob1 = adjt_isol2coll(prob1, 'halo2');
[data_halo1, axidx_halo1] = coco_get_adjt_data(prob1, 'halo1.coll', 'data', 'axidx');
opt_halo1 = data_halo1.coll_opt;
[data_halo2, axidx_halo2] = coco_get_adjt_data(prob1, 'halo2.coll', 'data', 'axidx');
opt_halo2 = data_halo2.coll_opt;
% impose periodic BCs for two halo orbits
prob1 = coco_add_adjt(prob1, 'halo1_bc', 'aidx', [axidx_halo1(opt_halo1.x0_idx);axidx_halo1(opt_halo1.x1_idx)]);
prob1 = coco_add_adjt(prob1, 'halo2_bc', 'aidx', [axidx_halo2(opt_halo2.x0_idx);axidx_halo2(opt_halo2.x1_idx)]);
% Define inactive prameters
prob1 = coco_add_adjt(prob1, 'halo1_par_x0', 'd.x0_1', 'aidx', axidx_halo1(opt_halo1.x0_idx(2)));
prob1 = coco_add_adjt(prob1, 'halo1_par_T', 'd.T_1', 'aidx', axidx_halo1(opt_halo1.T_idx));
prob1 = coco_add_adjt(prob1, 'halo2_par_x0', 'd.x0_2', 'aidx', axidx_halo2(opt_halo2.x0_idx(2)));
prob1 = coco_add_adjt(prob1, 'halo2_par_T', 'd.T_2', 'aidx', axidx_halo2(opt_halo2.T_idx));

prob1 = adjt_isol2coll(prob1, 'connect');
[data_connect, axidx_connect] = coco_get_adjt_data(prob1, 'connect.coll', 'data', 'axidx');
opt_connect = data_connect.coll_opt;
% Time periods for each instances
prob1 = coco_add_adjt(prob1, 'T_bc', 'aidx', ...
  [axidx_connect(opt_connect.T0_idx);axidx_connect(opt_connect.T_idx)]);
% between instances
prob1 = coco_add_adjt(prob1, 'cont_bc', 'aidx',...
    [axidx_halo1(opt_halo1.x0_idx); axidx_connect(opt_connect.x0_idx);...
  axidx_connect(opt_connect.x1_idx); axidx_halo2(opt_halo2.x0_idx)]);
% objective
prob1 = coco_add_adjt(prob1, 'obj', 'd.obj', 'aidx', axidx_connect(opt_connect.p_idx));

% {'alpha' 'x0' 'T'} 
cont_args = {'obj', 'd.obj', 'd.halo1.coll.T0', 'd.halo2.coll.T0'...
     'd.T_1', 'd.T_2', 'd.x0_1', 'd.x0_2', 'p1', 'p11', 'p21', 'p2',...
     'p12', 'p22', 'p3', 'd.p4', 'd.p5', 'd.p6', 'd.p7',...
     'd.p8', 'd.p9', 'd.p10', 'd.p13', 'd.p14', 'd.p15', 'd.p16',...
     'd.p17', 'd.p18', 'd.p19', 'd.p20', 'd.p23', 'd.p24', 'd.p25',...
     'd.p26', 'd.p27', 'd.p28', 'd.p29', 'd.p30'};
coco(prob1, 'optcont1', [], 1, cont_args, [7.5,8.7]);

% check whether xa and xb are moved, tures out no movement
bd  = coco_bd_read('optcont1');
lab = 1:1:6;
for k=1:numel(lab)
sol1 = coll_read_solution('halo1','optcont1',lab(k));
sol2 = coll_read_solution('halo2','optcont1',lab(k));
sol3 = coll_read_solution('connect','optcont1',lab(k));
figure(1)
plot(sol1.tbp,sol1.xbp(:,1));hold on
figure(2)
plot(sol2.tbp,sol2.xbp(:,1));hold on
figure(3)
plot3(sol1.xbp(:,1),sol1.xbp(:,2),sol1.xbp(:,3));hold on
plot3(sol1.xbp(1,1),sol1.xbp(1,2),sol1.xbp(1,3),'ro');
figure(4)
plot3(sol2.xbp(:,1),sol2.xbp(:,2),sol2.xbp(:,3));hold on
plot3(sol2.xbp(1,1),sol2.xbp(1,2),sol2.xbp(1,3),'ro');
figure(5)
plot3(sol3.xbp(:,1),sol3.xbp(:,2),sol3.xbp(:,3));hold on
plot3(sol3.xbp(1,1),sol3.xbp(1,2),sol3.xbp(1,3),'ro');
end


%% switch branch to get non-trivial adjoint solutions
bd1   = coco_bd_read('optcont1');
BPlab = coco_bd_labs(bd1, 'BP');
lab   = min(BPlab);

prob = coco_prob();
prob = coco_set(prob, 'cont', 'NAdapt', 5, 'h_max', 5);
% zero problems
prob1 = ode_BP2coll(prob, 'halo1', 'optcont1', 'halo1', lab);
prob1 = ode_BP2coll(prob1, 'halo2', 'optcont1', 'halo2', lab);
prob1 = ode_BP2coll(prob1, 'connect', 'optcont1', 'connect', lab);

% extract data
[data_halo1, uidx_halo1] = coco_get_func_data(prob1, 'halo1.coll', 'data', 'uidx');
maps_halo1 = data_halo1.coll_seg.maps;
[data_halo2, uidx_halo2] = coco_get_func_data(prob1, 'halo2.coll', 'data', 'uidx');
maps_halo2 = data_halo2.coll_seg.maps;
[data_connect, uidx_connect] = coco_get_func_data(prob1, 'connect.coll', 'data', 'uidx');
maps_connect = data_connect.coll_seg.maps;
% impose periodic BCs for two halo orbits
prob1 = coco_add_func(prob1, 'halo1_bc', cr3bp_bc_funcs{:}, [], 'zero', 'uidx', [uidx_halo1(maps_halo1.x0_idx);uidx_halo1(maps_halo1.x1_idx)]);
prob1 = coco_add_func(prob1, 'halo2_bc', cr3bp_bc_funcs{:}, [], 'zero', 'uidx', [uidx_halo2(maps_halo2.x0_idx);uidx_halo2(maps_halo2.x1_idx)]);
% Define inactive prameters
prob1 = coco_add_pars(prob1, 'halo1_par_x0', uidx_halo1(maps_halo1.x0_idx(2)),'x0_1');
prob1 = coco_add_pars(prob1, 'halo1_par_T', uidx_halo1(maps_halo1.T_idx),'T_1');
prob1 = coco_add_pars(prob1, 'halo2_par_x0', uidx_halo2(maps_halo2.x0_idx(2)),'x0_2');
prob1 = coco_add_pars(prob1, 'halo2_par_T', uidx_halo2(maps_halo2.T_idx),'T_2');
% impose conditions for T0 and T of connected trajectories
prob1 = coco_add_func(prob1, 'T_bc', T_bc_funcs{:}, [], 'zero', ...
  'uidx', [uidx_connect(maps_connect.T0_idx);uidx_connect(maps_connect.T_idx)]);
% add links between these instances
[data_connect, uidx_connect] = coco_get_func_data(prob1, 'connect.coll', 'data', 'uidx');
maps_connect = data_connect.coll_seg.maps;
cont_bc_funcs = {@optcont_bc, @optcont_bc_du, @optcont_bc_dudu};
prob1 = coco_add_func(prob1, 'cont_bc', cont_bc_funcs{:}, [], 'zero', ...
  'uidx', [uidx_halo1(maps_halo1.x0_idx); uidx_connect(maps_connect.x0_idx);...
  uidx_connect(maps_connect.x1_idx); uidx_halo2(maps_halo2.x0_idx)]);
% add objetive
prob1 = coco_add_func(prob1, 'obj', obj_funcs{:}, odata, 'inactive', 'obj', 'uidx', uidx_connect(maps_connect.p_idx));

% branch switch data
chart = coco_read_solution('optcont1', lab, 'chart');
cdata = coco_get_chart_data(chart, 'lsol');


% adjoint systems
prob1 = adjt_BP2coll(prob1, 'halo1', 'optcont1', lab);
prob1 = adjt_BP2coll(prob1, 'halo2', 'optcont1', lab);
prob1 = adjt_BP2coll(prob1, 'connect', 'optcont1', lab);
% extract adjoint data
[data_halo1, axidx_halo1] = coco_get_adjt_data(prob1, 'halo1.coll', 'data', 'axidx');
opt_halo1 = data_halo1.coll_opt;
[data_halo2, axidx_halo2] = coco_get_adjt_data(prob1, 'halo2.coll', 'data', 'axidx');
opt_halo2 = data_halo2.coll_opt;
[data_connect, axidx_connect] = coco_get_adjt_data(prob1, 'connect.coll', 'data', 'axidx');
opt_connect = data_connect.coll_opt;

% impose periodic BCs for two halo orbits
[chart, aidx] = coco_read_adjoint('halo1_bc', 'optcont1', lab, 'chart', 'lidx');
prob1 = coco_add_adjt(prob1, 'halo1_bc', 'aidx', [axidx_halo1(opt_halo1.x0_idx);axidx_halo1(opt_halo1.x1_idx)],...
    'l0', chart.x, 'tl0', cdata.v(aidx));
[chart, aidx] = coco_read_adjoint('halo2_bc', 'optcont1', lab, 'chart', 'lidx');
prob1 = coco_add_adjt(prob1, 'halo2_bc', 'aidx', [axidx_halo2(opt_halo2.x0_idx);axidx_halo2(opt_halo2.x1_idx)],...
    'l0', chart.x, 'tl0', cdata.v(aidx));
% Define inactive prameters
[chart, aidx] = coco_read_adjoint('halo1_par_x0', 'optcont1', lab, 'chart', 'lidx');
prob1 = coco_add_adjt(prob1, 'halo1_par_x0', 'd.x0_1', 'aidx',...
    axidx_halo1(opt_halo1.x0_idx(2)), 'l0', chart.x, 'tl0', cdata.v(aidx));
[chart, aidx] = coco_read_adjoint('halo1_par_T', 'optcont1', lab, 'chart', 'lidx');
prob1 = coco_add_adjt(prob1, 'halo1_par_T', 'd.T_1', 'aidx',...
    axidx_halo1(opt_halo1.T_idx), 'l0', chart.x, 'tl0', cdata.v(aidx));
[chart, aidx] = coco_read_adjoint('halo2_par_x0', 'optcont1', lab, 'chart', 'lidx');
prob1 = coco_add_adjt(prob1, 'halo2_par_x0', 'd.x0_2', 'aidx',...
    axidx_halo2(opt_halo2.x0_idx(2)), 'l0', chart.x, 'tl0', cdata.v(aidx));
[chart, aidx] = coco_read_adjoint('halo2_par_T', 'optcont1', lab, 'chart', 'lidx');
prob1 = coco_add_adjt(prob1, 'halo2_par_T', 'd.T_2', 'aidx',...
    axidx_halo2(opt_halo2.T_idx), 'l0', chart.x, 'tl0', cdata.v(aidx));

% Time periods for connected instance
[chart, aidx] = coco_read_adjoint('T_bc', 'optcont1', lab, 'chart', 'lidx');
prob1 = coco_add_adjt(prob1, 'T_bc', 'aidx', ...
  [axidx_connect(opt_connect.T0_idx);axidx_connect(opt_connect.T_idx)],...
  'l0', chart.x, 'tl0', cdata.v(aidx));
% between instances
[chart, aidx] = coco_read_adjoint('cont_bc', 'optcont1', lab, 'chart', 'lidx');
prob1 = coco_add_adjt(prob1, 'cont_bc', 'aidx',...
    [axidx_halo1(opt_halo1.x0_idx); axidx_connect(opt_connect.x0_idx);...
  axidx_connect(opt_connect.x1_idx); axidx_halo2(opt_halo2.x0_idx)],...
  'l0', chart.x, 'tl0', cdata.v(aidx));
% objective
[chart, aidx] = coco_read_adjoint('obj', 'optcont1', lab, 'chart', 'lidx');
prob1 = coco_add_adjt(prob1, 'obj', 'd.obj', 'aidx', axidx_connect(opt_connect.p_idx),...
    'l0', chart.x, 'tl0', cdata.v(aidx));
% computational domain
dobj_int = [chart.x(1) 1.1];
prob1 = coco_add_event(prob1, 'opt', 'BP', 'd.obj', '>', 1);
cont_args = {'d.obj', 'obj', 'd.halo1.coll.T0', 'd.halo2.coll.T0'...
     'd.T_1', 'd.T_2', 'd.x0_1', 'd.x0_2', 'p1', 'p11', 'p21', 'p2',...
     'p12', 'p22', 'p3', 'd.p4', 'd.p5', 'd.p6', 'd.p7',...
     'd.p8', 'd.p9', 'd.p10', 'd.p13', 'd.p14', 'd.p15', 'd.p16',...
     'd.p17', 'd.p18', 'd.p19', 'd.p20', 'd.p23', 'd.p24', 'd.p25',...
     'd.p26', 'd.p27', 'd.p28', 'd.p29', 'd.p30'};

coco(prob1, 'optcont2', [], 1, cont_args, dobj_int); 


%% continue to let d.p=0
pid_ini = [1 11 21 2 12 22 3];
pid = zeros(1,23);
pid(1:3:22) = 3:10;
pid(2:3:23) = 13:20;
pid(3:3:24) = 23:30;
pid = pid(2:end);
nrun = numel(pid);

for k=1:nrun
    prunid = sprintf('optcont%d',k+1);
    crunid = sprintf('optcont%d',k+2);
    
    bd2 = coco_bd_read(prunid);
    lab = coco_bd_labs(bd2, 'opt');

    % zero problems
    prob1 = ode_coll2coll(prob, 'halo1', prunid, 'halo1', lab);
    prob1 = ode_coll2coll(prob1, 'halo2', prunid, 'halo2', lab);
    prob1 = ode_coll2coll(prob1, 'connect', prunid, 'connect', lab);   
    % extract data
    [data_halo1, uidx_halo1] = coco_get_func_data(prob1, 'halo1.coll', 'data', 'uidx');
    maps_halo1 = data_halo1.coll_seg.maps;
    [data_halo2, uidx_halo2] = coco_get_func_data(prob1, 'halo2.coll', 'data', 'uidx');
    maps_halo2 = data_halo2.coll_seg.maps;
    [data_connect, uidx_connect] = coco_get_func_data(prob1, 'connect.coll', 'data', 'uidx');
    maps_connect = data_connect.coll_seg.maps;
    % impose periodic BCs for two halo orbits
    prob1 = coco_add_func(prob1, 'halo1_bc', cr3bp_bc_funcs{:}, [], 'zero', 'uidx', [uidx_halo1(maps_halo1.x0_idx);uidx_halo1(maps_halo1.x1_idx)]);
    prob1 = coco_add_func(prob1, 'halo2_bc', cr3bp_bc_funcs{:}, [], 'zero', 'uidx', [uidx_halo2(maps_halo2.x0_idx);uidx_halo2(maps_halo2.x1_idx)]);
    % Define inactive prameters
    prob1 = coco_add_pars(prob1, 'halo1_par_x0', uidx_halo1(maps_halo1.x0_idx(2)),'x0_1');
    prob1 = coco_add_pars(prob1, 'halo1_par_T', uidx_halo1(maps_halo1.T_idx),'T_1');
    prob1 = coco_add_pars(prob1, 'halo2_par_x0', uidx_halo2(maps_halo2.x0_idx(2)),'x0_2');
    prob1 = coco_add_pars(prob1, 'halo2_par_T', uidx_halo2(maps_halo2.T_idx),'T_2');
    % impose conditions for T0 and T of connected trajectories
    prob1 = coco_add_func(prob1, 'T_bc', T_bc_funcs{:}, [], 'zero', ...
      'uidx', [uidx_connect(maps_connect.T0_idx);uidx_connect(maps_connect.T_idx)]);
    % add links between these instances
    [data_connect, uidx_connect] = coco_get_func_data(prob1, 'connect.coll', 'data', 'uidx');
    maps_connect = data_connect.coll_seg.maps;
    cont_bc_funcs = {@optcont_bc, @optcont_bc_du, @optcont_bc_dudu};
    prob1 = coco_add_func(prob1, 'cont_bc', cont_bc_funcs{:}, [], 'zero', ...
      'uidx', [uidx_halo1(maps_halo1.x0_idx); uidx_connect(maps_connect.x0_idx);...
      uidx_connect(maps_connect.x1_idx); uidx_halo2(maps_halo2.x0_idx)]);
    % add objetive
    prob1 = coco_add_func(prob1, 'obj', obj_funcs{:}, odata, 'inactive', 'obj', 'uidx', uidx_connect(maps_connect.p_idx));


    % adjoint systems
    prob1 = adjt_BP2coll(prob1, 'halo1', prunid, lab);
    prob1 = adjt_BP2coll(prob1, 'halo2', prunid, lab);
    prob1 = adjt_BP2coll(prob1, 'connect', prunid, lab);
    % extract adjoint data
    [data_halo1, axidx_halo1] = coco_get_adjt_data(prob1, 'halo1.coll', 'data', 'axidx');
    opt_halo1 = data_halo1.coll_opt;
    [data_halo2, axidx_halo2] = coco_get_adjt_data(prob1, 'halo2.coll', 'data', 'axidx');
    opt_halo2 = data_halo2.coll_opt;
    [data_connect, axidx_connect] = coco_get_adjt_data(prob1, 'connect.coll', 'data', 'axidx');
    opt_connect = data_connect.coll_opt;
    % impose periodic BCs for two halo orbits
    chart = coco_read_adjoint('halo1_bc', prunid, lab, 'chart');
    prob1 = coco_add_adjt(prob1, 'halo1_bc', 'aidx', [axidx_halo1(opt_halo1.x0_idx);axidx_halo1(opt_halo1.x1_idx)],...
        'l0', chart.x);
    chart = coco_read_adjoint('halo2_bc', prunid, lab, 'chart');
    prob1 = coco_add_adjt(prob1, 'halo2_bc', 'aidx', [axidx_halo2(opt_halo2.x0_idx);axidx_halo2(opt_halo2.x1_idx)],...
        'l0', chart.x);
    % Define inactive prameters
    chart = coco_read_adjoint('halo1_par_x0', prunid, lab, 'chart');
    prob1 = coco_add_adjt(prob1, 'halo1_par_x0', 'd.x0_1', 'aidx',...
        axidx_halo1(opt_halo1.x0_idx(2)), 'l0', chart.x);
    chart = coco_read_adjoint('halo1_par_T', prunid, lab, 'chart');
    prob1 = coco_add_adjt(prob1, 'halo1_par_T', 'd.T_1', 'aidx',...
        axidx_halo1(opt_halo1.T_idx), 'l0', chart.x);
    chart = coco_read_adjoint('halo2_par_x0', prunid, lab, 'chart');
    prob1 = coco_add_adjt(prob1, 'halo2_par_x0', 'd.x0_2', 'aidx',...
        axidx_halo2(opt_halo2.x0_idx(2)), 'l0', chart.x);
    chart = coco_read_adjoint('halo2_par_T', prunid, lab, 'chart');
    prob1 = coco_add_adjt(prob1, 'halo2_par_T', 'd.T_2', 'aidx',...
        axidx_halo2(opt_halo2.T_idx), 'l0', chart.x);
    % Time periods for connected instance
    chart = coco_read_adjoint('T_bc', prunid, lab, 'chart');
    prob1 = coco_add_adjt(prob1, 'T_bc', 'aidx', ...
      [axidx_connect(opt_connect.T0_idx);axidx_connect(opt_connect.T_idx)],...
      'l0', chart.x);
    % between instances
    chart = coco_read_adjoint('cont_bc', prunid, lab, 'chart');
    prob1 = coco_add_adjt(prob1, 'cont_bc', 'aidx',...
        [axidx_halo1(opt_halo1.x0_idx); axidx_connect(opt_connect.x0_idx);...
      axidx_connect(opt_connect.x1_idx); axidx_halo2(opt_halo2.x0_idx)],...
      'l0', chart.x);

    % objective
    chart = coco_read_adjoint('obj', prunid, lab, 'chart');
    prob1 = coco_add_adjt(prob1, 'obj', 'd.obj', 'aidx', axidx_connect(opt_connect.p_idx),...
        'l0', chart.x);

    % computational domain
    chart = coco_read_adjoint('connect.coll.pars', prunid, lab, 'chart');
    dp_pid = chart.x(pid(k));
    dpid  = sprintf('d.p%d',pid(k));
    if dp_pid>0
      dp_int = [-0.1 dp_pid];
      prob1 = coco_add_event(prob1, 'opt', 'BP', dpid, '<', 0);
    else
      dp_int = [dp_pid 0.1];
      prob1 = coco_add_event(prob1, 'opt', 'BP', dpid, '>', 0);
    end
    
    cont_args = cell(1,3*ncheb+8);
    cont_args{1} = dpid;
    cont_args{2} = 'obj';
    cont_args{3} = 'd.halo1.coll.T0';
    cont_args{4} = 'd.halo2.coll.T0';
    cont_args{5} = 'd.x0_1';
    cont_args{6} = 'd.x0_2';
    cont_args{7} = 'd.T_1';
    cont_args{8} = 'd.T_2';

    for kid=1:7
        cont_args{8+kid} = sprintf('p%d',pid_ini(kid));
    end
    for kid=1:k
        cont_args{15+kid} = sprintf('p%d',pid(kid));
    end
    for kid=k+1:nrun
        cont_args{15+kid} = sprintf('d.p%d',pid(kid));
    end  
    coco(prob1, crunid, [], cont_args, dp_int);
end

% the last run is cruid = optcont25, obj=7.0925

%% ---- we drive the movement of x0 -----%%
%% first release x0_1 but fix x0_2 --> free-fixed cases, cost is increased to 7.2487
k = k+1;
prunid = sprintf('optcont%d',k+1);
crunid = sprintf('optcont%d',k+2);

bd2 = coco_bd_read(prunid);
lab = coco_bd_labs(bd2, 'opt');
prob = coco_prob();
prob = coco_set(prob, 'cont', 'NAdapt', 10, 'h_max', 50);
prob = coco_set(prob, 'cont', 'ItMX', 150);
% zero problems
prob1 = ode_coll2coll(prob, 'halo1', prunid, 'halo1', lab);
prob1 = ode_coll2coll(prob1, 'halo2', prunid, 'halo2', lab);
prob1 = ode_coll2coll(prob1, 'connect', prunid, 'connect', lab);   
% extract data
[data_halo1, uidx_halo1] = coco_get_func_data(prob1, 'halo1.coll', 'data', 'uidx');
maps_halo1 = data_halo1.coll_seg.maps;
[data_halo2, uidx_halo2] = coco_get_func_data(prob1, 'halo2.coll', 'data', 'uidx');
maps_halo2 = data_halo2.coll_seg.maps;
[data_connect, uidx_connect] = coco_get_func_data(prob1, 'connect.coll', 'data', 'uidx');
maps_connect = data_connect.coll_seg.maps;
% impose periodic BCs for two halo orbits
prob1 = coco_add_func(prob1, 'halo1_bc', cr3bp_bc_funcs{:}, [], 'zero', 'uidx', [uidx_halo1(maps_halo1.x0_idx);uidx_halo1(maps_halo1.x1_idx)]);
prob1 = coco_add_func(prob1, 'halo2_bc', cr3bp_bc_funcs{:}, [], 'zero', 'uidx', [uidx_halo2(maps_halo2.x0_idx);uidx_halo2(maps_halo2.x1_idx)]);
% Define inactive prameters
prob1 = coco_add_pars(prob1, 'halo1_par_x0', uidx_halo1(maps_halo1.x0_idx(2)),'x0_1');
prob1 = coco_add_pars(prob1, 'halo1_par_T', uidx_halo1(maps_halo1.T_idx),'T_1');
prob1 = coco_add_pars(prob1, 'halo2_par_x0', uidx_halo2(maps_halo2.x0_idx(2)),'x0_2');
prob1 = coco_add_pars(prob1, 'halo2_par_T', uidx_halo2(maps_halo2.T_idx),'T_2');
% impose conditions for T0 and T of connected trajectories
prob1 = coco_add_func(prob1, 'T_bc', T_bc_funcs{:}, [], 'zero', ...
  'uidx', [uidx_connect(maps_connect.T0_idx);uidx_connect(maps_connect.T_idx)]);
% add links between these instances
[data_connect, uidx_connect] = coco_get_func_data(prob1, 'connect.coll', 'data', 'uidx');
maps_connect = data_connect.coll_seg.maps;
cont_bc_funcs = {@optcont_bc, @optcont_bc_du, @optcont_bc_dudu};
prob1 = coco_add_func(prob1, 'cont_bc', cont_bc_funcs{:}, [], 'zero', ...
  'uidx', [uidx_halo1(maps_halo1.x0_idx); uidx_connect(maps_connect.x0_idx);...
  uidx_connect(maps_connect.x1_idx); uidx_halo2(maps_halo2.x0_idx)]);
% add objetive
prob1 = coco_add_func(prob1, 'obj', obj_funcs{:}, odata, 'inactive', 'obj', 'uidx', uidx_connect(maps_connect.p_idx));


% adjoint systems
prob1 = adjt_BP2coll(prob1, 'halo1', prunid, lab);
prob1 = adjt_BP2coll(prob1, 'halo2', prunid, lab);
prob1 = adjt_BP2coll(prob1, 'connect', prunid, lab);
% extract adjoint data
[data_halo1, axidx_halo1] = coco_get_adjt_data(prob1, 'halo1.coll', 'data', 'axidx');
opt_halo1 = data_halo1.coll_opt;
[data_halo2, axidx_halo2] = coco_get_adjt_data(prob1, 'halo2.coll', 'data', 'axidx');
opt_halo2 = data_halo2.coll_opt;
[data_connect, axidx_connect] = coco_get_adjt_data(prob1, 'connect.coll', 'data', 'axidx');
opt_connect = data_connect.coll_opt;
% impose periodic BCs for two halo orbits
chart = coco_read_adjoint('halo1_bc', prunid, lab, 'chart');
prob1 = coco_add_adjt(prob1, 'halo1_bc', 'aidx', [axidx_halo1(opt_halo1.x0_idx);axidx_halo1(opt_halo1.x1_idx)],...
    'l0', chart.x);
chart = coco_read_adjoint('halo2_bc', prunid, lab, 'chart');
prob1 = coco_add_adjt(prob1, 'halo2_bc', 'aidx', [axidx_halo2(opt_halo2.x0_idx);axidx_halo2(opt_halo2.x1_idx)],...
    'l0', chart.x);
% Define inactive prameters
chart = coco_read_adjoint('halo1_par_x0', prunid, lab, 'chart');
prob1 = coco_add_adjt(prob1, 'halo1_par_x0', 'd.x0_1', 'aidx',...
    axidx_halo1(opt_halo1.x0_idx(2)), 'l0', chart.x);
chart = coco_read_adjoint('halo1_par_T', prunid, lab, 'chart');
prob1 = coco_add_adjt(prob1, 'halo1_par_T', 'd.T_1', 'aidx',...
    axidx_halo1(opt_halo1.T_idx), 'l0', chart.x);
chart = coco_read_adjoint('halo2_par_x0', prunid, lab, 'chart');
prob1 = coco_add_adjt(prob1, 'halo2_par_x0', 'd.x0_2', 'aidx',...
    axidx_halo2(opt_halo2.x0_idx(2)), 'l0', chart.x);
chart = coco_read_adjoint('halo2_par_T', prunid, lab, 'chart');
prob1 = coco_add_adjt(prob1, 'halo2_par_T', 'd.T_2', 'aidx',...
    axidx_halo2(opt_halo2.T_idx), 'l0', chart.x);
% Time periods for connected instance
chart = coco_read_adjoint('T_bc', prunid, lab, 'chart');
prob1 = coco_add_adjt(prob1, 'T_bc', 'aidx', ...
  [axidx_connect(opt_connect.T0_idx);axidx_connect(opt_connect.T_idx)],...
  'l0', chart.x);
% between instances
chart = coco_read_adjoint('cont_bc', prunid, lab, 'chart');
prob1 = coco_add_adjt(prob1, 'cont_bc', 'aidx',...
    [axidx_halo1(opt_halo1.x0_idx); axidx_connect(opt_connect.x0_idx);...
  axidx_connect(opt_connect.x1_idx); axidx_halo2(opt_halo2.x0_idx)],...
  'l0', chart.x);

% objective
chart = coco_read_adjoint('obj', prunid, lab, 'chart');
prob1 = coco_add_adjt(prob1, 'obj', 'd.obj', 'aidx', axidx_connect(opt_connect.p_idx),...
    'l0', chart.x);

% computational domain
chart = coco_read_adjoint('halo1_par_x0', prunid, lab, 'chart');
dp_pid = chart.x(1);
dpid  = 'd.x0_1';
if dp_pid>0
  dp_int = [-0.1 dp_pid];
  prob1 = coco_add_event(prob1, 'opt', 'BP', dpid, '<', 0);
else
  dp_int = [dp_pid 0.1];
  prob1 = coco_add_event(prob1, 'opt', 'BP', dpid, '>', 0);
end


cont_args = {'d.x0_1', 'obj', 'd.halo1.coll.T0', 'd.halo2.coll.T0'...
     'd.T_1', 'd.T_2', 'x0_1', 'd.x0_2', 'p1', 'p11', 'p21', 'p2',...
     'p12', 'p22', 'p3', 'p13', 'p23', 'p4', 'p14',...
     'p24', 'p5', 'p15', 'p25', 'p6', 'p16', 'p26',...
     'p7', 'p17', 'p27', 'p8', 'p18', 'p28', 'p9',...
     'p19', 'p29', 'p10', 'p20', 'p30'};
coco(prob1, crunid, [], cont_args, dp_int);


%% then release x0_2 to yield free-free case
k = k+1;
prunid = sprintf('optcont%d',k+1);
crunid = sprintf('optcont%d',k+2);

bd2 = coco_bd_read(prunid);
lab = coco_bd_labs(bd2, 'opt');
% prob = coco_set(prob, 'cont', 'NAdapt', 10, 'h_max', 50);
% prob = coco_set(prob, 'cont', 'ItMX', 150);
% zero problems
prob1 = ode_coll2coll(prob, 'halo1', prunid, 'halo1', lab);
prob1 = ode_coll2coll(prob1, 'halo2', prunid, 'halo2', lab);
prob1 = ode_coll2coll(prob1, 'connect', prunid, 'connect', lab);   
% extract data
[data_halo1, uidx_halo1] = coco_get_func_data(prob1, 'halo1.coll', 'data', 'uidx');
maps_halo1 = data_halo1.coll_seg.maps;
[data_halo2, uidx_halo2] = coco_get_func_data(prob1, 'halo2.coll', 'data', 'uidx');
maps_halo2 = data_halo2.coll_seg.maps;
[data_connect, uidx_connect] = coco_get_func_data(prob1, 'connect.coll', 'data', 'uidx');
maps_connect = data_connect.coll_seg.maps;
% impose periodic BCs for two halo orbits
prob1 = coco_add_func(prob1, 'halo1_bc', cr3bp_bc_funcs{:}, [], 'zero', 'uidx', [uidx_halo1(maps_halo1.x0_idx);uidx_halo1(maps_halo1.x1_idx)]);
prob1 = coco_add_func(prob1, 'halo2_bc', cr3bp_bc_funcs{:}, [], 'zero', 'uidx', [uidx_halo2(maps_halo2.x0_idx);uidx_halo2(maps_halo2.x1_idx)]);
% Define inactive prameters
prob1 = coco_add_pars(prob1, 'halo1_par_x0', uidx_halo1(maps_halo1.x0_idx(2)),'x0_1');
prob1 = coco_add_pars(prob1, 'halo1_par_T', uidx_halo1(maps_halo1.T_idx),'T_1');
prob1 = coco_add_pars(prob1, 'halo2_par_x0', uidx_halo2(maps_halo2.x0_idx(2)),'x0_2');
prob1 = coco_add_pars(prob1, 'halo2_par_T', uidx_halo2(maps_halo2.T_idx),'T_2');
% impose conditions for T0 and T of connected trajectories
prob1 = coco_add_func(prob1, 'T_bc', T_bc_funcs{:}, [], 'zero', ...
  'uidx', [uidx_connect(maps_connect.T0_idx);uidx_connect(maps_connect.T_idx)]);
% add links between these instances
[data_connect, uidx_connect] = coco_get_func_data(prob1, 'connect.coll', 'data', 'uidx');
maps_connect = data_connect.coll_seg.maps;
cont_bc_funcs = {@optcont_bc, @optcont_bc_du, @optcont_bc_dudu};
prob1 = coco_add_func(prob1, 'cont_bc', cont_bc_funcs{:}, [], 'zero', ...
  'uidx', [uidx_halo1(maps_halo1.x0_idx); uidx_connect(maps_connect.x0_idx);...
  uidx_connect(maps_connect.x1_idx); uidx_halo2(maps_halo2.x0_idx)]);
% add objetive
prob1 = coco_add_func(prob1, 'obj', obj_funcs{:}, odata, 'inactive', 'obj', 'uidx', uidx_connect(maps_connect.p_idx));

% adjoint systems
prob1 = adjt_BP2coll(prob1, 'halo1', prunid, lab);
prob1 = adjt_BP2coll(prob1, 'halo2', prunid, lab);
prob1 = adjt_BP2coll(prob1, 'connect', prunid, lab);
% extract adjoint data
[data_halo1, axidx_halo1] = coco_get_adjt_data(prob1, 'halo1.coll', 'data', 'axidx');
opt_halo1 = data_halo1.coll_opt;
[data_halo2, axidx_halo2] = coco_get_adjt_data(prob1, 'halo2.coll', 'data', 'axidx');
opt_halo2 = data_halo2.coll_opt;
[data_connect, axidx_connect] = coco_get_adjt_data(prob1, 'connect.coll', 'data', 'axidx');
opt_connect = data_connect.coll_opt;
% impose periodic BCs for two halo orbits
chart = coco_read_adjoint('halo1_bc', prunid, lab, 'chart');
prob1 = coco_add_adjt(prob1, 'halo1_bc', 'aidx', [axidx_halo1(opt_halo1.x0_idx);axidx_halo1(opt_halo1.x1_idx)],...
    'l0', chart.x);
chart = coco_read_adjoint('halo2_bc', prunid, lab, 'chart');
prob1 = coco_add_adjt(prob1, 'halo2_bc', 'aidx', [axidx_halo2(opt_halo2.x0_idx);axidx_halo2(opt_halo2.x1_idx)],...
    'l0', chart.x);
% Define inactive prameters
chart = coco_read_adjoint('halo1_par_x0', prunid, lab, 'chart');
prob1 = coco_add_adjt(prob1, 'halo1_par_x0', 'd.x0_1', 'aidx',...
    axidx_halo1(opt_halo1.x0_idx(2)), 'l0', chart.x);
chart = coco_read_adjoint('halo1_par_T', prunid, lab, 'chart');
prob1 = coco_add_adjt(prob1, 'halo1_par_T', 'd.T_1', 'aidx',...
    axidx_halo1(opt_halo1.T_idx), 'l0', chart.x);
chart = coco_read_adjoint('halo2_par_x0', prunid, lab, 'chart');
prob1 = coco_add_adjt(prob1, 'halo2_par_x0', 'd.x0_2', 'aidx',...
    axidx_halo2(opt_halo2.x0_idx(2)), 'l0', chart.x);
chart = coco_read_adjoint('halo2_par_T', prunid, lab, 'chart');
prob1 = coco_add_adjt(prob1, 'halo2_par_T', 'd.T_2', 'aidx',...
    axidx_halo2(opt_halo2.T_idx), 'l0', chart.x);
% Time periods for connected instance
chart = coco_read_adjoint('T_bc', prunid, lab, 'chart');
prob1 = coco_add_adjt(prob1, 'T_bc', 'aidx', ...
  [axidx_connect(opt_connect.T0_idx);axidx_connect(opt_connect.T_idx)],...
  'l0', chart.x);
% between instances
chart = coco_read_adjoint('cont_bc', prunid, lab, 'chart');
prob1 = coco_add_adjt(prob1, 'cont_bc', 'aidx',...
    [axidx_halo1(opt_halo1.x0_idx); axidx_connect(opt_connect.x0_idx);...
  axidx_connect(opt_connect.x1_idx); axidx_halo2(opt_halo2.x0_idx)],...
  'l0', chart.x);

% objective
chart = coco_read_adjoint('obj', prunid, lab, 'chart');
prob1 = coco_add_adjt(prob1, 'obj', 'd.obj', 'aidx', axidx_connect(opt_connect.p_idx),...
    'l0', chart.x);

% computational domain
chart = coco_read_adjoint('halo2_par_x0', prunid, lab, 'chart');
dp_pid = chart.x(1);
dpid  = 'd.x0_2';
if dp_pid>0
  dp_int = [-0.1 dp_pid];
  prob1 = coco_add_event(prob1, 'opt', 'BP', dpid, '<', 0);
else
  dp_int = [dp_pid 0.1];
  prob1 = coco_add_event(prob1, 'opt', 'BP', dpid, '>', 0);
end


cont_args = {'d.x0_2', 'obj', 'd.halo1.coll.T0', 'd.halo2.coll.T0'...
     'd.T_1', 'd.T_2', 'x0_1', 'x0_2', 'p1', 'p11', 'p21', 'p2',...
     'p12', 'p22', 'p3', 'p13', 'p23', 'p4', 'p14',...
     'p24', 'p5', 'p15', 'p25', 'p6', 'p16', 'p26',...
     'p7', 'p17', 'p27', 'p8', 'p18', 'p28', 'p9',...
     'p19', 'p29', 'p10', 'p20', 'p30'};
coco(prob1, crunid, [], cont_args, dp_int);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% first release x0_2 but fix x0_1 --> fixed-free cases, cost is increased to 7.1428
k = 24;
prunid = sprintf('optcont%d',k+1);
crunid = sprintf('optcont%d',k+4); % save to optcont28

bd2 = coco_bd_read(prunid);
lab = coco_bd_labs(bd2, 'opt');

% zero problems
prob1 = ode_coll2coll(prob, 'halo1', prunid, 'halo1', lab);
prob1 = ode_coll2coll(prob1, 'halo2', prunid, 'halo2', lab);
prob1 = ode_coll2coll(prob1, 'connect', prunid, 'connect', lab);   
% extract data
[data_halo1, uidx_halo1] = coco_get_func_data(prob1, 'halo1.coll', 'data', 'uidx');
maps_halo1 = data_halo1.coll_seg.maps;
[data_halo2, uidx_halo2] = coco_get_func_data(prob1, 'halo2.coll', 'data', 'uidx');
maps_halo2 = data_halo2.coll_seg.maps;
[data_connect, uidx_connect] = coco_get_func_data(prob1, 'connect.coll', 'data', 'uidx');
maps_connect = data_connect.coll_seg.maps;
% impose periodic BCs for two halo orbits
prob1 = coco_add_func(prob1, 'halo1_bc', cr3bp_bc_funcs{:}, [], 'zero', 'uidx', [uidx_halo1(maps_halo1.x0_idx);uidx_halo1(maps_halo1.x1_idx)]);
prob1 = coco_add_func(prob1, 'halo2_bc', cr3bp_bc_funcs{:}, [], 'zero', 'uidx', [uidx_halo2(maps_halo2.x0_idx);uidx_halo2(maps_halo2.x1_idx)]);
% Define inactive prameters
prob1 = coco_add_pars(prob1, 'halo1_par_x0', uidx_halo1(maps_halo1.x0_idx(2)),'x0_1');
prob1 = coco_add_pars(prob1, 'halo1_par_T', uidx_halo1(maps_halo1.T_idx),'T_1');
prob1 = coco_add_pars(prob1, 'halo2_par_x0', uidx_halo2(maps_halo2.x0_idx(2)),'x0_2');
prob1 = coco_add_pars(prob1, 'halo2_par_T', uidx_halo2(maps_halo2.T_idx),'T_2');
% impose conditions for T0 and T of connected trajectories
prob1 = coco_add_func(prob1, 'T_bc', T_bc_funcs{:}, [], 'zero', ...
  'uidx', [uidx_connect(maps_connect.T0_idx);uidx_connect(maps_connect.T_idx)]);
% add links between these instances
[data_connect, uidx_connect] = coco_get_func_data(prob1, 'connect.coll', 'data', 'uidx');
maps_connect = data_connect.coll_seg.maps;
cont_bc_funcs = {@optcont_bc, @optcont_bc_du, @optcont_bc_dudu};
prob1 = coco_add_func(prob1, 'cont_bc', cont_bc_funcs{:}, [], 'zero', ...
  'uidx', [uidx_halo1(maps_halo1.x0_idx); uidx_connect(maps_connect.x0_idx);...
  uidx_connect(maps_connect.x1_idx); uidx_halo2(maps_halo2.x0_idx)]);
% add objetive
prob1 = coco_add_func(prob1, 'obj', obj_funcs{:}, odata, 'inactive', 'obj', 'uidx', uidx_connect(maps_connect.p_idx));


% adjoint systems
prob1 = adjt_BP2coll(prob1, 'halo1', prunid, lab);
prob1 = adjt_BP2coll(prob1, 'halo2', prunid, lab);
prob1 = adjt_BP2coll(prob1, 'connect', prunid, lab);
% extract adjoint data
[data_halo1, axidx_halo1] = coco_get_adjt_data(prob1, 'halo1.coll', 'data', 'axidx');
opt_halo1 = data_halo1.coll_opt;
[data_halo2, axidx_halo2] = coco_get_adjt_data(prob1, 'halo2.coll', 'data', 'axidx');
opt_halo2 = data_halo2.coll_opt;
[data_connect, axidx_connect] = coco_get_adjt_data(prob1, 'connect.coll', 'data', 'axidx');
opt_connect = data_connect.coll_opt;
% impose periodic BCs for two halo orbits
chart = coco_read_adjoint('halo1_bc', prunid, lab, 'chart');
prob1 = coco_add_adjt(prob1, 'halo1_bc', 'aidx', [axidx_halo1(opt_halo1.x0_idx);axidx_halo1(opt_halo1.x1_idx)],...
    'l0', chart.x);
chart = coco_read_adjoint('halo2_bc', prunid, lab, 'chart');
prob1 = coco_add_adjt(prob1, 'halo2_bc', 'aidx', [axidx_halo2(opt_halo2.x0_idx);axidx_halo2(opt_halo2.x1_idx)],...
    'l0', chart.x);
% Define inactive prameters
chart = coco_read_adjoint('halo1_par_x0', prunid, lab, 'chart');
prob1 = coco_add_adjt(prob1, 'halo1_par_x0', 'd.x0_1', 'aidx',...
    axidx_halo1(opt_halo1.x0_idx(2)), 'l0', chart.x);
chart = coco_read_adjoint('halo1_par_T', prunid, lab, 'chart');
prob1 = coco_add_adjt(prob1, 'halo1_par_T', 'd.T_1', 'aidx',...
    axidx_halo1(opt_halo1.T_idx), 'l0', chart.x);
chart = coco_read_adjoint('halo2_par_x0', prunid, lab, 'chart');
prob1 = coco_add_adjt(prob1, 'halo2_par_x0', 'd.x0_2', 'aidx',...
    axidx_halo2(opt_halo2.x0_idx(2)), 'l0', chart.x);
chart = coco_read_adjoint('halo2_par_T', prunid, lab, 'chart');
prob1 = coco_add_adjt(prob1, 'halo2_par_T', 'd.T_2', 'aidx',...
    axidx_halo2(opt_halo2.T_idx), 'l0', chart.x);
% Time periods for connected instance
chart = coco_read_adjoint('T_bc', prunid, lab, 'chart');
prob1 = coco_add_adjt(prob1, 'T_bc', 'aidx', ...
  [axidx_connect(opt_connect.T0_idx);axidx_connect(opt_connect.T_idx)],...
  'l0', chart.x);
% between instances
chart = coco_read_adjoint('cont_bc', prunid, lab, 'chart');
prob1 = coco_add_adjt(prob1, 'cont_bc', 'aidx',...
    [axidx_halo1(opt_halo1.x0_idx); axidx_connect(opt_connect.x0_idx);...
  axidx_connect(opt_connect.x1_idx); axidx_halo2(opt_halo2.x0_idx)],...
  'l0', chart.x);

% objective
chart = coco_read_adjoint('obj', prunid, lab, 'chart');
prob1 = coco_add_adjt(prob1, 'obj', 'd.obj', 'aidx', axidx_connect(opt_connect.p_idx),...
    'l0', chart.x);

% computational domain
chart = coco_read_adjoint('halo2_par_x0', prunid, lab, 'chart');
dp_pid = chart.x(1);
dpid  = 'd.x0_2';
if dp_pid>0
  dp_int = [-0.1 dp_pid];
  prob1 = coco_add_event(prob1, 'opt', 'BP', dpid, '<', 0);
else
  dp_int = [dp_pid 0.1];
  prob1 = coco_add_event(prob1, 'opt', 'BP', dpid, '>', 0);
end

cont_args = {'d.x0_2', 'obj', 'd.halo1.coll.T0', 'd.halo2.coll.T0'...
     'd.T_1', 'd.T_2', 'd.x0_1', 'x0_2', 'p1', 'p11', 'p21', 'p2',...
     'p12', 'p22', 'p3', 'p13', 'p23', 'p4', 'p14',...
     'p24', 'p5', 'p15', 'p25', 'p6', 'p16', 'p26',...
     'p7', 'p17', 'p27', 'p8', 'p18', 'p28', 'p9',...
     'p19', 'p29', 'p10', 'p20', 'p30'};
coco(prob1, crunid, [], cont_args, dp_int);



%% then we further release x0_1 --> free-free cases, cost is increased to 7.2978
k = 27;
prunid = sprintf('optcont%d',k+1);
crunid = sprintf('optcont%d',k+2);

bd2 = coco_bd_read(prunid);
lab = coco_bd_labs(bd2, 'opt');
% zero problems
prob1 = ode_coll2coll(prob, 'halo1', prunid, 'halo1', lab);
prob1 = ode_coll2coll(prob1, 'halo2', prunid, 'halo2', lab);
prob1 = ode_coll2coll(prob1, 'connect', prunid, 'connect', lab);   
% extract data
[data_halo1, uidx_halo1] = coco_get_func_data(prob1, 'halo1.coll', 'data', 'uidx');
maps_halo1 = data_halo1.coll_seg.maps;
[data_halo2, uidx_halo2] = coco_get_func_data(prob1, 'halo2.coll', 'data', 'uidx');
maps_halo2 = data_halo2.coll_seg.maps;
[data_connect, uidx_connect] = coco_get_func_data(prob1, 'connect.coll', 'data', 'uidx');
maps_connect = data_connect.coll_seg.maps;
% impose periodic BCs for two halo orbits
prob1 = coco_add_func(prob1, 'halo1_bc', cr3bp_bc_funcs{:}, [], 'zero', 'uidx', [uidx_halo1(maps_halo1.x0_idx);uidx_halo1(maps_halo1.x1_idx)]);
prob1 = coco_add_func(prob1, 'halo2_bc', cr3bp_bc_funcs{:}, [], 'zero', 'uidx', [uidx_halo2(maps_halo2.x0_idx);uidx_halo2(maps_halo2.x1_idx)]);
% Define inactive prameters
prob1 = coco_add_pars(prob1, 'halo1_par_x0', uidx_halo1(maps_halo1.x0_idx(2)),'x0_1');
prob1 = coco_add_pars(prob1, 'halo1_par_T', uidx_halo1(maps_halo1.T_idx),'T_1');
prob1 = coco_add_pars(prob1, 'halo2_par_x0', uidx_halo2(maps_halo2.x0_idx(2)),'x0_2');
prob1 = coco_add_pars(prob1, 'halo2_par_T', uidx_halo2(maps_halo2.T_idx),'T_2');
% impose conditions for T0 and T of connected trajectories
prob1 = coco_add_func(prob1, 'T_bc', T_bc_funcs{:}, [], 'zero', ...
  'uidx', [uidx_connect(maps_connect.T0_idx);uidx_connect(maps_connect.T_idx)]);
% add links between these instances
[data_connect, uidx_connect] = coco_get_func_data(prob1, 'connect.coll', 'data', 'uidx');
maps_connect = data_connect.coll_seg.maps;
cont_bc_funcs = {@optcont_bc, @optcont_bc_du, @optcont_bc_dudu};
prob1 = coco_add_func(prob1, 'cont_bc', cont_bc_funcs{:}, [], 'zero', ...
  'uidx', [uidx_halo1(maps_halo1.x0_idx); uidx_connect(maps_connect.x0_idx);...
  uidx_connect(maps_connect.x1_idx); uidx_halo2(maps_halo2.x0_idx)]);
% add objetive
prob1 = coco_add_func(prob1, 'obj', obj_funcs{:}, odata, 'inactive', 'obj', 'uidx', uidx_connect(maps_connect.p_idx));


% adjoint systems
prob1 = adjt_BP2coll(prob1, 'halo1', prunid, lab);
prob1 = adjt_BP2coll(prob1, 'halo2', prunid, lab);
prob1 = adjt_BP2coll(prob1, 'connect', prunid, lab);
% extract adjoint data
[data_halo1, axidx_halo1] = coco_get_adjt_data(prob1, 'halo1.coll', 'data', 'axidx');
opt_halo1 = data_halo1.coll_opt;
[data_halo2, axidx_halo2] = coco_get_adjt_data(prob1, 'halo2.coll', 'data', 'axidx');
opt_halo2 = data_halo2.coll_opt;
[data_connect, axidx_connect] = coco_get_adjt_data(prob1, 'connect.coll', 'data', 'axidx');
opt_connect = data_connect.coll_opt;
% impose periodic BCs for two halo orbits
chart = coco_read_adjoint('halo1_bc', prunid, lab, 'chart');
prob1 = coco_add_adjt(prob1, 'halo1_bc', 'aidx', [axidx_halo1(opt_halo1.x0_idx);axidx_halo1(opt_halo1.x1_idx)],...
    'l0', chart.x);
chart = coco_read_adjoint('halo2_bc', prunid, lab, 'chart');
prob1 = coco_add_adjt(prob1, 'halo2_bc', 'aidx', [axidx_halo2(opt_halo2.x0_idx);axidx_halo2(opt_halo2.x1_idx)],...
    'l0', chart.x);
% Define inactive prameters
chart = coco_read_adjoint('halo1_par_x0', prunid, lab, 'chart');
prob1 = coco_add_adjt(prob1, 'halo1_par_x0', 'd.x0_1', 'aidx',...
    axidx_halo1(opt_halo1.x0_idx(2)), 'l0', chart.x);
chart = coco_read_adjoint('halo1_par_T', prunid, lab, 'chart');
prob1 = coco_add_adjt(prob1, 'halo1_par_T', 'd.T_1', 'aidx',...
    axidx_halo1(opt_halo1.T_idx), 'l0', chart.x);
chart = coco_read_adjoint('halo2_par_x0', prunid, lab, 'chart');
prob1 = coco_add_adjt(prob1, 'halo2_par_x0', 'd.x0_2', 'aidx',...
    axidx_halo2(opt_halo2.x0_idx(2)), 'l0', chart.x);
chart = coco_read_adjoint('halo2_par_T', prunid, lab, 'chart');
prob1 = coco_add_adjt(prob1, 'halo2_par_T', 'd.T_2', 'aidx',...
    axidx_halo2(opt_halo2.T_idx), 'l0', chart.x);
% Time periods for connected instance
chart = coco_read_adjoint('T_bc', prunid, lab, 'chart');
prob1 = coco_add_adjt(prob1, 'T_bc', 'aidx', ...
  [axidx_connect(opt_connect.T0_idx);axidx_connect(opt_connect.T_idx)],...
  'l0', chart.x);
% between instances
chart = coco_read_adjoint('cont_bc', prunid, lab, 'chart');
prob1 = coco_add_adjt(prob1, 'cont_bc', 'aidx',...
    [axidx_halo1(opt_halo1.x0_idx); axidx_connect(opt_connect.x0_idx);...
  axidx_connect(opt_connect.x1_idx); axidx_halo2(opt_halo2.x0_idx)],...
  'l0', chart.x);

% objective
chart = coco_read_adjoint('obj', prunid, lab, 'chart');
prob1 = coco_add_adjt(prob1, 'obj', 'd.obj', 'aidx', axidx_connect(opt_connect.p_idx),...
    'l0', chart.x);

% computational domain
chart = coco_read_adjoint('halo1_par_x0', prunid, lab, 'chart');
dp_pid = chart.x(1);
dpid  = 'd.x0_1';
if dp_pid>0
  dp_int = [-0.1 dp_pid];
  prob1 = coco_add_event(prob1, 'opt', 'BP', dpid, '<', 0);
else
  dp_int = [dp_pid 0.1];
  prob1 = coco_add_event(prob1, 'opt', 'BP', dpid, '>', 0);
end

cont_args = {'d.x0_1', 'obj', 'd.halo1.coll.T0', 'd.halo2.coll.T0'...
     'd.T_1', 'd.T_2', 'x0_1', 'x0_2', 'p1', 'p11', 'p21', 'p2',...
     'p12', 'p22', 'p3', 'p13', 'p23', 'p4', 'p14',...
     'p24', 'p5', 'p15', 'p25', 'p6', 'p16', 'p26',...
     'p7', 'p17', 'p27', 'p8', 'p18', 'p28', 'p9',...
     'p19', 'p29', 'p10', 'p20', 'p30'};
coco(prob1, crunid, [], cont_args, dp_int);




%% plot optimal results
crunid = 'optcont25'
bd6 = coco_bd_read(crunid);
% id1 = [1 11 21 2 12 22 3];
% id2 = setdiff((1:30),id1)
% control input
p = coco_bd_col(bd6, par_args);
idx = coco_bd_idxs(bd6, 'opt');
popt = p(:,idx);
save orbit_10_terms_fixed_fixed.mat popt
tf = 2;
tp = (0:0.01:tf)';   % physical domain [0,tf]
t  = (tp-tf/2)/tf*2; % standard domain [-1,1] 




% control input
lab = coco_bd_labs(bd6, 'opt');
sol = coll_read_solution('connect', crunid, lab);
save connected_orbits_fixed_fixed.mat sol
x  = sol.xbp;
t  = sol.tbp;
tp = (t+1)*tf/2;
u = zeros(3,numel(t));
for i=1:3
    for j=1:ncheb
        u(i,:) = u(i,:) + p((i-1)*ncheb+j,idx)*mychebyshevT(j-1,t)';
    end
end
figure(6)
plot(tp,u(1,:),'r-'); hold on
plot(tp,u(2,:),'b-');
plot(tp,u(3,:),'k-');
xlabel('t');
ylabel('control input');
legend('u1','u2','u3');
figure(7)
plot(tp,x(:,1),'r-.'); hold on
plot(tp,x(:,2),'b-.');
plot(tp,x(:,3),'k-.');

figure(5)
plot3(x(:,1),x(:,2),x(:,3),'b.')

toc