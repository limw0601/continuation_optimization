% minimize (x-2)^2+2*(y-1)^2 s.t. G1:=x+4*y-3<=0, G2:=-x+y<=0
% Starting point (x0,y0)=(1,-2) such that G1(x0,y0)<0 & G2(x0,y0)<0

%% First run to detect a fold point
prob = coco_prob();
prob = coco_set(prob, 'cont', 'NPR', 2);
prob = coco_add_func(prob, 'ineq1', @G1_func, @G1_dfunc, [],...
    'inequality', 'G1', 'u0', [1 -2]);       % G1 inequality
prob = coco_add_func(prob, 'ineq2', @G2_func, @G2_dfunc, [],...
    'inequality', 'G2', 'uidx', [1 2]);      % G2 inequality
prob = coco_add_func(prob, 'obj', @obj_func, @obj_dfunc, [],...
    'inactive', 'obj', 'uidx', [1,2]);       % Objective
prob = coco_add_pars(prob, 'atfic_mon', 2, 'y', 'inactive'); % Psi_2

% adjoint
prob = coco_add_adjt(prob, 'ineq1', 'ncp.G1');
prob = coco_add_adjt(prob, 'ineq2', 'ncp.G2', 'aidx', [1,2]);
prob = coco_add_adjt(prob, 'obj', 'd.obj', 'aidx', [1,2]);
prob = coco_add_adjt(prob, 'atfic_mon', 'd.y', 'aidx', 2);
coco(prob, 'run1', [], 1, {'obj' 'd.obj' 'd.y' 'ncp.G1' 'ncp.G2' 'G1' 'G2' 'y'},[0 19]);


%% Switch to secondary branch and continue along it until d.obj=1
bd  = coco_bd_read('run1');
lab = coco_bd_labs(bd, 'BP');
prob = coco_prob();
prob = coco_set(prob, 'cont', 'NPR', 2);
chart = coco_read_solution('ineq1', 'run1', lab,'chart');
cdata = coco_get_chart_data(chart, 'lsol');
prob  = coco_add_func(prob, 'ineq1', @G1_func, @G1_dfunc, [],...
    'inequality', 'G1', 'u0', chart.x);
prob  = coco_add_func(prob, 'ineq2', @G2_func, @G2_dfunc, [],...
    'inequality', 'G2', 'uidx', [1 2]);
prob  = coco_add_func(prob, 'obj', @obj_func, @obj_dfunc, [],...
    'inactive', 'obj', 'uidx', [1 2]);
prob  = coco_add_pars(prob, 'atfic_mon', 2, 'y', 'inactive');

% adjoint
[chart, lidx] = coco_read_adjoint('ineq1', 'run1', lab, 'chart', 'lidx');
prob  = coco_add_adjt(prob, 'ineq1', 'ncp.G1', 'l0', chart.x, 'tl0', cdata.v(lidx));
[chart, lidx] = coco_read_adjoint('ineq2', 'run1', lab, 'chart', 'lidx');
prob  = coco_add_adjt(prob, 'ineq2', 'ncp.G2', 'aidx', [1,2], 'l0', chart.x, 'tl0', cdata.v(lidx));
[chart, lidx] = coco_read_adjoint('obj', 'run1', lab, 'chart', 'lidx');
prob = coco_add_adjt(prob, 'obj', 'd.obj', 'aidx', [1,2], 'l0', chart.x, 'tl0', cdata.v(lidx));
[chart, lidx] = coco_read_adjoint('atfic_mon', 'run1', lab, 'chart', 'lidx');
prob = coco_add_adjt(prob, 'atfic_mon', 'd.y', 'aidx', 2, 'l0', chart.x, 'tl0', cdata.v(lidx));
coco(prob, 'run2', [], 1, {'d.obj' 'obj' 'd.y' 'ncp.G1' 'ncp.G2' 'G1' 'G2'},[0 1]);

%% Last continuation run until d.y=0
bd  = coco_bd_read('run2');
lab = coco_bd_labs(bd, 'EP');
lab = max(lab);
prob  = coco_prob();
prob  = coco_set(prob, 'cont', 'NPR', 3);
chart = coco_read_solution('ineq1', 'run2',lab,'chart');
prob  = coco_add_func(prob,'ineq1', @G1_func, @G1_dfunc, [],...
    'inequality', 'G1', 'u0', chart.x);
prob  = coco_add_func(prob,'ineq2', @G2_func, @G2_dfunc, [],...
    'inequality', 'G2', 'uidx', [1 2]);
prob  = coco_add_func(prob,'obj', @obj_func, @obj_dfunc, [],...
    'inactive', 'obj', 'uidx', [1 2]);
prob  = coco_add_pars(prob,'atfic_mon', 2, 'y', 'inactive');

% adjoint
chart = coco_read_adjoint('ineq1', 'run2', lab, 'chart');
prob  = coco_add_adjt(prob, 'ineq1', 'ncp.G1', 'l0', chart.x);
chart = coco_read_adjoint('ineq2', 'run2', lab, 'chart');
prob  = coco_add_adjt(prob, 'ineq2', 'ncp.G2', 'aidx', [1,2], 'l0', chart.x);
chart = coco_read_adjoint('obj', 'run2', lab, 'chart');
prob  = coco_add_adjt(prob, 'obj', 'd.obj', 'aidx', [1,2], 'l0', chart.x);
chart = coco_read_adjoint('atfic_mon', 'run2', lab, 'chart');
prob  = coco_add_adjt(prob, 'atfic_mon', 'd.y', 'aidx', 2, 'l0', chart.x);
coco(prob, 'run3', [], 1, {'d.y' 'obj' 'y' 'ncp.G1' 'ncp.G2' 'd.obj' 'G1' 'G2'},[0 12]);

% print optimal point
bd  = coco_bd_read('run3');
lab = coco_bd_labs(bd, 'EP');
lab = max(lab);
chart = coco_read_solution('obj','run3',lab,'chart');
disp('Optimal point:');
chart.x

%% Plot projection of continuation paths
x     = [];
sigma = [];
for k=1:15
    chart = coco_read_solution('obj','run1',k,'chart');
    x  = [x chart.x(1:2)];
    chart = coco_read_adjoint('ineq1','run1',k,'chart');
    s1 = chart.x;
    chart = coco_read_adjoint('atfic_mon','run1',k,'chart');
    s2 = chart.x;
    chart = coco_read_adjoint('obj','run1',k,'chart');
    s3 = chart.x;
    s  = [s1;s2;s3];
    sigma = [sigma s];
end
x1 = x(1,:);
x2 = x(2,:);
mu = (x1-2).^2+2*(x2-1).^2;
sigma1 = sigma(1,:);
eta_y  = sigma(2,:);
sigma3 = sigma(3,:);
figure(1);
plot3(x1,x2,mu,'k-','LineWidth',2); hold on
bd  = coco_bd_read('run1');
idx1 = coco_bd_labs(bd, 'FP');
plot3(x1(idx1), x2(idx1), mu(idx1), 'ro', 'MarkerFaceColor', 'r')
figure(2)
plot3(sigma1, eta_y, sigma3,'k-','LineWidth',2); hold on
plot3(sigma1(idx1), eta_y, sigma3(idx1), 'ro', 'MarkerFaceColor', 'r')


x     = [];
sigma = [];
for k=1:19
    chart = coco_read_solution('obj','run2',k,'chart');
    x  = [x chart.x(1:2)];
    chart = coco_read_adjoint('ineq1','run2',k,'chart');
    s1 = chart.x;
    chart = coco_read_adjoint('atfic_mon','run2',k,'chart');
    s2 = chart.x;
    chart = coco_read_adjoint('obj','run2',k,'chart');
    s3 = chart.x;
    s  = [s1;s2;s3];
    sigma = [sigma s];
end
x1 = x(1,:);
x2 = x(2,:);
mu = (x1-2).^2+2*(x2-1).^2;
sigma1 = sigma(1,:);
eta_y  = sigma(2,:);
sigma3 = sigma(3,:);
figure(1);
plot3(x1,x2,mu,'k-','LineWidth',2); hold on
bd  = coco_bd_read('run2');
idx1 = coco_bd_labs(bd, 'EP');
idx1 = max(idx1);
figure(2)
plot3(sigma1, eta_y, sigma3,'k-','LineWidth',2); hold on
plot3(sigma1(idx1), eta_y(idx1), sigma3(idx1), 'bo', 'MarkerFaceColor', 'b')


x     = [];
sigma = [];
for k=1:23
    chart = coco_read_solution('obj','run3',k,'chart');
    x  = [x chart.x(1:2)];
    chart = coco_read_adjoint('ineq1','run3',k,'chart');
    s1 = chart.x;
    chart = coco_read_adjoint('atfic_mon','run3',k,'chart');
    s2 = chart.x;
    chart = coco_read_adjoint('obj','run3',k,'chart');
    s3 = chart.x;
    s  = [s1;s2;s3];
    sigma = [sigma s];
end
x1 = x(1,:);
x2 = x(2,:);
mu = (x1-2).^2+2*(x2-1).^2;
sigma1 = sigma(1,:);
eta_y  = sigma(2,:);
sigma3 = sigma(3,:);
figure(1);
plot3(x1,x2,mu,'k-','LineWidth',2); hold on
bd  = coco_bd_read('run3');
idx1 = coco_bd_labs(bd, 'EP');
idx1 = max(idx1);
plot3(x1(idx1), x2(idx1), mu(idx1), 'ko', 'MarkerFaceColor', 'k');
plot3(x1(18), x2(18), mu(18), 'ms', 'MarkerFaceColor', 'm');
set(gca,'LineWidth',1.2);
set(gca,'FontSize',14);
xlabel('$x$','interpreter','latex','FontSize',14);
ylabel('$y$','interpreter','latex','FontSize',14);
zlabel('$\mu_1$','interpreter','latex','FontSize',14);
grid on 
box on
[y, z] = meshgrid(-2:0.1:0.5,-5:0.5:25);  
Xv = @(y) 3 - 4*y;
surf(Xv(y),y,z,'FaceColor', 0.9*[1 1 1], 'FaceAlpha', 0.5, 'LineStyle', 'none');
Xv = @(y) y;
surf(Xv(y),y,z,'FaceColor', 0.9*[1 1 1], 'FaceAlpha', 0.5, 'LineStyle', 'none');
axis([0.0 3.5 -2.5 0.5 -5 25])

figure(2)
plot3(sigma1, eta_y, sigma3,'k-','LineWidth',2); hold on
plot3(sigma1(idx1), eta_y(idx1), sigma3(idx1), 'ko', 'MarkerFaceColor', 'k');
plot3(sigma1(18), eta_y(18), sigma3(18), 'ms', 'MarkerFaceColor', 'm');
set(gca,'LineWidth',1.2);
set(gca,'FontSize',14);
xlabel('$\sigma_1$','interpreter','latex','FontSize',14);
ylabel('$\eta_2$','interpreter','latex','FontSize',14);
zlabel('$\eta_1$','interpreter','latex','FontSize',14);
grid on 
box on
axis([-0.2 1 0 15 0 1.2])
