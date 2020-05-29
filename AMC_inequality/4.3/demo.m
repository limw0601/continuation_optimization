% minimize (x-2)^2+2*(y-1)^2 s.t. G1:=x+4*y-3<=0, G2:=-x+y<=0
% Starting point (x0,y0)=(1,2) such that G1(x0,y0)>0 & G2(x0,y0)>0

%% First continuation run
prob = coco_prob();
prob = coco_set(prob, 'cont', 'NPR', 1);
prob = coco_add_func(prob, 'ineq1', @G1_func, @G1_dfunc, [],...
    'inequality', 'G1', 'u0', [1 2]);       % G1 inequality
prob = coco_add_func(prob, 'ineq2', @G2_func, @G2_dfunc, [],...
    'inequality', 'G2', 'uidx', [1 2]);     % G2 inequality
prob = coco_add_func(prob, 'obj', @obj_func, @obj_dfunc, [],...
    'inactive', 'obj', 'uidx', [1 2]);      % Objective

% adjoint
prob = coco_add_adjt(prob, 'ineq1', 'ncp.G1');
prob = coco_add_adjt(prob, 'ineq2', 'ncp.G2', 'aidx', [1,2]);
prob = coco_add_adjt(prob, 'obj', 'd.obj', 'aidx', [1,2]);
coco(prob, 'run1', [], 1, {'d.obj' 'obj' 'ncp.G1' 'ncp.G2' 'G1' 'G2'},[0 1]);

%% Second run to drive ncp.G1 (kappa1) to zero
bd  = coco_bd_read('run1');
lab = coco_bd_labs(bd, 'EP');
lab = max(lab);
prob  = coco_prob();
prob  = coco_set(prob, 'cont', 'NPR', 2);
chart = coco_read_solution('ineq1', 'run1', lab, 'chart');
prob  = coco_add_func(prob,'ineq1', @G1_func, @G1_dfunc, [],...
    'inequality', 'G1', 'u0', chart.x);
prob  = coco_add_func(prob,'ineq2', @G2_func, @G2_dfunc, [],...
    'inequality', 'G2', 'uidx', [1 2]);
prob  = coco_add_func(prob,'obj', @obj_func, @obj_dfunc, [],...
    'inactive', 'obj', 'uidx', [1 2]);

% adjoint
chart = coco_read_adjoint('ineq1', 'run1', lab, 'chart');
prob  = coco_add_adjt(prob, 'ineq1', 'ncp.G1', 'l0', chart.x);
chart = coco_read_adjoint('ineq2', 'run1', lab, 'chart');
prob  = coco_add_adjt(prob, 'ineq2', 'ncp.G2', 'aidx', [1,2], 'l0', chart.x);
chart = coco_read_adjoint('obj', 'run1', lab, 'chart');
prob  = coco_add_adjt(prob, 'obj', 'd.obj', 'aidx', [1,2], 'l0', chart.x);
coco(prob, 'run2', [], 1, {'ncp.G1' 'obj' 'd.obj' 'ncp.G2' 'G1' 'G2'},[0 12]);

%% Last run to drive ncp.G2 (kappa2) to zero
bd  = coco_bd_read('run2');
lab = coco_bd_labs(bd, 'EP');
lab = max(lab);
prob  = coco_prob();
prob  = coco_set(prob, 'cont', 'NPR', 1);
chart = coco_read_solution('ineq1','run2',lab,'chart');
prob  = coco_add_func(prob, 'ineq1', @G1_func, @G1_dfunc, [],...
    'inequality', 'G1', 'u0', chart.x);
prob  = coco_add_func(prob, 'ineq2', @G2_func, @G2_dfunc, [],...
    'inequality', 'G2', 'uidx', [1,2]);
prob  = coco_add_func(prob, 'obj', @obj_func, @obj_dfunc, [],...
    'inactive', 'obj', 'uidx', [1,2]);

% adjoint
chart = coco_read_adjoint('ineq1', 'run2', lab, 'chart');
prob  = coco_add_adjt(prob, 'ineq1', 'ncp.G1', 'l0', chart.x);
chart = coco_read_adjoint('ineq2', 'run2', lab, 'chart');
prob  = coco_add_adjt(prob, 'ineq2', 'ncp.G2', 'aidx', [1,2], 'l0', chart.x);
chart = coco_read_adjoint('obj', 'run2', lab, 'chart');
prob  = coco_add_adjt(prob, 'obj', 'd.obj', 'aidx', [1,2], 'l0', chart.x);
coco(prob, 'run3', [], 1, {'ncp.G2' 'obj' 'd.obj' 'ncp.G1' 'G1' 'G2'},[0 2]);

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
for k=1:12
    chart = coco_read_solution('obj','run1',k,'chart');
    x  = [x chart.x(1:2)];
    chart = coco_read_adjoint('ineq1','run1',k,'chart');
    s1 = chart.x;
    chart = coco_read_adjoint('ineq2','run1',k,'chart');
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
sigma2 = sigma(2,:);
sigma3 = sigma(3,:);
figure(1);
plot3(x1,x2,mu,'k-','LineWidth',2); hold on
bd  = coco_bd_read('run1');
idx1 = coco_bd_labs(bd, 'EP');
idx1 = max(idx1);
h1 = x1+4*x2-3;
kappa1 = (sigma1.^2+h1.^2).^0.5-sigma1+h1;
h2 = -x1+x2;
kappa2 = (sigma2.^2+h2.^2).^0.5-sigma2+h2;
plot3(x1(idx1), x2(idx1), mu(idx1), 'yo', 'MarkerFaceColor', 'y')
figure(2)
plot(sigma1,sigma2,'k-','LineWidth',2); hold on
plot(sigma1(idx1), sigma2(idx1), 'yo', 'MarkerFaceColor', 'y')


x     = [];
sigma = [];
for k=1:16
    chart = coco_read_solution('obj','run2',k,'chart');
    x  = [x chart.x(1:2)];
    chart = coco_read_adjoint('ineq1','run2',k,'chart');
    s1 = chart.x;
    chart = coco_read_adjoint('ineq2','run2',k,'chart');
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
sigma2 = sigma(2,:);
sigma3 = sigma(3,:);
figure(1);
plot3(x1,x2,mu,'k-','LineWidth',2); hold on
bd  = coco_bd_read('run2');
idx1 = coco_bd_labs(bd, 'EP');
idx1 = max(idx1);
h1 = x1+4*x2-3;
kappa1 = (sigma1.^2+h1.^2).^0.5-sigma1+h1;
h2 = -x1+x2;
kappa2 = (sigma2.^2+h2.^2).^0.5-sigma2+h2;
plot3(x1(idx1), x2(idx1), mu(idx1), 'bo', 'MarkerFaceColor', 'b')
figure(2)
plot(sigma1,sigma2,'k-','LineWidth',2); hold on
plot(sigma1(idx1), sigma2(idx1), 'bo', 'MarkerFaceColor', 'b')


x     = [];
sigma = [];
for k=1:9
    chart = coco_read_solution('obj','run3',k,'chart');
    x  = [x chart.x(1:2)];
    chart = coco_read_adjoint('ineq1','run3',k,'chart');
    s1 = chart.x;
    chart = coco_read_adjoint('ineq2','run3',k,'chart');
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
sigma2 = sigma(2,:);
sigma3 = sigma(3,:);
figure(1);
plot3(x1,x2,mu,'k-','LineWidth',2); hold on
bd  = coco_bd_read('run3');
idx1 = coco_bd_labs(bd, 'EP');
idx1 = max(idx1);
h1 = x1+4*x2-3;
kappa1 = (sigma1.^2+h1.^2).^0.5-sigma1+h1;
h2 = -x1+x2;
kappa2 = (sigma2.^2+h2.^2).^0.5-sigma2+h2;
plot3(x1(idx1), x2(idx1), mu(idx1), 'ko', 'MarkerFaceColor', 'k')
set(gca,'LineWidth',1.2);
set(gca,'FontSize',14);
xlabel('$x$','interpreter','latex','FontSize',14);
ylabel('$y$','interpreter','latex','FontSize',14);
zlabel('$\mu_1$','interpreter','latex','FontSize',14);
grid on 
box on
[y, z] = meshgrid(0:0.1:2,0:0.1:3);  
Xv = @(y) 3 - 4*y;
surf(Xv(y),y,z,'FaceColor', 0.8*[1 1 1], 'FaceAlpha', 0.5, 'LineStyle', 'none');
Xv = @(y) y;
surf(Xv(y),y,z,'FaceColor', 0.9*[1 1 1], 'FaceAlpha', 0.5, 'LineStyle', 'none');
axis([1 2 0 2 0 3])

figure(2)
plot(sigma1,sigma2,'k-','LineWidth',2); hold on
plot(sigma1(idx1), sigma2(idx1), 'ko', 'MarkerFaceColor', 'k')
set(gca,'LineWidth',1.2);
set(gca,'FontSize',14);
xlabel('$\sigma_1$','interpreter','latex','FontSize',14);
ylabel('$\sigma_2$','interpreter','latex','FontSize',14);
grid on 
box on