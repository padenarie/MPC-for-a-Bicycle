clc
clear all
close all
addpath('functions/');

cvx_solver sedumi
cvx_clear
cvx_precision best

%% Initialize system benchmark model
load('Benchmark_model.mat')

% Parameters to be tuned
h = 0.1;                    % Sampling time [s]
v = 2;                      % Bicycle velocity [m/s]
N = 16;                      % Control horizon [-]
steps = 400;                 

% initialize states
phi_0 = 0 * pi/180;        % Initial roll angle in radians
delta_0 = 0 * pi/180;       % Initial steering angle in radians
phi_dot_0 = 0 * pi/180;    % Initial roll angle speed in radians/s
delta_dot_0 = 0 * pi/180;   % Initial steering angle speed in radians/s

% Weighting matrices
% Weight matrices
Q = [10000 0 0 0            % Weigthing matrix for the states
      0 10000 0 0
      0 0 1 0
      0 0 0 1];
 
R = 1*eye(2);          % Weighting matrix for the inputs


%% Discretize system
[bic_sys, matrices] = Bicycle_model(v);
check_controllability(bic_sys);

% Disturbance dynamics
rho_air = 1.293;              % Air density [kg/m^3]
m_bike = m_rf+m_ff+m_fw+m_rw; % Mass of bicycle and rider [kg]
A_bike = 0.5;                 % Contact surface of wind on bicycle and rider [m^2]
v_max = 4;                    % Maximum velocity of the wind [m/s]

% Sin for defining the wind disturbance
ysin(1) = 0;
for k = 1:steps/2
    ysin(k+1) = sin(2*pi/(steps/3)*k);
    if abs(ysin(k+1)) < 10^(-8)
        ysin(k+1) = 0;
    end
    phi_dotdot(1,k) = (1/(m_bike*abs(z_rf)))*0.5*rho_air*A_bike*(v_max*ysin(1,k))^2;
end


%% Discretize system
T = steps*h;
B_d_cont = [0  0
            0  1
            1  0
            0  0];
B_d_aug = [bic_sys.B B_d_cont];
bic_sys_d_aug = ss(bic_sys.A, B_d_aug, bic_sys.C,[]);

sysd = c2d(bic_sys_d_aug,h,'ZOH');

% LTI
LTI.A = sysd.A;
LTI.B = sysd.B(:,1:2);
LTI.C = sysd.C;

B_d = sysd.B(:,3:4);

A = sysd.A;
B = sysd.B(:,1:2);
C = sysd.C;

% Dimensions
dim.nx = size(A,1);
dim.ny = size(C,1);
dim.nu = size(B,2);
dim.nd = 2;
dim.N = 16;

% Augmented Dimensions
dimaug.nx = dim.nx + dim.nd;
dimaug.ny = dim.ny;
dimaug.nu = dim.nu;
dimaug.nd = dim.nd;
dimaug.N = dim.N;

%% MPC regulator
% Initialize matrices for speed
x = zeros(dim.nx,steps+1);
u = zeros(dim.nu,steps);
u_store_dif = zeros(dim.nu*dim.N,steps);
y = zeros(dim.ny,steps);

t_time = zeros(1,steps+1);
x_aug_hat = zeros(dimaug.nx,steps+1);
y_aug_hat = zeros(dimaug.ny,steps);

u_ref = zeros(dim.nu,steps);
x_ref = zeros(dim.nx,steps+1);
y_ref = zeros(dim.ny,steps);
d_hat = zeros(dimaug.nd,steps+1);

% Initial states
x0 = [phi_0; delta_0; phi_dot_0 ; delta_dot_0];  
LTI.x0 = x0;
x(:,1) = x0;

weight.Q = Q;
weight.R = R;

% Terminal cost
[P,K_LQR,~] = idare(A,B,Q,R);    % Solving the DARE to find the optimal terminal weighting matrix
weight.P = P;
check_posdef(P)

% Prediction matrices
predmod = predmodgen(LTI,dim);

% Cost matrices
[H,F_x0less] = optimatrix(predmod,dim,weight);

% State constraints 
F_constr = kron(eye(dim.nx),[1; -1]);    
phi_max = pi/3;
e_state = phi_max*[0.5; 0.5; 0.5; 0.5; 0.42; 0.42; 0.84; 0.84];

x_const_set = constraintmaker(dim.N,F_constr,e_state);

% input consraints
E_constr = kron(eye(dim.nu),[1; -1]);
max_lean_moment = 40/29.2*(m_rf+m_ff+m_fw+m_rw);
max_steer_torque = 5;
e_input = [max_lean_moment; max_lean_moment; max_steer_torque; max_steer_torque];

u_const_set = constraintmaker(dim.N,E_constr,e_input);

constraint.F = F_constr;
constraint.es = e_state;
constraint.E = E_constr;
constraint.ei = e_input;

% OTS constraints
A_ots = blkdiag(E_constr,F_constr);
b_ots = [e_state;e_input];

%% Disturbance
dist_phi = [phi_dotdot zeros(1,steps/2)];
road_noise = 0.001*0.84*phi_max*randn(1,steps);
dist_steer = road_noise;
d = [dist_phi; dist_steer];

%% Disturbance matrices
C_d = [1  0
       0  1];

I_dis = eye(dimaug.nd,dimaug.nu);
  
aug.A = [A                       B_d;
         zeros(dimaug.nd,dim.nx) I_dis];
     
aug.B = [B; zeros(dimaug.nd,dimaug.nu)];

aug.C = [C C_d];

A_aug = aug.A;
B_aug = aug.B;
C_aug = aug.C;

%% Observer
% check rank of augmented system
I_A = eye(dim.nx);
test = [(I_A-A) -B_d 
            C    C_d];
rnk_aug = rank(test);

% Weight tuning
weight_aug.Q = 1*eye(dimaug.nx); weight_aug.Q(5,5) = 10000; weight_aug.Q(6,6) = 10; 
weight_aug.R = 1*eye(dimaug.nu);

Q_aug_obs = weight_aug.Q;
R_aug_obs = weight_aug.R;

[~,L,~,info] = idare(A_aug',C_aug',Q_aug_obs,R_aug_obs);
L = L';

test = A_aug - L*C_aug;
labda = eig(test);          % All real parts of eigenvalues are smaller than 1

%% OTS matrices
Aeq_OTS = [(eye(dim.nx) - A) -B
                 C   zeros(dim.nu)];
       
Q_OTS = 1*eye(dim.nx);
R_OTS = 1*eye(dim.nu);

H_OTS = blkdiag(Q_OTS,R_OTS);

%% Xf invariant set calc
A_LQR = A - B*K_LQR;
C_xf = [K_LQR; eye(dim.nx)];
s = size(A_ots,1);
[Xf_set_H, Xf_set_h, kstar] = calcInvariantXf(A_LQR,C_xf,A_ots,b_ots,s,dim);

% check if the x-location is within X_f
inSet = all(Xf_set_H*x0 <= Xf_set_h);

%% MPC
for i = 1:steps
    
    % Online OTS
    beq_OTS = [B_d*d_hat(:,i) 
               y_ref(:,i) - C_d*d_hat(:,i)];
    
    opts = optimoptions('quadprog','Display','off');
    [xu,~,exitflag] = quadprog(H_OTS,zeros(dim.nx+dim.nu,1),[],[],Aeq_OTS,beq_OTS,[],[],[],opts);
    
%     xu_store(:,i) = xu;
    x_ref(:,i) = xu(1:dim.nx);
    u_ref(:,i) = xu(dim.nx+1:end);
      
    t_time(i+1) = t_time(i) + h;
    x0_est_dif = x_aug_hat(1:dim.nx,i) - x_ref(:,i);
%     x0_est_dif_store(:,i) = x0_est_dif;
    F = x0_est_dif'*F_x0less;
    
    cvx_begin quiet
        variable u_pred_dif(dim.nu*dim.N)
        minimize ( (1/2)*u_pred_dif'*H*u_pred_dif + F*u_pred_dif )
    
        % input constraints
        u_const_set.E*(u_pred_dif + repmat(u_ref(:,i),dim.N,1))<= u_const_set.e
    
        % state constraints
        x_const_set.E*(predmod.T(dim.nx+1:end,:)*x0_est_dif + predmod.S(dim.nx+1:end,:)*(u_pred_dif + repmat(u_ref(:,i),dim.N,1))) <= x_const_set.e
        
        % Terminal constraint
        Xf_set_H*(predmod.T(end-dim.nx+1:end,:)*x0_est_dif + predmod.S(end-dim.nx+1:end,:)*(u_pred_dif + repmat(u_ref(:,i),dim.N,1))) <= Xf_set_h;
        
    cvx_end
    
%     u_store_pre(:,i) = u_pred_dif;

    u_store_dif(:,i) = u_pred_dif;
    u(:,i) = u_pred_dif(1:2,1) + u_ref(:,i);
    
    x(:,i+1) = A*x(:,i) + B*u(:,i) + B_d*d(:,i);
    y(:,i) = C*x(:,i) + C_d*d(:,i);
    
    y_aug_hat(:,i) = C_aug*x_aug_hat(:,i);
    x_aug_hat(:,i+1) = A_aug*x_aug_hat(:,i) + B_aug*u(:,i) + L*(y(:,i) - y_aug_hat(:,i));
    
    d_hat(:,i+1) = x_aug_hat(dim.nx+1:end,i+1);
    
    % errors
    ed(:,i) = d(:,i)-d_hat(:,i);
    e(:,i) = y(:,i) - y_aug_hat(:,i);
    
    Vf(i) = 0.5*x(:,i)'*P*x(:,i);
    l(i) = 0.5*x(:,i)'*Q*x(:,i);
end

% Reshape data
% d_store = [d;d_hat(1,1:end-1)];
% x_store = [x;x_aug_hat(1:dim.nx,1:end)];
% x_aug_hat_deg = x_aug_hat(1:dim.nx,:)*180/pi;
% x_deg = x*180/pi;
% d_hat_deg = d_hat*180/pi;
% d_deg = d*180/pi;
% y_deg = y*180/pi;

% Stability plots
Vf_diff = Vf(2:end)-Vf(1:end-1);

Vfkp1 = Vf(2:end);
Vfk = Vf(1:end-1);

lQ = l(2:end);

kt = t_time(1:end-2);


%% Plot
t_time = linspace(0,T,steps+1);

figure(1)
hold on
yline(0.5*phi_max*180/pi);
yline(-0.5*phi_max*180/pi);
stairs(t_time(1,1:end-1),y_deg(1,:))
stairs(t_time(1,1:end-1),y_deg(2,:))
legend 'boundup' 'boundlow' 'roll angle' 'steering angle'
hold off

figure(2)
hold on
stairs(t_time(1,1:end-1),u(1,:))
stairs(t_time(1,1:end-1),u(2,:))
legend 'roll torque' 'steering torque' 
hold off

figure(3)
hold on
stairs(t_time(1,1:end),x_aug_hat_deg(1,:))
stairs(t_time(1,1:end),x_deg(1,:))
stairs(t_time(1,1:end),x_aug_hat_deg(2,:))
stairs(t_time(1,1:end),x_deg(2,:))
yline(0.5*phi_max*180/pi);
yline(-0.5*phi_max*180/pi);

legend('roll angle est','roll angle','steering angle est', 'steering angle')
hold off

figure(4)
hold on
stairs(t_time(1,1:end-1),d_hat_deg(1,1:end-1))
stairs(t_time(1,1:end-1),d_deg(1,:))
stairs(t_time(1,1:end-1),d_hat_deg(2,1:end-1))
stairs(t_time(1,1:end-1),d_deg(2,:))
legend 'est dist phi' 'dist phi' 'est dist steer' 'dist steer'
hold off

figure(5)
subplot(3,1,1)
hold on
grid on
stairs(t_time(1,1:end-1),x(1,1:end-1),'LineWidth',2)
stairs(t_time(1,1:end-1),x_aug_hat(1,1:end-1),'--','LineWidth',1)
stairs(t_time(1,1:end-1),x(2,1:end-1),'LineWidth',2)
stairs(t_time(1,1:end-1),x_aug_hat(2,1:end-1),'--','LineWidth',1)
ylabel('$\phi$ and $\delta$ [rad]','interpreter','latex')
legend('$\phi$','$\phi_{est}$','$\delta$','$\delta_{est}$','interpreter','latex')

subplot(3,1,2)
hold on
grid on
stairs(t_time(1,1:end-1),x(3,1:end-1),'LineWidth',2)
stairs(t_time(1,1:end-1),x_aug_hat(3,1:end-1),'--','LineWidth',1)
stairs(t_time(1,1:end-1),x(4,1:end-1),'LineWidth',2)
stairs(t_time(1,1:end-1),x_aug_hat(4,1:end-1),'--','LineWidth',1)
ylabel('$\dot\phi$ and $\dot\delta$ [rad/s]','interpreter','latex')
legend('$\dot{\phi}$','$\dot{\phi_{est}}$', '$\dot{\delta}$','$\dot{\delta_{est}}$','interpreter','latex') 
hold off

subplot(3,1,3)
hold on
grid on
stairs(t_time(1,1:end-1),d(1,1:end),'LineWidth',2)
stairs(t_time(1,1:end-1),d_hat(1,1:end-1),'--','LineWidth',1)
xlabel('Time [s]')
ylabel('$\dot\phi$ and $\dot\delta$ [rad/s]','interpreter','latex')
legend( '$d$', '$d_{est}$','interpreter','latex')
hold off

% figure(6)
% hold on
% stairs(t_time(1,1:end-1),ed(1,1:end))
% stairs(t_time(1,1:end-1),ed(2,1:end))
% stairs(t_time(1,1:end-1),e(1,1:end))
% stairs(t_time(1,1:end-1),e(2,1:end))
% legend 'dist error phi' 'dist error steer' 'error est output phi' 'error est output steer'
% hold off

figure(7)
hold on;
stairs(kt,-lQ);
stairs(kt,Vf_diff);
% stairs(kt,-0.3*lQ);
% stairs(kt,Vfk-lQ);
% stairs(kt,Vfkp1-(Vfk-lQ));
grid on;

% figure(8)
% hold on
% stairs(t_time(1,1:end-1),u_store_ref(1,:))
% stairs(t_time(1,1:end-1),u_store_ref(2,:))
% hold off
