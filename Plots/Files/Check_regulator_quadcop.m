%% QUADROTOR BALANCING PENDULUM MODEL PREDICTIVE CONTROL SIMULATION
%
% MATLAB simulation of the paper A Flying Inverted Pendulum by Markus Hehn 
% and Raffaello D'Andrea using a Model Predictive Controller

%% INIT
clc
clear
addpath('functions/');
cvx_clear




%% DEFINE CONSTANTS
load('Benchmark_model.mat')
% cvx_solver sedumi
cvx_precision best

v = 2;

[bic_sys, matrices] = Bicycle_model(v);
check_controllability(bic_sys);


%% DISCRETIZE SYSTEM

h = 0.1;

sysd = c2d(bic_sys,h,'tustin');

% LTI
LTI.A = sysd.A;
LTI.B = sysd.B;
LTI.C = sysd.C;

A = sysd.A;
B = sysd.B;
C = sysd.C;

% Dimensions
dim.nx = size(A,1);
dim.nu = size(B,2);
dim.N = 25;

%% MODEL PREDICTIVE CONTROL
steps = 60;
T = steps*h;
% initial state
phi_0 = 10 * pi/180;                  % Initial roll angle in radians
delta_0 = 0 * pi/180;                 % Initial steering angle in radians
phi_dot_0 = 10 * pi/180;                  % Initial roll angle speed in radians/s
delta_dot_0 = 0 * pi/180;                 % Initial steering angle speed in radians/s
x0 = [phi_0; delta_0; phi_dot_0 ; delta_dot_0];          % Initial state
LTI.x0 = x0;
x(:,1) = x0;

% desired reference (x,y,z,yaw)
r = [zeros(1,steps);     % x reference
     zeros(1,steps)];    % yaw reference

% B_ref relates reference to states x_ref = B_ref*r

B_ref = zeros(4,2);

x = zeros(length(A(:,1)),T);    % state trajectory
u = zeros(length(B(1,:)),T);    % control inputs
y = zeros(length(C(:,1)),T);    % measurements 
t = zeros(1,T);                 % time vector

Vf = zeros(1,T);                % terminal cost sequence
l = zeros(1,T);                 % stage cost sequence

x(:,1) = x0';

% Define MPC Control Problem

% MPC cost function
%          N-1
% V(u_N) = Sum 1/2[ x(k)'Qx(k) + u(k)'Ru(k) ] + x(N)'Sx(N) 
%          k = 0

% tuning weights
Q = 100000*eye(size(A,1));            % state cost
R = eye(length(B(1,:)));    % input cost

% terminal cost = unconstrained optimal cost (Lec 5 pg 6)
[S,K_LQR,~] = idare(A,B,Q,R);       % terminal cost % OLD: S = 10*eye(size(A));

% S = Q;
% S = 1000*eye(16);

% prediction horizon
N = 25; 

Qbar = kron(Q,eye(N));
Rbar = kron(R,eye(N));
Sbar = S;

LTI.A = A;
LTI.B = B;
LTI.C = C;

dim.N = N;
dim.nx = size(A,1);
dim.nu = size(B,2);
dim.ny = size(C,1);

[P,Z,W] = predmodgen_quad(LTI,dim);
              
H = (Z'*Qbar*Z + Rbar + 2*W'*Sbar*W);
d = (x0'*P'*Qbar*Z + 2*x0'*(A^N)'*Sbar*W)';
 
%%

% State constraints 
% F_constr * statevector <= e_state
F_constr = kron(eye(dim.nx),[1; -1]);    
phi_max = pi/3;
e_state = phi_max*[0.5; 0.5; 0.5; 0.5; 0.42; 0.42; 0.84; 0.84];

x_const_set = constraintmaker(dim.N,F_constr,e_state);

% input consraints
% E_constr * inputvector <= e_input
E_constr = kron(eye(dim.nu),[1; -1]);
max_lean_moment = (abs(z_rf)*sin(phi_max))*g*(0.2)*m_rf;
max_steer_torque = 5;
e_input = [max_lean_moment; max_lean_moment; max_steer_torque; max_steer_torque];

u_const_set = constraintmaker(dim.N,E_constr,e_input);

constraint.F = F_constr;
constraint.es = e_state;
constraint.E = E_constr;
constraint.ei = e_input;


for k = 1:1:steps
    t(k) = (k-1)*h;
%     if ( mod(t(k),1) == 0 )
%         fprintf('t = %d sec \n', t(k));
%     end
    
    % determine reference states based on reference input r
    x_ref = B_ref*r(:,k);
    x0 = x(:,k) - x_ref;
    d = (x0'*P'*Qbar*Z + 2*x0'*(A^N)'*Sbar*W)';
    
    % compute control action
    cvx_begin quiet
        variable u_N(2*N)
        minimize ( (1/2)*quad_form(u_N,H) + d'*u_N )
        % input constraints
        % input constraints
        
        u_const_set.E*u_N <= u_const_set.e
    
        % state constraints
        (P*x0 + Z*u_N) <= repmat(phi_max*[0.5; 0.5; 0.42; 0.84],N,1); 
        (P*x0 + Z*u_N) >= -repmat(phi_max*[0.5; 0.5; 0.42; 0.84],N,1);
    cvx_end
    
    u(:,k) = u_N(1:2); % MPC control action
    
    % apply control action
    x(:,k+1) = A*x(:,k) + B*u(:,k); % + B_ref*r(:,k);
    y(:,k) = C*x(:,k);
    
    % stability analysis
%     Q = 10*eye(16);
%     R = 0.1*eye(4);
    
    [X,eigvals,K] = dare(A,B,Q,R);
    Vf(k) = 0.5*x(:,k)'*X*x(:,k);
    l(k) = 0.5*x(:,k)'*Q*x(:,k) + 0.1*u(:,k)'*R*u(:,k);
end

% states_trajectory: Nx16 matrix of trajectory of 16 states
states_trajectory = y';

stairs(t,x(1,1:end-1))

%% PLOT RESULTS
% 
% % Stability plots
% Vf_diff = Vf(2:end)-Vf(1:end-1);
% 
% Vfkp1 = Vf(2:end);
% Vfk = Vf(1:end-1);
% 
% lQ = l(2:end);
% % lQ = l(1:end-1);
% 
% kt = t(1:end-1);

% figure(123);
% clf;
% hold on;
% % stairs(kt,Vfkp1);
% % stairs(kt,Vfk);
% % stairs(kt,lQ);
% stairs(kt,Vf_diff);
% stairs(kt,-0.3*lQ);
% % stairs(kt,Vfk-lQ);
% % stairs(kt,Vfkp1-(Vfk-lQ));
% grid on;

% legend('$V_f(f(x,u)) - V_f(x)$','$-l(x,u)$','interpreter','latex');

% legend('Vf','Vf(k+1)-Vf(k)','stage l(k)','Vf - l');

% % show 3D simulation
% X = states_trajectory(:,[3 9 13 11 5 15 1 7]);
% visualize_quadrotor_trajectory(states_trajectory(:,[3 9 13 11 5 15 1 7]));
% 
% saved_data.t = t;
% saved_data.x = states_trajectory;
% saved_data.u = u;

%% Basic Plots
% plot 2D results fo state trajectories
% plot_2D_plots(t, states_trajectory);
% % 
% % % plot the inputs
% plot_inputs(t,u,u_limit);

%% Comparison Plots

% plot_comparison_S_different(); % 543
% plot_comparison_R_different(); % 544
% plot_comparison_Q_different(); % 546
% 
% plot_comparison_R_inputs();    % 589
% 
% %%
% plot_comparison_MPC_LQR();     % 567
% 
% plot_comparison_horizon();      % 789
