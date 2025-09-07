clc
clear all

w = 1.02;           % Wheel base [m]
t = 0.08;           % trail [m]
alph = atan(3);    % Head angle
g = 9.81;           % Gravity [N/kg]

% Rear wheel
R_rw = 0.3;         % Radius [m]
m_rw = 2;           % Mass [kg]
A_xx = 0.06;        % Mass moments of inertia [km^2]
A_yy = 0.12;        % Mass moments of inertia [km^2]
A_zz = 0.06;        % Mass moments of inertia [km^2]

% Rear frame
x_rf = 0.3;         % Pos. centre of mass [m]
y_rf = 0;           % Pos. centre of mass [m]
z_rf = -0.9;        % Pos. centre of mass [m]
m_rf = 85;          % Mass [kg]
B_xx = 9.2;         % Mass moments of interia [km^2]
B_xz = 2.4;         % Mass moments of interia [km^2]
B_yy = 11;          % Mass moments of interia [km^2]
B_zz = 2.8;         % Mass moments of interia [km^2]

MMI_rf = [B_xx 0 B_xz;...
          0   B_yy 0;
          B_xz 0 B_zz];
      
% Front frame
x_ff = 0.9;         % Pos. centre of mass [m]
y_ff = 0;           % Pos. centre of mass [m]
z_ff = -0.7;        % Pos. centre of mass [m]
m_ff = 4;           % Mass [kg]
C_xx = 0.0546;      % Mass moments of interia [km^2]
C_xz = -0.0162;     % Mass moments of interia [km^2]
C_yy = 0.06;        % Mass moments of interia [km^2]
C_zz = 0.0114;      % Mass moments of interia [km^2]

MMI_ff = [C_xx 0 C_xz;...
          0   C_yy 0;
          C_xz 0 C_zz];

% Front wheel
R_fw = 0.35;        % Radius [m]
m_fw = 3;           % Mass [kg]
D_xx = 0.14;        % Mass moments of inertia [km^2]
D_yy = 0.28;        % Mass moments of inertia [km^2]
D_zz = 0.14;        % Mass moments of inertia [km^2]


save Benchmark_model 