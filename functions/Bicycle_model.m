function [Bic_sys, Matrices] = Bicycle_model(v)

load('Benchmark_model.mat');

m_t = m_ff + m_rf + m_rw + m_fw;                % Total mass [kg]
x_t = (x_rf*m_rf + x_ff*m_ff + w*m_fw)/m_t;
z_t = (-R_rw*m_rw + z_rf*m_rf + z_ff*m_ff - R_fw*m_fw)/m_t;

% Products of Inertia
T_xx = A_xx + B_xx + C_xx + D_xx + m_rw*R_rw^2 + m_rf*z_rf^2 + m_ff*z_ff^2 + m_fw*R_fw^2;
T_xz = B_xz + C_xz - m_rf*x_rf*z_rf - m_ff*x_ff*z_ff + m_fw*w*R_fw;
T_zz = A_zz + B_zz + C_zz + D_zz + m_rf*x_rf^2 + m_ff*x_ff^2 + m_fw*w^2;

m_f = m_ff + m_fw;
x_f = (x_ff*m_ff + w*m_fw)/m_f;
z_f = (z_ff*m_ff - R_fw*m_fw)/m_f;

% Products of Intertia
F_xx = C_xx + D_xx + m_ff*(z_ff - z_f)^2 + m_fw*(R_fw + z_f)^2;
F_xz = C_xz - m_ff*(x_ff - x_f)*(z_ff - z_f) + m_fw*(w-x_f)*(R_fw + z_f);
F_zz = C_zz + D_zz + m_ff*(x_ff - x_f)^2 + m_fw*(w-x_f)^2;

% Angle of steering axis
lab = pi/2 - alph;
LABDA = [sin(lab); 0; cos(lab)];

% perpendicular distance centre of mass of front assembly ahead of steering
% axis
u = (x_f - w - t)*cos(lab) - z_f*sin(lab);

% Products of Inertia
F_lablab = m_f*u^2 + F_xx*sin(lab)^2 + 2*F_xz*sin(lab)*cos(lab) + F_zz*cos(lab)^2;
F_labx = -m_f*u*z_f + F_xx*sin(lab) + F_xz*cos(lab);
F_labz = m_f*u*x_f + F_xz*sin(lab) + F_zz*cos(lab);

% ratio of mechanical trail to wheelbase
f = t*cos(lab)/w;

% Angular momentum
S_r = A_yy/R_rw;
S_f = D_yy/R_fw;
S_t = S_r + S_f;

% Static Moment term
S_u = m_f*u + f*m_t*x_t;

% Mass matrix
M(1,1) = T_xx;
M(1,2) = F_labx + f*T_xz;
M(2,1) = M(1,2);
M(2,2) = F_lablab + 2*f*F_labz + f^2*T_zz;
Matrices.M = M;

% Velocity independent Stifness Matrix
K0(1,1) = g*m_t*z_t;
K0(1,2) = -g*S_u;
K0(2,1) = K0(1,2);
K0(2,2) = -g*S_u*sin(lab);
Matrices.K0 = K0;

% Velocity dependent Stifness Matrix
K2(1,1) = 0;
K2(1,2) = (S_t - m_t*z_t)*cos(lab)/w;
K2(2,1) = 0;
K2(2,2) = (S_u + S_f*sin(lab))*cos(lab)/w;
Matrices.K2 = K2;

% Damping Matrix
C1(1,1) = 0;
C1(1,2) = f*S_t + S_f*cos(lab) + T_xz*cos(lab)/w - f*m_t*z_t;
C1(2,1) = -(f*S_t + S_f*cos(lab));
C1(2,2) = F_labz*cos(lab)/w + f*(S_u + T_zz*cos(lab)/w);
Matrices.C1 = C1;

A_matrix = [zeros(2,2)                 eye(2);...
            -inv(M)*(K0 + K2*v^2)      -inv(M)*C1*v];

B_matrix = [zeros(2,2);...
            inv(M)];
        
C_matrix = [eye(2) zeros(2,2)];

% Continuous time state-space
Bic_sys = ss(A_matrix,B_matrix,C_matrix,0);


end
