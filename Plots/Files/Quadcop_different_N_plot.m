clear all
close all

t = linspace(0,6,61);

load('Quadcop_different_N\Quadcop_states_v_2_Q_100000_R_1_N_8.mat')
x_Q_100000_R_1_N_8 = x;

load('Quadcop_different_N\Quadcop_states_v_2_Q_100000_R_1_N_16.mat')
x_Q_100000_R_1_N_16 = x;

load('Quadcop_different_N\Quadcop_states_v_2_Q_100000_R_1_N_20.mat')
x_Q_100000_R_1_N_20 = x;

load('Quadcop_different_N\Quadcop_states_v_2_Q_100000_R_1_N_25.mat')
x_Q_100000_R_1_N_25 = x;

figure(1)
hold on
stairs(t,x_Q_100000_R_1_N_8(1,:))
stairs(t,x_Q_100000_R_1_N_16(1,:))
stairs(t,x_Q_100000_R_1_N_20(1,:))
stairs(t,x_Q_100000_R_1_N_25(1,:))
legend('N = 8','N = 16','N = 20','N = 25')
hold off

%% 
load('Quadcop_different_N_new_constr\Quadcop_states_v_2_Q_100000_R_1_N_8.mat')
x_Q_100000_R_1_N_8 = x;

load('Quadcop_different_N_new_constr\Quadcop_states_v_2_Q_100000_R_1_N_16.mat')
x_Q_100000_R_1_N_16 = x;

% load('Quadcop_different_N_new_constr\Quadcop_states_v_2_Q_100000_R_1_N_20.mat')
% x_Q_100000_R_1_N_20 = x;

load('Quadcop_different_N_new_constr\Quadcop_states_v_2_Q_100000_R_1_N_25.mat')
x_Q_100000_R_1_N_25 = x;

figure(2)
hold on
stairs(t,x_Q_100000_R_1_N_8(1,:))
stairs(t,x_Q_100000_R_1_N_16(1,:))
% stairs(t,x_Q_100000_R_1_N_20(1,:))
stairs(t,x_Q_100000_R_1_N_25(1,:))
legend('N = 8','N = 16','N = 20','N = 25')
hold off