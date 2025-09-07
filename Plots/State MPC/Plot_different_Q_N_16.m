clear all
close all

t = linspace(0,6,61);

load("Data\MPC_reg_states_N_16\States_N_16_v_2_Q_1_R_1")
x_1 = x;
load("Data\MPC_reg_states_N_16\States_N_16_v_2_Q_10_R_1")
x_10 = x;
load("Data\MPC_reg_states_N_16\States_N_16_v_2_Q_100_R_1")
x_100 = x;
load("Data\MPC_reg_states_N_16\States_N_16_v_2_Q_1000_R_1")
x_1000 = x;
load("Data\MPC_reg_states_N_16\States_N_16_v_2_Q_10000_R_1")
x_10000 = x;


hold on
grid on
stairs(t,x_1(1,:))
% stairs(t,x_10(1,:))
stairs(t,x_100(1,:))
stairs(t,x_1000(1,:))
stairs(t,x_10000(1,:))
ylabel('$\phi$ [rad]','interpreter','latex')
xlabel('Time [s]')
legend('Q = 1, R = 1','Q = 100, R = 1','Q = 1000, R = 1','Q = 10000, R = 1')
hold off