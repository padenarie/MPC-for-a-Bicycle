close all
clear all

steps = 400;
t = linspace(0,40,steps+1);

load("MPC_output_different_N\MPC_output_N_8")
x_N_8 = x;

load("MPC_output_different_N\MPC_output_N_16")
x_N_16 = x;

load("MPC_output_different_N\MPC_output_N_25")
x_N_25 = x;


hold on
grid on
stairs(t,x_N_8(1,:))
stairs(t,x_N_16(1,:))
stairs(t,x_N_25(1,:))
xlabel('Time [s]')
ylabel('Roll angle $\phi$ [rad]','Interpreter','latex')


