%%% Control Engineering - Tracking and regulation 
%%% Coursework PartB

%% B1
% Define system parameters
M=1;
L=1;
F=1;
g=9.81;

omega=1;
alpha=1;

% Define system matrix from partA
A=[0,1,0,0;
   0,-F/M,0,0;
   0,0,0,1;
   0,F/(M*L),g/L,0];

B=[0;1/M;0;-1/(M*L)];

% Disturbance input matrix for d1 
P_dist=[0;1/M;0;-1/(M*L)];
% Output matrix for regulation. Error e = Cx + Qw.
Ce=[1,0,0,0];

% Define Exosystem Matrix from A5
S=[0,0,0;
   0,0,omega;
   0,-omega,0];

% Expanded P matrix 
P_exp=[P_dist, zeros(4, 2)]; 

% Q matrix to extract reference from w 
Q=[0, -1, 0];

% Design State Feedback Gain K. We need to place poles for (A + B*K) in the left half plane.
% We choose poles that are faster than the open-loop dynamics but not too aggressive.
desired_poles=[-2,-2.5,-3,-3.5]; 
K_reverse = place(A, B, desired_poles);
K = -K_reverse;

disp('State Feedback Gain K:');
disp(K);

% Verify Stability
eigenvalues_closed_loop = eig(A + B*K);
disp('Closed-loop Eigenvalues:');
disp(eigenvalues_closed_loop);

% Column 1 Constant Disturbance
M_const=[A,B;Ce,0];
RHS_const=[-P_exp(:,1);-Q(1)];
sol_1=M_const \ RHS_const;
Pi_1=sol_1(1:4);
Gamma_1=sol_1(5);

% Columns 2 and 3 Sinusoidal Reference
n = 4;
BigMat=[A,omega*eye(n),B,zeros(n,1);
        -omega*eye(n),A,zeros(n,1), B;
        Ce,zeros(1,n),0,0;
        zeros(1,n),Ce,0,0];
BigRHS=[-P_exp(:,2);-P_exp(:,3);-Q(2);-Q(3)];
sol_sin=BigMat \ BigRHS;

Pi_2=sol_sin(1:4);
Pi_3=sol_sin(5:8);
Gamma_2= sol_sin(9);
Gamma_3=sol_sin(10);

Pi=[Pi_1, Pi_2, Pi_3];
Gamma=[Gamma_1, Gamma_2, Gamma_3];

% Calculate L
L = Gamma - K * Pi;

fprintf('Feedforward Gain L:\n');
disp(L);

%% B2
% Simulation setup 
T_final=100;           % 2 periods of disturbance
tspan=[0 T_final];
x0=[0; 0; 0; 0];       

%Define Signal Functions
% Disturbance d1(t)
get_d1=@(t) 0.5 * square(2*pi*t/50);

% Reference d2(t)
get_d2=@(t) 1 * sin(omega*t);
get_d2_dot=@(t) 1 * omega * cos(omega*t); 

% dx/dt = Ax + Bu + P*d1, u = Kx + Lw
dynamics=@(t, x) linear_system_dynamics(t, x, A, B, P_dist, K, L, get_d1, get_d2, get_d2_dot, omega);

% Run simulation
[t_sim, x_sim]=ode45(dynamics, tspan, x0);

% Recompute u, d1, d2 for all time points to plot them for plotting
num_steps=length(t_sim);
u_sim=zeros(num_steps, 1);
d1_sim=zeros(num_steps, 1);
d2_sim=zeros(num_steps, 1);
s_sim=x_sim(:, 1);  % Cart position

for i = 1:num_steps
    ti = t_sim(i);
    xi = x_sim(i, :)';
    
    % Reconstruct Exosystem State w(t)
    d1_val= get_d1(ti);
    d2_val= get_d2(ti);
    d2_dot_val = get_d2_dot(ti);
    
    % w = [d1; d2; d2_dot/omega]
    w_val = [d1_val; d2_val; d2_dot_val/omega];
    
    % Control Input u = Kx + Lw
    u_sim(i) = K * xi + L * w_val;
    
    d1_sim(i) =d1_val;
    d2_sim(i) = d2_val;
end

% Plotting
figure(1); clf;

% Output s(t) vs Reference d2(t)
subplot(3,1,1);
plot(t_sim, s_sim, 'b', 'LineWidth',1); 
hold on;
plot(t_sim, d2_sim, 'r--', 'LineWidth',1);
title('B2) Output Tracking: Cart Position s(t) vs Reference d_2(t)');
xlabel('Time (s)'); 
ylabel('Position (m)');
legend('Output s(t)', 'Reference d_2(t)');
grid on;

% Control Input u(t)
subplot(3,1,2);
plot(t_sim, u_sim, 'k', 'LineWidth',1);
title('Control Input \mu(t)');
xlabel('Time (s)'); ylabel('Force (N)');
grid on;

% Disturbance d1(t)
subplot(3,1,3);
plot(t_sim, d1_sim, 'm', 'LineWidth',1);
title('External Disturbance d_1(t) (Square Wave)');
xlabel('Time (s)'); ylabel('Force (N)');
ylim([-0.6 0.6]);
grid on;

%% B3
L_len=1; 

% Simulation for Nonlinear System
dynamics_nonlin=@(t, x) nonlinear_system_dynamics(t, x, M, L_len, F, g, K, L, get_d1, get_d2, get_d2_dot, omega);

[t_nl, x_nl]=ode45(dynamics_nonlin, tspan, x0);

d1_nl=get_d1(t_nl);
d2_nl=get_d2(t_nl);
d2_dot_nl=get_d2_dot(t_nl);

% W_nl: [d1; d2; d2_dot/omega]
W_nl=[d1_nl'; d2_nl'; d2_dot_nl'/omega];

% Calculate Control Input u 
u_nl=(K * x_nl' + L * W_nl)';

% Linear vs Nonlinear
figure(2); clf;

% Output s(t) Comparison
subplot(3,1,1);
plot(t_sim, x_sim(:,1), 'b', 'LineWidth', 1); 
hold on; 
plot(t_nl, x_nl(:,1), 'r--', 'LineWidth', 1);          % Nonlinear
plot(t_sim, d2_sim, 'g:', 'LineWidth', 1);               % Reference
title('B3) Comparison: Linear (Blue) vs Nonlinear (Red) Output s(t)');
xlabel('Time (s)'); ylabel('Position (m)');
legend('Linear Model', 'Nonlinear Model', 'Reference');
grid on;

% Pendulum Angle phi(t)
subplot(3,1,2);
plot(t_sim, x_sim(:,3), 'b', 'LineWidth', 1); 
hold on;
plot(t_nl, x_nl(:,3), 'r--', 'LineWidth', 1);
title('Pendulum Angle \phi(t) (rad)');
xlabel('Time (s)'); ylabel('Angle (rad)');
legend('Linear Model', 'Nonlinear Model');
grid on;

% Control Input u(t)
subplot(3,1,3);
plot(t_sim, u_sim, 'b', 'LineWidth', 1); 
hold on;
plot(t_nl, u_nl, 'r--', 'LineWidth', 1);
title('Control Input \mu(t)');
xlabel('Time (s)'); ylabel('Force (N)');
grid on;

% RMS
error_nl=x_nl(:,1) - d2_nl;
rms_error=sqrt(mean(error_nl.^2));
fprintf('RMS Tracking Error (Nonlinear): %.4f\n', rms_error);

