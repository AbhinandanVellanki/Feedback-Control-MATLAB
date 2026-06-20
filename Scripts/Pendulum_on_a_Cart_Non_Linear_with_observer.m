clear; % Remove all previous variables
clc;
close all;

% Define symbolic variables
syms xc_dot phi_dot xc_ddot phi_ddot phi u real
syms gamma beta alpha mu D real

% Define the matrix M
M = [gamma, -beta*cos(phi);
    -beta*cos(phi), alpha];

% Define the matrix N based on the equations of motion
N = [u - beta*phi_dot^2*sin(phi) - mu*xc_dot;
    D*sin(phi)];

% Solve for state update variables [xc_ddot, phi_ddot]
state_updates = M \ N;

% Substitute numerical values
gamma_val = 2;
beta_val = 1;
alpha_val = 1;
D_val = 1;
mu_val = 3;

% Substitute into accelerations
state_updates = subs(state_updates, ...
    {gamma, beta, alpha, D, mu}, ...
    {gamma_val, beta_val, alpha_val, D_val, mu_val});

%Create the non-linear state space model
F = [xc_dot;phi_dot;state_updates(1); state_updates(2)];

% Define the linear system
A = [0 0 1 0; 0 0 0 1; 0 1 -3 0; 0 2 -3 0];
B = [0; 0; 1; 1];

% Output matrix to convert output into position in inches
C = [39.37 0 0 0];

% Controllability Test
Q_c = [B A*B A*A*B A*A*A*B];
if rank(Q_c) == min(size(Q_c))
    disp('Controllable')
else
    return
end

% Observability Test
Q_o = [C C*A C*A*A C*A*A*A]';
if rank(Q_o) == min(size(Q_o))
    disp('Observable')
else
    return
end

% Define cost function elements for LQR
Q_u = 1;
Q_x = [100 0 0 0; 0 1 0 0; 0 0 100 0; 0 0 0 1];

% Perform LQR to obtain Feedback Gain K
[K_c, s, cl_eigs] =lqr(A,B,Q_x,Q_u);

% Obtain observer gain K_o using eigenvalue placement
K_o = place(A',C',4*cl_eigs)';

% Compute Feedforward Gain K_f
K_f = -inv(C*inv(A-B*K_c)*B);

% Define Simulation Parameters
x_actual0 = [0;0;0;0]; % True initial state
x_estimated0 = [0.01;0.01;-0.03;0.01]; % Initial State Estimate
dt = 0.01;
tspan = 0:dt:20;

% Define desired output
period = 100;  % period of the wave
amplitude = 40;  % amplitude of the wave

% Create the desired output function
function y = square_wave(t, period, amplitude)
    % Calculate the square wave
    normalized_t = mod(t, period) / period;
    y = amplitude * (normalized_t < 0.5) - amplitude/2;
end

% Create function handle for desired output function
y_d = @(t) square_wave(t, period, amplitude);

% Defining the ODE equation
function dxdt = odefun(t, x, K_c, K_f, K_o, C, y_d)
    % Get current state and estimate
    x_actual = x(1:4);
    x_estimated = x(5:8);
    
    % True Output - usually measured by sensor
    y_actual = C*x_actual;
    
    % Input based on estimated state and reference signal
    input = -K_c * x_estimated + K_f*y_d(t);
    
    % Actual system dynamics
    dx_actual = [x_actual(3,:);...
        x_actual(4,:);...
       -(-sin(x_actual(2,:))*x_actual(4,:)^2 + input - 3*x_actual(3,:) + cos(x_actual(2,:))*sin(x_actual(2,:)))/(cos(x_actual(2,:))^2 - 2);...
       -(-cos(x_actual(2,:))*sin(x_actual(2,:))*x_actual(4,:)^2 + 2*sin(x_actual(2,:)) + input*cos(x_actual(2,:)) - 3*x_actual(3,:)*cos(x_actual(2,:)))/(cos(x_actual(2,:))^2 - 2)];
    
    % Observer (estimated system dynamics)
    dx_estimated = [x_estimated(3,:);...
        x_estimated(4,:);...
       -(-sin(x_estimated(2,:))*x_estimated(4,:)^2 + input - 3*x_estimated(3,:) + cos(x_estimated(2,:))*sin(x_estimated(2,:)))/(cos(x_estimated(2,:))^2 - 2);...
       -(-cos(x_estimated(2,:))*sin(x_estimated(2,:))*x_estimated(4,:)^2 + 2*sin(x_estimated(2,:)) + input*cos(x_estimated(2,:)) - 3*x_estimated(3,:)*cos(x_estimated(2,:)))/(cos(x_estimated(2,:))^2 - 2)];

    % Add weighted difference in true output and output based on estimated state
    dx_estimated = dx_estimated + K_o*(y_actual - C*x_estimated);
    
    % update of state and state estimate
    dxdt = [dx_actual; dx_estimated];
end

% Solve the ODE for the actual system using ode45
[t, x] = ode45(@(t, x) odefun(t, x, K_c, K_f, K_o, C, y_d), tspan, [x_actual0; x_estimated0]);

% Extract true states and estimated states
x_true = x(:, 1:4);
x_hat = x(:, 5:8);

% Convert position to inches in state and estimated state
x_true (:, 1) = (C(1) * x_true(:, 1));
x_hat (:, 1) = (C(1) * x_hat(:, 1));


% State variable info
state_variables = ["Position", "Pendulum Angle", "Velocity", "Pendulum Velocity"];
units = ["Inches", "Radians", "Metre/Second", "Radian/Second"];

% Create a figure with 4 subplots
figure;
for i = 1:4
    subplot(4, 1, i);
    plot(t, x_true(:, i), 'b', 'LineWidth', 2);
    hold on;
    plot(t, x_hat(:, i), 'r', 'LineWidth', 1);
    title([state_variables(i)]);
    xlabel('Time (s)');
    ylabel(units(i));
    legend(['x', num2str(i)], ['x', num2str(i), '_{hat}'], 'Location', 'northeast');
    grid on;
end

% Set common title
sgtitle('True State vs Estimated State');

% Plot estimation error
% figure;
% plot(t, x_true - x_hat);
% xlabel('Time');
% ylabel('Estimation Error');
% legend('e1', 'e2', 'e3', 'e4');
% title('Estimation Error');
% grid on;