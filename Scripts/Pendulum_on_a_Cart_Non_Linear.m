clear;

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
Q = [B A*B A*A*B A*A*A*B];
if rank(Q) == min(size(Q))
    disp('Controllable')
end

% Define cost function elements for LQR
Q_u = 1;
Q_x = [100 0 0 0; 0 1 0 0; 0 0 10 0; 0 0 0 1];

% Perform LQR to obtain Feedback Gain K
K=lqr(A,B,Q_x,Q_u);

% Compute Feedforward Gain K_f
K_f = -inv(C*inv(A-B*K)*B);

% Define Simulation Parameters
x0 = [0;0;0;0];
dt = 0.01;
tspan = 0:dt:200;

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
function dxdt = odefun(t, x, K, y_d, K_f)
input = -K * x + K_f*y_d(t);
dxdt = [x(3,:);...
        x(4,:);...
       -(-sin(x(2,:))*x(4,:)^2 + input - 3*x(3,:) + cos(x(2,:))*sin(x(2,:)))/(cos(x(2,:))^2 - 2);...
       -(-cos(x(2,:))*sin(x(2,:))*x(4,:)^2 + 2*sin(x(2,:)) + input*cos(x(2,:)) - 3*x(3,:)*cos(x(2,:)))/(cos(x(2,:))^2 - 2)];
end

% Solve the ODE using ode45
[t, x] = ode45(@(t, x) odefun(t, x, K, y_d, K_f), tspan, x0);

% Apply Output Matrix
y = C * x';
v = y_d(tspan);

% Plot the position for tracking
figure;
plot(t, y, "LineWidth",2, "Color",[1, 0, 0]);
hold on
plot(t, v, "LineWidth",2, "Color",[0, 0, 1])
title("Tracking Controller Implemented");
xlabel('Time (s)');
ylabel('Position (in)');
legend('System Position', 'Desired Position')
hold off