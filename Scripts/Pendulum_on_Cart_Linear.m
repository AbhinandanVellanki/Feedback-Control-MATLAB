clear;

% Define the linear system
A = [0 0 1 0; 0 0 0 1; 0 1 -3 0; 0 2 -3 0];
B = [0; 0; 1; 1];
C = eye(4);
D = zeros(4,1);

% Controllability Test
Q = [B A*B A*A*B A*A*A*B];
if rank(Q) == min(size(Q))
    disp('Controllable')
end

% Define cost function elements for LQR
Q_u = 10;
Q_x = [1 0 0 0; 0 5 0 0; 0 0 1 0; 0 0 0 5];


% Perform LQR to obtain K
K=lqr(A,B,Q_x,Q_u);

% Define Simulation Parameters

x0 = [0;1.1;0;0];
tspan = 0:0.01:30;
dt = 0.01;

% Define the input function u(t)
u = @(t) -K * x; % Feedback control

% Define the ODE function with feedback
function dxdt = odefun(t, x, A, B, K)
    u = -K * x; % Calculate control input
    dxdt = A * x + B * u;
end


% Solve the ODE using ode45
[t, x] = ode45(@(t, x) odefun(t, x, A, B, K), tspan, x0);

% State variable info
state_variables = ["Position", "Pendulum Angle", "Velocity", "Pendulum Velocity"];
units = ["Metres", "Radians", "Metre/Second", "Radian/Second"];

% Plot the results
for i = 1:length(state_variables)
    subplot(length(state_variables),1,i);
    plot(t, x(:,i), 'LineWidth', 2); % Thick lines
    title([state_variables(i)]);
    xlabel('Time (s)');
    ylabel(units(i));
    grid on;
end