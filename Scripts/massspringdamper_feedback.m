clear

m = 1;
u = 0.5;
k = 5;
dt = 0.1;
t = 0:dt:10;

x = zeros(2, length(t));
x(:,1)=[1; 1];

A = [0 1; -k/m -u/m];
B = [0; 1/m];
poles = [-1+1i, -1-1i];
K = place(A, B, poles);
disp(size(K));
disp(size(B));
Acl = A-B*K;


for i = 2:length(t)
    x(:,i) = expm(Acl*dt)*x(:,i-1) ;
end

% State variable info
state_variables = ["Position", "Velcocity"];
units = ["Metres", "Metres/Second"];

% Plot the results
for i = 1:length(state_variables)
    subplot(length(state_variables),1,i);
    plot(t, x(i,:), 'LineWidth', 2); % Thick lines
    title([state_variables(i)]);
    xlabel('Time (s)');
    ylabel(units(i));
    grid on;
end