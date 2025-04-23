% Parameters
L = 10;           % Length of the rod
T_final = 2;      % Final time
alpha = 0.01;     % Thermal diffusivity
Nx = 50;          % Number of spatial points
Nt = 500;         % Number of time steps
x = linspace(0, L, Nx);  % Spatial grid
t = linspace(0, T_final, Nt); % Time grid

% Initial condition: u(x, 0) = sin(pi*x)
u_initial = sin(pi * x);  % Initial temperature distribution

% Eigenvalue problem: X''(x) + lambda X(x) = 0
% The solution of this equation is X_n(x) = sin(n*pi*x / L), where n is an integer
n_values = 1:10;  % First few eigenvalues
lambda_values = (n_values * pi / L).^2;  % Corresponding eigenvalues

% Time-dependent solution: T_n(t) = exp(-alpha * lambda_n * t)
% The general solution is a sum of all these modes
u = zeros(Nx, Nt);  % Solution matrix (space x time)

for n = 1:length(n_values)
    lambda_n = lambda_values(n);
    X_n = sin(n * pi * x / L);  % Spatial part
    T_n = exp(-alpha * lambda_n * t);  % Temporal part
    
    % Add contribution of the n-th mode to the solution
    % Using Fourier coefficients for initial condition
    A_n = trapz(x, u_initial .* X_n) * 2 / L;  % Fourier coefficient
    u = u + A_n * X_n' * T_n;  % Build the solution
end

% Display the values of u(x,t) at specific points (for example, at t = T_final)
disp('Temperature values u(x,t) at specific time points:');

% Example: Display u(x,t) at the final time step
disp('Temperature at t = T_final:');
disp(u(:, Nt));  % Display the temperature at the final time (t = T_final)

% Optionally, display u(x,t) at other time steps as well (for example, at T/2)
disp('Temperature at t = T_final/2:');
disp(u(:, round(Nt/2)));  % Display the temperature at middle time step (T/2)

% Optionally, display u(x,t) at the first time step (initial condition)
disp('Initial temperature distribution at t = 0:');
disp(u(:, 1));  % Display the initial temperature at t = 0

% Plot the results for different times
figure;
subplot(3,1,1);
plot(x, u(:,1), 'r', 'LineWidth', 2);
xlabel('Position (x)');
ylabel('Temperature (u)');
title('Temperature Distribution at t = 0 (Initial Condition)');

subplot(3,1,2);
plot(x, u(:, round(Nt/2)), 'g', 'LineWidth', 2);
xlabel('Position (x)');
ylabel('Temperature (u)');
title(['Temperature Distribution at t = ', num2str(T_final/2)]);

subplot(3,1,3);
plot(x, u(:, Nt), 'b', 'LineWidth', 2);
xlabel('Position (x)');
ylabel('Temperature (u)');
title(['Temperature Distribution at t = ', num2str(T_final)]);

% Create a 3D surface plot to visualize the entire solution
figure;
surf(t, x, u);
xlabel('Time (t)');
ylabel('Position (x)');
zlabel('Temperature (u)');
title('Heat Equation Solution Using Separation of Variables');