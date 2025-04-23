% ------------------------------------------------------------
% Solving the 1D Wave Equation using Separation of Variables
% ------------------------------------------------------------
%
% Wave Equation:
%     ∂²u/∂t² = c² * ∂²u/∂x²
%
% Domain:
%     0 < x < 1,   t > 0
%
% Boundary Conditions:
%     u(0, t) = 0           % Fixed end at x = 0
%     u(1, t) = 0           % Fixed end at x = 1
%
% Initial Conditions:
%     u(x, 0) = sin(pi * x)        % Initial displacement
%     ∂u/∂t (x, 0) = 0             % Initial velocity
%
% Solution method: Separation of variables and Fourier sine series


% Parameters
L = 1;           % Length of the string
c = 1;           % Wave speed
T = 2;           % Final time
N = 100;         % Number of spatial points
M = 100;         % Number of time points
num_terms = 50;  % Number of terms in Fourier series

% Spatial and temporal grid
x = linspace(0, L, N);
t = linspace(0, T, M);
[X, T_grid] = meshgrid(x, t);

% Initial condition functions
f = @(x) sin(pi*x);             % Initial displacement
g = @(x) 0.*x;                  % Initial velocity

% Preallocate solution matrix
u = zeros(M, N);

% Build the solution using separation of variables
for n = 1:num_terms
    lambda_n = n*pi/L;
    A_n = 2/L * trapz(x, f(x).*sin(lambda_n*x));     % Fourier coefficient for f(x)
    B_n = 2/(lambda_n*c*L) * trapz(x, g(x).*sin(lambda_n*x)); % Fourier coeff for g(x)

    u = u + (A_n * cos(c*lambda_n*T_grid) + B_n * sin(c*lambda_n*T_grid)) .* sin(lambda_n*X);
end


% Plot the solution as animation
for k = 1:M
    plot(x, u(k, :), 'b', 'LineWidth', 2);
    axis([0 L -1.5 1.5]);
    xlabel('x');
    ylabel('u(x,t)');
    title(['Time t = ', num2str(t(k), '%.2f')]);
    grid on;
    drawnow;
end