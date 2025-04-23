% ------------------------------------------------------------
% Solving Laplace's Equation using Separation of Variables
% ------------------------------------------------------------
%
% Laplace's Equation:
%     ∂²u/∂x² + ∂²u/∂y² = 0
%
% Domain:
%     0 < x < 1,
%     0 < y < 1
%
% Boundary Conditions:
%     u(x, 0) = 0          % Bottom edge
%     u(x, 1) = sin(pi*x)  % Top edge
%     u(0, y) = 0          % Left edge
%     u(1, y) = 0          % Right edge
%
% Solution method: Separation of variables with Fourier sine series


% Domain setup
a = 1; b = 1;        % Domain size
Nx = 50; Ny = 50;    % Number of grid points
x = linspace(0, a, Nx);
y = linspace(0, b, Ny);
[X, Y] = meshgrid(x, y);

% Number of terms in the series
Nterms = 50;

% Initialize solution
U = zeros(Ny, Nx);  % meshgrid format: rows = y, cols = x

% Solve using separation of variables
for n = 1:Nterms
    lambda_n = n * pi;

    % Compute Fourier coefficient Bn numerically
    integrand = f(x) .* sin(lambda_n * x);    % f(x) * sin(n*pi*x)
    Bn = 2 * trapz(x, integrand) / sinh(lambda_n * b);

    % Add term to solution
    U = U + Bn * sin(lambda_n * X) .* sinh(lambda_n * Y);
end
disp(U)

% Plot the solution
surf(X, Y, U, 'EdgeColor', 'none');
xlabel('x');
ylabel('y');
zlabel('u(x,y)');
title('Solution of Laplace Equation using Separation of Variables');
colorbar;
view(45, 30);