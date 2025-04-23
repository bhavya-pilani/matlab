function linear_phase_portrait()
    clc; clear; close all;

    % ===== Define Linear System dx/dt = A * x =====
    A = [-1 0;0 1]; % Example system matrix

    % Define the function for ODE solver
    f = @(t, x) A * x; % dx/dt = A * x

    % ===== Generate Phase Portrait =====
    [x, y] = meshgrid(-5:0.5:5, -5:0.5:5); % Grid points
    u = zeros(size(x)); 
    v = zeros(size(y));

    for i = 1:numel(x)
        dx = f(0, [x(i); y(i)]); % Compute derivatives at each point
        u(i) = dx(1);
        v(i) = dx(2);
    end

    % Plot vector field
    figure; hold on;
    quiver(x, y, u, v, 'b'); % Vector field arrows

    % ===== Solve ODE for Multiple Initial Conditions =====
    tspan = [0, 5]; % Time span
    init_conditions = [0.2 0.5; -0.2 0.5; -0.2 -0.5; 0.2 -0.5;
                       0.5 0.2; -0.5 0.2; -0.5 -0.2; 0.5 -0.2];

    for i = 1:size(init_conditions, 1)
        [T, X] = ode45(f, tspan, init_conditions(i, :)');
        plot(X(:,1), X(:,2), 'r', 'LineWidth', 1.5); % Trajectories
    end

    % Formatting
    title('Phase Portrait for dx/dt = A * x');
    xlabel('x_1'); ylabel('x_2');
    grid on; axis equal;
    legend('Vector Field', 'Trajectories');
end