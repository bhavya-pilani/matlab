%Write a program in MATLAB to solve non-linear PDE using Charpit’s method. 
%NAME : DINKAR SHARMA
%ROLL NO. : 23/MC/050

% Define symbols
syms x(t) y(t) z(t) p(t) q(t) X Y

% Step 1: Define the PDE
% Example PDE: p^2 + q^2 - z = 0
F = p^2 + q^2 - z;

% Step 2: Charpit's Equations
dx_dt = diff(F, p);                 % dx/dt = 2*p
dy_dt = diff(F, q);                 % dy/dt = 2*q
dz_dt = p*diff(F, p) + q*diff(F, q); % dz/dt = 2*p^2 + 2*q^2
dp_dt = -diff(F, x);                % dp/dt = 0 (since F has no x)
dq_dt = -diff(F, y);                % dq/dt = 0 (since F has no y)

% Step 3: Solve the system of ODEs using dsolve
sol = dsolve([diff(x, t) == dx_dt, ...
              diff(y, t) == dy_dt, ...
              diff(z, t) == dz_dt, ...
              diff(p, t) == dp_dt, ...
              diff(q, t) == dq_dt]);

% Step 4: Display Results
disp('Charpit''s solution:');
disp('x(t) ='); disp(sol.x);
disp('y(t) ='); disp(sol.y);
disp('z(t) ='); disp(sol.z);
disp('p(t) ='); disp(sol.p);
disp('q(t) ='); disp(sol.q);

% Step 5: Substitute back to get general solution if needed
% Eliminate parameters to get solution in form z = f(x, y)

% Example: assume initial conditions to simplify
% Assume p = a, q = b constants (since dp/dt = dq/dt = 0)
syms a b t
x_t = 2*a*t + X;     % from dx/dt = 2*p = 2*a
y_t = 2*b*t + Y;     % from dy/dt = 2*q = 2*b
z_t = 2*a^2*t + 2*b^2*t;  % from dz/dt = 2*p^2 + 2*q^2

% Eliminate t from x_t and y_t
t1 = (x_t - X)/(2*a);
t2 = (y_t - Y)/(2*b);

% From either t1 or t2, substitute in z_t
% For consistency, use t1:
z_gen = subs(z_t, t, t1);
z_gen = simplify(z_gen);

disp('General solution z(x, y):');
disp(z_gen);