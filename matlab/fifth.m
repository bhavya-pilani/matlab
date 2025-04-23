% MATLAB Code to solve the given PDE using Lagrange's Method

% Define symbolic variables
syms x y z p q c1 c2 U V

% Define the given PDE in the form P*p + Q*q = R
lhs = ((y^2)*z/x)*p + (x)*z*q;  % LHS side
rhs = y^2;                        % RHS side

% Extract coefficients P, Q, and R explicitly
P = coeffs(lhs, p, 'All');  
Q = coeffs(lhs, q, 'All');  
R = rhs;   % Extracting R

% Ensure P and Q are scalars
if length(P) > 1, P = P(1); end
if length(Q) > 1, Q = Q(1); end

% Display extracted coefficients
disp('P = '), disp(P);
disp('Q = '), disp(Q);
disp('R = '), disp(R);

% solve Characteristic Equations
% dx/P = dy/Q = dz/R
% Solve for U (integrating first subsidiary equation)
U = int(Q, x) - int(P, y) + c1; 
U = simplify(U);
disp('U = '), disp(U);

% Solve for V (integrating second subsidiary equation)
V = int(R, x) - int(P, z) + c2;
V = simplify(V);
disp('V = '), disp(V);

% Display final solution
disp('The general integral is given by: ');
fprintf('f(U, V) = f(%s, %s) = 0\n', char(U), char(V));