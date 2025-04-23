% Clear the screen
clc;

% Define x as a symbolic variable
syms x;

% Input Section
y = input('Enter the nonlinear equation: ');
x0 = input('Enter the initial guess : ');
e = input('Enter the tolerable error : ');
n = input('Enter the maximum number of iterations: ');

% Compute the derivative
dy = diff(y, x);

% Initialize
iter = 0;
fprintf('Iter\t x\t\t f(x)\n');

% Newton-Raphson Iteration
while iter < n
    f_val = eval(subs(y, x, x0));
    df_val = eval(subs(dy, x, x0));
    
    if df_val == 0
        disp('Derivative is zero. Method fails.');
        break;
    end
    
    x1 = x0 - f_val / df_val;
    
    fprintf('%d\t%f\t%f\n', iter, x0, f_val);
    
    if abs(x1 - x0) < e
        break;
    end
    
    x0 = x1;
    iter = iter + 1;
end

% Display result
if iter == n
    fprintf('\nMaximum iterations reached. Method may not have converged.\n');
else
    fprintf('\nRoot is: %f\n', x1);
end
