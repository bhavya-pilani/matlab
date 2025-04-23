% Clear the screen
clc;

% Define x as a symbolic variable
syms x;

% Input Section
y = input('Enter the nonlinear equation: ');
a = input('Enter the first guess : ');
b = input('Enter the second guess : ');
e = input('Enter the tolerable error : ');
n = input('Enter the maximum number of iterations: ');

% Evaluate functional values at initial guesses
fa = eval(subs(y, x, a));
fb = eval(subs(y, x, b));

% Initialize
iter = 0;
fprintf('Iter\t a\t\t b\t\t c\t\t f(c)\n');

% Perform Secant Method
while iter < n
    if fb - fa == 0
        disp('Division by zero detected in the secant formula.');
        break;
    end

    c = b - fb * (b - a) / (fb - fa);
    fc = eval(subs(y, x, c));
    
    fprintf('%d\t%f\t%f\t%f\t%f\n', iter, a, b, c, fc);
    
    if abs(fc) < e
        break;
    end
    
    % Update points
    a = b;
    fa = fb;
    b = c;
    fb = fc;
    
    iter = iter + 1;
end

% Display result
if iter == n
    fprintf('\nMaximum iterations reached. Method may not have converged.\n');
else
    fprintf('\nRoot is: %f\n', c);
end
