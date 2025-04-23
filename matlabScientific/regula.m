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

% Check if initial values bracket the root
if fa * fb > 0
    disp('Given initial values do not bracket the root.');
else
    % Initialize
    iter = 0;
    c = a; % Initial assignment

    % Perform Regula Falsi Method
    while iter < n
        c_old = c;
        c = (a * fb - b * fa) / (fb - fa);
        fc = eval(subs(y, x, c));
        
        fprintf('%d\t%f\t%f\t%f\t%f\n', iter, a, b, c, fc);
        
        if abs(fc) < e || abs(c - c_old) < e
            break;
        end
        
        if fa * fc < 0
            b = c;
            fb = fc;
        else
            a = c;
            fa = fc;
        end
        
        iter = iter + 1;
    end

    % Display result
    if iter == n
        fprintf('\nMaximum iterations reached. Method may not have converged.\n');
    else
        fprintf('\nRoot is: %f\n', c);
    end
end
