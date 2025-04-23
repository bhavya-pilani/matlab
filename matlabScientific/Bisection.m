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
    % Initialize variables
    c = (a + b) / 2;
    fc = eval(subs(y, x, c));

    % Perform Bisection Method
    iter = 0;
    while abs(fc) > e && iter < n
        fprintf('%f\t%f\t%f\t%f\n', a, b, c, fc);
        if fa * fc < 0
            b = c;
        else
            a = c;
        end
        c = (a + b) / 2;
        fc = eval(subs(y, x, c));
        iter = iter + 1;
    end

    % Display the root or an error message if maximum iterations reached
    if iter == n
        fprintf('\nMaximum iterations reached. Method did not converge.\n');
    else
        fprintf('\nRoot is: %f\n', c);
    end
end
