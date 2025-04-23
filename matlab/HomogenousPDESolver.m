function HomogenousPDESolver()
    clc;  close all;

    syms D D2 x y Z f(x) a b c m;
    
    % Define the homogeneous differential equation where D is derivative
    % w.r.t x and D2 is derivative w.r.t. y.
    lhs = (D^3 - 6*D^2 + 11*D - 6*D2^2)*Z;
    rhs = exp(5*x + 6*y);
    
    % Substituting D and D2 with m to solve the characteristic equation
    eq = subs(lhs, D, m);
    eq = subs(eq, D2, 1); % Assume D2 = 1 for simplicity
    disp('Characteristic Equation:');
    disp(eq);
    
    % Solve the characteristic equation for m
    value = solve(eq, m);
    disp('Roots of the characteristic equation:');
    disp(value);
    
    % Forming the complementary function (CF)
    if length(value) == 3
        if value(1) ~= value(2) && value(2) ~= value(3) && value(1) ~= value(3)
            CF = a * exp(value(1) * x) + b * exp(value(2) * x) + c * exp(value(3) * x);
            disp('Complementary Function (CF):');
            disp(CF);
        elseif value(1) == value(2) && value(2) == value(3)
            % If all roots are equal
            CF = (a + b*x + c*x^2) * exp(value(1) * x);
            disp('Complementary Function (CF) when all roots are equal:');
            disp(CF);
        else
            % When two roots are equal
            CF = (a + b*x) * exp(value(1) * x) + c * exp(value(2) * x);
            disp('Complementary Function (CF) when two roots are equal:');
            disp(CF);
        end
    end
    
    % Finding values of constants X & Y
    X = diff(log(rhs), x);
    Y = diff(log(rhs), y);
    disp('X:');
    disp(X);
    disp('Y:');
    disp(Y);

    % Express lhs in terms of Z
    lhs = lhs / Z;
    
    % Finding the particular integral (PI)
    eqn = subs(lhs, D, X);
    eqn = subs(eqn, D2, Y);
    PI = rhs / eqn;
    disp('Particular Integral (PI):');
    disp(PI);
    
    % Displaying the final solution
    answer = CF + PI;
    disp('Final Solution:');
    disp(answer);
end