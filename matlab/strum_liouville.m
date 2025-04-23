function [e_value, e_function, non_zero] = strum_liouville(L)
    syms y(x) lambda n L
    
    disp('Solving for lambda > 0...');
    assume(lambda > 0);
    solution = dsolve(diff(y,2) + lambda * y == 0);
    e_function = solution;
    disp('General solution:');
    disp(e_function);

    diff_sol = diff(solution, x);
    vals = solve(subs(diff_sol, x, 0) == 0, subs(diff_sol, x, L) == 0);
    % Assuming vals is a struct with fields (like vals.n or vals.lambda)
non_zero = simplify(struct2array(vals)); 

    
    disp('Non-zero values in the function:');
    disp(non_zero);
    
    e_value = (n * pi / L)^2;
    
    disp('Solving for lambda = 0...');
    solution_zero = dsolve(diff(y,2) == 0);
    disp('General solution for lambda = 0:');
    disp(solution_zero);

    disp('Solving for lambda < 0...');
    assume(lambda < 0);
    solution_neg = dsolve(diff(y,2) + lambda * y == 0);
    disp('General solution for lambda < 0:');
    disp(solution_neg);
    diff_sol_neg = diff(solution_neg, x);
    
    vals_neg = solve(subs(diff_sol_neg, x, 0) == 0, subs(diff_sol_neg, x, L) == 0);
    disp('Non-trivial solutions for lambda < 0:');
    disp(vals_neg);
end