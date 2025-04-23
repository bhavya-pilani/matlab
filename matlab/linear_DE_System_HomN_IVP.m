function [sols] = linear_DE_System_HomN_IVP(~)
    syms t lambda
    % Input the system matrix, non-homogeneous term, and initial conditions
    A = input('Enter the coefficient matrix A: ');
    F = input('Enter the non-homogeneous term (vector or matrix) F(t): ');
    C0 = input('Enter the initial condition vector (e.g., [1; 0; 1]): ');
    n = length(A);
    
    % Solve for eigenvalues and eigenvectors
    [V, D] = eig(A);
    eigenvalues = diag(D);
    const = sym('c', [n, 1]); % Constants for homogeneous solution
    unique_eigenvalues = unique(eigenvalues);
    mults = histc(eigenvalues, unique_eigenvalues); % Count multiplicities
    generalized_V = sym(zeros(n, n)); % Store all eigenvectors and generalized eigenvectors
    i = 1;
    
    % Handle repeated eigenvalues
    if length(unique_eigenvalues) ~= length(eigenvalues)
        for idx = 1:length(unique_eigenvalues)
            lambda_i = unique_eigenvalues(idx);
            mult = mults(idx);
            ch_mat = A - lambda_i * eye(n);
            
            % Compute eigenvectors and generalized eigenvectors
            for j = 1:mult
                if j == 1
                    % Basic eigenvector
                    e_vector = V(:, i);
                else
                    % Generalized eigenvector
                    e_vector = (ch_mat^(j-1)) \ generalized_V(:, i + j - 2);
                end
                
                % Scale by t^(j-1)
                generalized_V(:, i + j - 1) = e_vector * t^(j-1);
            end
            i = i + mult;
        end
    else
        generalized_V = V; % No repeated eigenvalues
    end
    
    % Homogeneous solution
    hom_sols = sym(zeros(n, 1));
    for k = 1:n
        hom_sols = hom_sols + exp(eigenvalues(k) * t) * generalized_V(:, k) * const(k);
    end
    
    % Particular solution for non-homogeneous system
    % Variation of parameters
    W = generalized_V * diag(exp(eigenvalues * t)); % Fundamental matrix
    inv_W = simplify(inv(W)); % Inverse of fundamental matrix
    int_term = int(inv_W * F, t); % Integral of (W^-1 * F)
    part_sols = simplify(W * int_term); % Particular solution
    
    % Combine homogeneous and particular solutions
    sols = simplify(hom_sols + part_sols);
    
    % Apply initial conditions
    sols_at_t0 = subs(sols, t, 0); % Solution at t = 0
    eqns = sols_at_t0 == C0; % Equations for constants
    vals = solve(eqns, const); % Solve for constants
    
    % Substitute constants into the solution
    for k = 1:n
        sols = subs(sols, const(k), vals.(sprintf('c%d', k)));
    end
    
    % Simplify the final solution
    sols = simplify(sols);
end