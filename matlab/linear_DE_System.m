function [sols] = linear_DE_System(~)
    syms t lambda
    % Input the system matrix and the non-homogeneous term
    A = input('Enter the coefficient matrix A: ');
    F = input('Enter the non-homogeneous term (vector or matrix) F(t): ');
    n = length(A);
    
    % Solve for eigenvalues and eigenvectors
    [V, D] = eig(A);
    eigenvalues = diag(D);
    const = reshape(sym('c%d', [1 n]), n, 1);
    unique_eigenvalues = unique(eigenvalues);
    mults = histc(eigenvalues, unique_eigenvalues); % Count multiplicities
    sols = sym('x%d', [1 n]); 
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
                    e_vector = (ch_mat^(j-1))\ generalized_V(:, i + j - 2);
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
    generalized_V
   
    hom_sols = sym(zeros(n, 1));
    for k = 1:n
        hom_sols = hom_sols + exp(eigenvalues(k) * t) * (generalized_V(:, k).' * const);
    end
    hom_sols
    
    % Particular solution for non-homogeneous system
    % Variation of parameters
    W = generalized_V * exp(D * t) % Fundamental matrix
    inv_W = simplify(inv(W)) % Inverse of fundamental matrix
    int_term = int(inv_W * F, t) % Integral of (W^-1 * F)
    part_sols = simplify(W * int_term) % Particular solution
    
    % Combine homogeneous and particular solutions
    sols = simplify(hom_sols + part_sols);
end