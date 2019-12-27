%% Function Description
% The function performs LU decomposition of a matrix A
% Inputs: A (a matrix)
% Outputs: L and U (the LU decomposition of A)

%% Function code
function [L, U] = lu_nopivot(A)

% get the size of A
n = size(A, 1); 

% Initialize L
L = eye(n); 

% iteratively calculate L and U (A will become U)
for k = 1 : n
    
    % this means we can't continue without pivoting
    if (A(k,k) == 0) 
        Error('Pivoting is needed!');
    end

    L(k + 1 : n, k) = A(k + 1 : n, k) / A(k, k);
    
    % performing Gaussian Elimination
    for l = k + 1 : n
        A(l, :) = A(l, :) - L(l, k) * A(k, :);
    end
end

% set U as A
U = A;

end