%% Function description
% The function performs forward substitution to solve Ly = b
% Inputs: L (a lower triangular matrix), b (a vector)
% Outputs: y (the solution vector)

%% Function code
function y = forwardSub(L,b)

% fetch the size of vector b
n = length(b);

% initialize y vector
y = zeros(n,1);

% Solve for y using forward substitution
y(1) = b(1)/L(1,1);
for j=2:n
    y(j)=(b(j) -sum(L(j,1:j-1)'.*y(1:j-1)))/L(j,j);
end