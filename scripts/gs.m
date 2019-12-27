%% Function Description
% The function performs and returns the QR decomposition of a matrix
% It uses the standard Gram-Schmidt procedure
% Input: Matrix A
% Output: Q and R (the QR decomposition of A

%% Function code
function [Q,R] = gs(A)

% get the size of the matrix
n = size(A,1);

% initialize the Q and R matrices
Q = zeros(n,n);
R = zeros(n,n);

% compute first column of Q directly
R(1,1) = norm(A(:,1));
Q(:,1) = A(:,1)/R(1,1);

% iteratively apply GS to calculate Q and R
for j=2:n
    
    Q(:,j) = A(:,j);
    
    for i=1:j-1
        R(i,j) = Q(:,i)'*A(:,j);
        Q(:,j) = Q(:,j) - R(i,j)*Q(:,i);
    end
   
    R(j,j) = norm(Q(:,j));
    Q(:,j) = Q(:,j)/R(j,j);
    
end