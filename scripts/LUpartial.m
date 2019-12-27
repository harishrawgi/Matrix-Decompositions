%% Function description
% Function calculates the determinant and inverse of a matrix using partial
% pivoting. It also returns the LU decomposition and permutation matrix.
% Inputs: A (a matrix)
% Outputs: det (the determinant of the matrix)
%          A_inverse (the inverse of the matrix)

%% Function code
function [det, A_inverse] = LUpartial(A)

% fetch the size of A
n = size(A,1);

% variable to keep track of number of permutations, used while computing
% the determinant of the matrix
num_permutations = 0;

% The permutation matrix
P = eye(n);

% the LU decomposition matrices
L = zeros(n,n);
U = A;

%% LU decomposition of A using partial pivoting
for j=1:n
    
    % set the pivot index initially to current iterate
    pivot = j;
    
    % find the correct pivot by taking the max
    for i=j:n
        if abs(U(i,j)) >= abs(U(pivot,j))
            pivot= i;
        end
    end
    
    % error if the matrix is singular
    if U(pivot,j) == 0
        error('Matrix is singular');
    end
    
    % if we got a pivot which requires permutation then increment the count
    if ~(pivot==j)
        num_permutations = num_permutations + 1;
    end

    % Perform the permutation on L
    t=L(j,1:j-1);
    L(j,1:j-1)=L(pivot,1:j-1);
    L(pivot,1:j-1)=t;
    
    % Perform the permutation on U
    t=U(j,j:end);
    U(j,j:end)=U(pivot,j:end);
    U(pivot,j:end)=t;
    
    % Also perform the permutation on P, this is to generate the
    % permutation matrix
    t = P(pivot,:);
    P(pivot,:)  = P(j,:);
    P(j,:) = t;

    
    % Perform the reduction
    L(j,j)=1;
    for i=(1+j):size(U,1)
       c= U(i,j)/U(j,j);
       U(i,j:n)=U(i,j:n)-U(j,j:n)*c;
       L(i,j)=c;
    end
    
end

%% Determinant calculation

% calculate the magnitude of the determinant
det = 1;
for i=1:n
    det = det*U(i,i);
end

% give the sign acc to number of permutations performed
det = det*(-1)^num_permutations;

%% Inverse calculation using LU decomposition

% we'll need the identity matrix in inverse calculation
I = eye(n);

% initialize the A_inverse matrix
A_inverse = zeros(n);

% Use the LU decomposition to calculate the inverse
for i=1:n
    
    % get the b vector for this iteration
    b = I(:,i);
    
    % Solve Ly = Pb
    % forwardSub is a function that does forward substitution
    y = forwardSub(L,P*b);
    
    % Solve Ux = y
    % x is the column of A_inverse for this iteration
    % backwardSub is a function that does backward substitution
    A_inverse(:,i) = backwardSub(U,y);
end

end