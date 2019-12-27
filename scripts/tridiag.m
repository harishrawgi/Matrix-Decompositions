%% Function description
%Function to solve a tridiagonal system of equations
%Inputs: A (a tridiagonal matrix), b (a vector)
%Outputs: x (the solution to the system)

%% Function code
function x = tridiag(A,b)

% get the size of the matrix 
[m,n] = size(A);

%% Get the LU decomposition without pivoting
L = eye(n);
U = eye(n);

for i=1:n-1
    U(i,i+1) = A(i,i+1);
end

U(1,1) = A(1,1);

for k=2:n
    L(k,k-1) = A(k,k-1)/U(k-1,k-1);
    U(k,k) = A(k,k) - L(k,k-1)*A(k-1,k);
end

% print the L and U of the LU decomposition of A
fprintf("\nThe L and U matrices are below respectively:\n");
disp(L);
disp(U);

%% Solve the tridiagnonal system of equations using the computed LU decomposition

% y denotes the solution of Ly = b
% x denotes the solution to Ux = y
y = zeros(n,1);
x = zeros(n,1);

% solving for y
y(1) = b(1);
for k=2:n
    y(k) = b(k) - L(k,k-1)*y(k-1);
end

%solving for x
x(n) = y(n)/U(n,n);
for k=(n-1):-1:1
    x(k) = (y(k) - U(k,k+1)*x(k+1))/U(k,k);
end
    
fprintf("\nThe solution to the set of tridiagonal equations is x:\n");
disp(x);

end