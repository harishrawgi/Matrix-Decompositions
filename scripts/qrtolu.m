%% Function description
% The function converts QR decomposition to LU decomposition
% Inputs: Q and R (the QR decomposition of a matrix)
% Outputs: L and U (the LU decomposition of a matrix)

%% Function code
function [L ,U] = qrtolu(Q,R)

% compute the LU decomposition of Q

% get the LU decomposition of Q
[L_Q, U_Q] = lu_nopivot(Q);

% the L matrix same as L_Q as it is already a lower triangular matrix
L  = L_Q;

% use U_Q along with R to get U of final LU decomposition
U = U_Q*R;

% print the computed L and U
fprintf("\nThe computed L and U are below respectively:\n");
disp(L);
disp(U);

end