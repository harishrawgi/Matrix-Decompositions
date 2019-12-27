%% Function description
% The function converts LU decomposition to QR decomposition
% Inputs: L and U (the LU decomposition of a matrix)
% Outputs: Q and R (the QR decomposition of a matrix)

%% Function code
function [Q,R] = lutoqr(L,U)

% get the QR decomposition of matrix L using gs (or mgs)
[Q_L, R_L] = mgs(L);

% the Q matrix is the same as it is already an orthogonal matrix
Q = Q_L;

% get the R matrix using R_L and U
R = R_L*U;

% print the computed Q and R
fprintf("\nThe computed Q and R are below respectively:\n");
disp(Q);
disp(R);

end