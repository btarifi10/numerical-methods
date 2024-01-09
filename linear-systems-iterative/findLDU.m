function [L,D,U] = findLDU(A)
% FIND L, D AND U MATRICES
% Written by Basheq Tarifi (1696842)

% This function takes as inputs a square matrix A and determines the
% stricly lower matrix L, diagonal matrix D and strictly upper matrix U
% such that L+D+U = A.

n = size(A,1);
D = zeros(n);
L = A;
U = A;
for j = 1:n
    D(j,j) = A(j,j);
    U(j,1:j) = zeros(1,j);
    L(j,j:n) = zeros(1,n-j+1);
end

end

