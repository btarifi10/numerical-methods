function x = LUsolve(A,b)
% SOLVING LINEAR SYSTEMS USING LU DECOMPOSITION
% Written by Basheq Tarifi (1696842)
% 
% This function takes as inputs an n*n matrix A and n*1 column vector b,
% and returns the n*1 column vector x.
%
% It solves the linear system Ax = b by LU Decomposition, using the
% function "LU.m" written by Basheq Tarifi to decompose the matrix
% A into matrices L and U such that PA = LU. It then solves LUx = Pb by
% setting y = Ux, solving for y in Ly = Pb and then finds x by solving
% Ux = y.


if size(A,1) ~= size(A,2)
    error('Please input a square matrix for A');
end

[P,L,U] = LU(A);

n = size(A,1);
b = P*b;

%Forward substitution to solve Ly = b
y = zeros(n,1);
y(1,1) = b(1,1)/L(1,1);
for k = 2:n
    y(k,1) = b(k,1) - (L(k,1:k-1)*y(1:k-1,1));
end

%Backward substitution to solve Ux = y
x = zeros(n,1);
x(n,1) = y(n,1)/U(n,n);
for k = n-1:-1:1
    x(k,1) = (y(k,1) - (U(k,k+1:n)*x(k+1:n,1)))/U(k,k);
end

end

