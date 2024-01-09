function [P,L,U] = LU(A)
% LU DECOMPOSITION OF MATRIX A
% Written by Basheq Tarifi (1696842)
%
% This function takes an n*n matrix A and decomposes it into an lower
% triangular matrix L and an upper triangular matrix U, by using partial
% pivoting and Gaussian Elimination. It also records the row interchanges
% made and stores this in the permutation matrix P, such that P*A = L*U

if size(A,1) ~= size(A,2)
    error('Please input a square matrix for A.');
end

n = size(A,1);
P = eye(n);

L = eye(n);
U = zeros(n);

for j = 1:n
    [~, indexMax] = max(abs(A(j:n,j)));
        % Perform partial pivoting as well as swapping rows of P matrix
        % Swap rows of completed elements below the diagonal in L
        if indexMax ~= 1
            tempRow = A(j+indexMax-1,:);
                iTempRow = P(j+indexMax-1,:);
                    LTemp = L(j+indexMax-1,1:j-1);
                    
            A(j+indexMax-1,:) = A(j,:);
                P(j+indexMax-1,:) = P(j,:);
                    L(j+indexMax-1,1:j-1) = L(j,1:j-1);
                    
            A(j,:) = tempRow;
                P(j,:) = iTempRow;
                    L(j,1:j-1) = LTemp;
        end
        % Check if any zero column
        if A(j:n,j) == zeros(n-j+1,1)
            error("Matrix is singular, no unique solution.");
        end
        
        m = A(j,j);
        % Perform row operations to reduce matrix to upper triangular form
        for i = (j+1):n
            L(i,j) = A(i,j)/m;
            A(i,:) = A(i,:) - (A(i,j)/m)*A(j,:);
        end
        
        % Check for a row of zeros
        if A(n,:) == zeros(1,n)
            error(['Matrix is singular, there is no unique solution '...
            ' and it cannot be decomposed into L and U matrices.']);
            return
        end
        
        U = A;
end

end

