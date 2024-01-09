function [X,U,C] = gaussMultipleSystems(A,B)
% GAUSS ELIMINATION WITH MULTIPLE SYSTEMS
% Written by Basheq Tarifi (student number 1696842)
% 
% This function takes as inputs an n*n matrix A and a n*m column vector B.
% It returns the n*m column vector X by solving the equation AX = B via
% Gauss Elimination using partial pivoting. It also returns error 
% messages in the event of infinitely many or no solutions.

if size(A,1) ~= size(A,2)
    error('Please input a square matrix for A');
elseif size(A,1) ~= size(B,1)
    error('Sizes of A and b do not match');
else
    n = size(A,1);
    %Augmented Matrix
    AB = [A B];
    for j = 1:(n-1)
        % Perform partial pivoting
        [~, indexMax] = max(abs(AB(j:n,j)));
        if indexMax ~= 1
            tempRow = AB(j+indexMax-1,:);
            AB(j+indexMax-1,:) = AB(j,:);
            AB(j,:) = tempRow;
        end
        if AB(j:n,j) == zeros(n-j+1,1)
            error("Matrix is singular, no unique solution");
        end
        m = AB(j,j);
        % Perform row operations to reduce matrix to upper triangular form
        for i = (j+1):n
            %Row Manipulation
            AB(i,:) = AB(i,:) - (AB(i,j)/m)*AB(j,:);
        end
    end
    
    %Seperate augmented matrix into final form of A and b
    U = AB(:,1:n);
    C = AB(:,n+1:end);

    if U(n,:) == zeros(1,n)
        error('Matrix is singular, no unique solution.');
    end
    % Perform backward substitution to find X vector
    X = zeros(size(C));
 
    X(n,:) = C(n,:)/U(n,n);
    for k = n-1:-1:1
        X(k,:) = (C(k,:) - (U(k,k+1:n)*X(k+1:n,:)))/U(k,k);
    end
end
end


