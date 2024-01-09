function [x,A] = modifiedGaussElimination(A,b)
% MODIFIED GAUSS ELIMINATION
% Written by Basheq Tarifi (student number 1696842)
% 
% This function is a modified version of "gaussElimination.m". It takes as
% inputs an n*n matrix A and a n*1 column vector b, and returns the n*1
% column vector x as well as the modified version of A, by solving the
% equation Ax = b via Gauss Elimination, reducing the augmented matrix
% to lower right triangular form and using pivoting where necessary.
% It also returns error messages in the event of infinitely many
% or no solutions.

if size(A,1) ~= size(A,2)
    error('Please input a square matrix for A');
end
if size(A,1) ~= size(b,1)
    error('Sizes of A and b do not match');
end
if size(b,2) ~= 1
    error('The input b needs to be a column vector');
end

    n = size(A,1);
    %Augmented Matrix
    Ab = [A b];
    for j = 1:n-1
        % Perform pivoting if diagonal entry is zero
        if Ab(n-j+1,j) == 0
            nz = find(Ab(1:n-j,j));
            tempRow = Ab(nz(end),:);
            Ab(nz(end),:) = Ab(n-j+1,:);
            Ab(n-j+1,:) = tempRow;
        end
        m = Ab(n-j+1,j);
        
        % Perform row operations to reduce matrix to lower right triangular form
        for i = (n-j):-1:1
            %Row Manipulation
            Ab(i,:) = Ab(i,:) - (Ab(i,j)/m)*Ab(n-j+1,:);
            %Check for row of zeros (singularity)
            if isempty(intersect(Ab(i,1:n),zeros(1,n),'rows')) == 0
                if Ab(i,n+1) == 0
                  error(['Matrix A is singular, and Ax = b has ' ...
                  'infinitely many solutions.']);
                else
                  error(['Matrix A is singular, and Ax = b ' ...
                      'has no solution.']);
                end
            end
        end
    end
    
    %Seperate augmented matrix into final form of A and b
    A_f = Ab(:,1:n);
    b_f = Ab(:,n+1);
    
    % Perform forward substitution to find x vector
    x = zeros(n,1);
    x(n,1) = b_f(1,1)/A_f(1,n);
    for k = n-1:-1:1
        x(k,1) = (b_f(n-k+1,1) - (A_f(n-k+1,k+1:n)*x(k+1:n,1)))/A_f(n-k+1,k);
    end
    
    A = A_f;
end


