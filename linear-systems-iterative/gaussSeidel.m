function [x,iterationCount] = gaussSeidel(A,b,x_0,tol)
% GAUSS-SEIDEL METHOD FOR ITERATIVE SOLUTION OF A LINEAR SYSTEM
% Written by Basheq Tarifi (1696842)

% This function takes as inputs a square matrix A and column vector b, as
% well as an initial guess for x, x_0, and a tolerance. It then calculates
% an approximation to x in the equation Ax = b, within a tolerance
% provided, via the Gauss-Seidel Method. It also provides the number of
% iterations required to get to this value.

x = NaN;  %output NaN if error
iterationCount = NaN;

if size(A,1) ~= size(A,2) % check size
    disp('Error: A needs to be a square matrix.')
    return
end
n = size(A,2);
if (size(x_0,1) ~= n) || (size(x_0,2) ~= 1)
    disp(['Error: Initial x vector x_0 needs to be a column vector with the'...
        ' same number of rows as A.'])
    return
end
if (size(b,1) ~= n) || (size(b,2) ~= 1)
    disp(['Error: b needs to be a column vector with the'...
        ' same number of rows as A.'])
    return
end

% check if convergence criteria is met
[L,D,U] = findLDU(A);

DnL = D+L;

invDnL = inv(DnL);

B_gs = -1*invDnL*U;

sr_B = max(abs(eig(B_gs)));

if sr_B >= 1
    disp(['Error: Iteration sequence does not converge because spectral '...
        'radius is not less than 1.'])
    return
end


x_k = zeros(n,1);
k = 0;
% begin iterations
while (1) % stopping condition is checked at the end
    for i = 1:n
        x_k(i,1) = (b(i,1) - (A(i,1:i-1)*x_k(1:i-1,1) + ...
            A(i,i+1:n)*x_0(i+1:n,1)))/A(i,i);
    end
    
    stopCriteria = max(abs(x_k - x_0))/max(abs(x_k));
    
    if (stopCriteria < tol) % if the error is below the tolerance, exit 
        x_k = x_0;
        break
    end
    % otherwise, continue and increase counter
    k = k+1;
    x_0 = x_k;
end

iterationCount = k;
x = x_k;
end

