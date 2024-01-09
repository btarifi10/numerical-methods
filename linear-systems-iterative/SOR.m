function [x,iterationCount] = SOR(A,b,x_0,tol)
% SUCCESSIVE OVER RELAXATION METHOD FOR ITERATIVE SOLUTION OF A
% LINEAR SYSTEM
% Written by Basheq Tarifi (1696842)

% This function takes as inputs a square matrix A and column vector b, as
% well as an initial guess for x, x_0, and a tolerance. It then calculates
% a weighting factor w, used to assist in the approximation of x in the
% equation Ax = b, within a tolerance provided. This is the SOR Method.
% It also provides the number of iterations required to get to this
% iteration of x.

x = NaN;  %output NaN if error
iterationCount = NaN;

if size(A,1) ~= size(A,2)
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

LnU = L+U;

invD = zeros(n);
for j = 1:n
    invD(j,j) = 1/D(j,j);
end

B_j = -1*invD*LnU;


sr_B = max(abs(eig(B_j)));

if sr_B >= 1
    disp(['Error: Iteration sequence does not converge because spectral '...
        'radius is not less than 1.'])
    return
end

% calculate optimal value for w
w = 2/(1+sqrt(1-sr_B^(2)));

x_k = zeros(n,1);
k = 0;
% begin iterations
while (1) % stopping condition is checked at the end
    for i = 1:n
        x_k(i,1) = x_0(i,1)+(w/A(i,i))*(b(i,1) - ...
          (A(i,1:i-1)*x_k(1:i-1,1) +  A(i,i:n)*x_0(i:n,1)));
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
