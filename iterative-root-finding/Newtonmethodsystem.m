function [x,iterationCount,e] = Newtonmethodsystem(F, Fprime, x0, tol)
% NEWTON'S METHOD FOR NON LINEAR SYSTEMS
% Basheq Tarifi (1696842)
%
% This function approximates the solution to a non linear equation or a non
% linear system of equations given by the function handle F, by using
% Newton's Method for non linear equations or systems. F is either a single
% function or a vector of functions. Fprime is also a function handle and
% represents either the derivative of the equation F or the Jacobian of the
% system given by F. The function uses an initial guess given by x0 and
% approximates the solution to within the tolerance specified by tol. The
% function returns the root, the number of iterations as well as the vector
% of error values (e) at each iteration.

n = size(x0,1);

N = 100;
xk = zeros(n,1);
k = 0;
e = [];
stopCriteria = 0;
while (stopCriteria > tol)||(k == 0)

    J = Fprime(x0);
    
    % Using Gaussian Elimination
    % s = gaussElimination(J,-1*F(x0));
    % xk = s + x0;
    
    xk = x0 - J\F(x0);
    
    stopCriteria = max(abs(xk - x0)./abs(xk));
    e(k+1) = stopCriteria;
   
    k = k +1;
    if k > N
        error("There have been more than " + N + ...
        " iterations and the sequence appears to diverge.")
    end
    x0 = xk;
end

x = x0;
iterationCount = k-1;
e = e(2:end);
end

