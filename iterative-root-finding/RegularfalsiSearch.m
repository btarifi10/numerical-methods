function [root,iterationCount,e] = RegularfalsiSearch(f,tol,I_0)
% REGULAR FALSI SEARCH 
% Basheq Tarifi (1696842)
%
% This function approximates the solution to a non linear equation given by
% the function handle f, according to the iterative Regular Falsi (False
% Position) Method. The solution lies within the interval given by I_0 and
% has a relative error less than the tolerance specified. The function
% returns the root, the number of iterations as well as the vector of error
% values (e) at each iteration.

a_0 = I_0(1);
b_0 = I_0(2);
c_0 = (a_0*f(b_0) - b_0*f(a_0))/(f(b_0) - f(a_0));
e = [];
k = 1;
N = 100;
stopCriteria = 0;
while (stopCriteria > tol)||(k==1)
    if (f(c_0) == 0)
        c_k = c_0;
        k = k+1;
        e(k) = 0;
        break
    end
 
    if (k == 1)
        if f(c_0)*f(a_0) < 0
            toLeft = true;
            b_k = c_0;
            a_k = a_0;
        elseif f(c_0)*f(b_0) < 0
            toLeft = false;
            a_k = c_0;
            b_k = b_0;
        end
    end
    
    if toLeft
        b_k = c_0;
    else
        a_k = c_0;
    end
    
    c_k = (a_k*f(b_k) - b_k*f(a_k))/(f(b_k) - f(a_k));
    
    stopCriteria = abs(c_k - c_0)/abs(c_k);
    e(k) = stopCriteria;
    
    k = k+1;
    if k > N
        error("There have been more than " + string(N) + ...
        " iterations and the sequence appears to diverge.")
    end
    c_0 = c_k;
end

root = c_0;
iterationCount = k-1;
end

