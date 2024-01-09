function o = diagDominant(A)
% DIAGONALLY DOMINANT CHECK
% Written by Basheq Tarifi (1696842)

% This function takes as inputs a square matrix A and determines whether or
% not it is diagonally dominant.

o = true;
A = abs(A);
n = size(A,1);
for i = 1:n
   tempSum = sum(A(i,1:i-1)) + sum(A(i,i+1:n));
   if A(i,i) < tempSum
       o = false;
       return
   end
end

end

