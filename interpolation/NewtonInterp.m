function [co,T] = NewtonInterp(x,y)
% NEWTON POLYNOMIAL INTERPOLATION
% Basheq Tarifi (1696842)
%
% This function takes x and y column vectors, where y is a function value
% at x. The output is the Divided Difference table of Newton's
% interpolating polynomial of degree n-1, where n is the length of the
% column vectors x and y.

n = length(x);

T = zeros(n,n-1);
D = [y T];
for j = 2:n
    for i = 1:n-j+1
       D(i,j) = (D(i+1,j-1) - D(i,j-1))/(x(j+i-1)-x(i));
    end
end

T = D;
co = T(1,:);

end

