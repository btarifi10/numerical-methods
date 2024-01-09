function [x,y] = RK4system(f,xspan,y0,N)
% FOURTH ORDER RUNGE-KUTTA METHOD TO SOLVE SYSTEMS OF FIRST ORDER
% ORDINARY DIFFERENTIAL EQUATIONS
% Basheq Tarifi (1696842)
%
% This function computes the solution to the system of ordinary
% differential equation given by dy/dx = f, where y is a vector and f is a
% vectorised function handle with the same number of components as y.
% xspan specifies the inital and final values of x and y0 specifies the
% initial values of all the components of y, given as a row vector of
% length n. N is the number of subintervals to be used. The outputs are an
% x vector of length N+1 and a y matrix of size N+1 by n, where each column
% is the numerical solution to the corresponding y component. This
% Runge-Kutta method has an error of the order of h^4, where h=(xf-x0)/N.

x0 = xspan(1);
xf = xspan(2);
h = (xf-x0)/N;
x = (x0:h:xf)';

y = zeros(length(x),length(y0));
y(1,:) = y0;

for i = 1:N
    k1 = h*f(x(i),y(i,:))';
    k2 = h*f(x(i) + h/2,y(i,:) + k1/2)';
    k3 = h*f(x(i) + h/2,y(i,:) + k2/2)';
    k4 = h*f(x(i) + h,y(i,:) + k3)';
    y(i+1,:) = y(i,:) + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
end

end

