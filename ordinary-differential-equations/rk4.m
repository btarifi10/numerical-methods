function [x,y] = rk4(f,x0,y0,N,xf)
% FOURTH ORDER RUNGE-KUTTA METHOD TO SOLVE ORDINARY DIFFERENTIAL EQUATIONS
% Basheq Tarifi (1696842)
%
% This function computes the solution to the ordinary differential equation
% given by dy/dx = f, where f is a function handle. x0 and y0 specify the
% initial conditions, and xf specifies the stopping point of x. N is the
% number of subintervals between x0 and xf. The outputs are x and y vectors
% of length N+1, which specify all the x values and the approximations to
% the solution y. This Runge-Kutta method has an error of the order of h^4,
% where h is given by (xf-x0)/N.

h = (xf-x0)/N;
x = (x0:h:xf)';
y = zeros(size(x));
y(1) = y0;

for i = 1:N
    k1 = h*f(x(i),y(i));
    k2 = h*f(x(i) + h/2,y(i) + k1/2);
    k3 = h*f(x(i) + h/2,y(i) + k2/2);
    k4 = h*f(x(i) + h,y(i) + k3);
    y(i+1) = y(i) + (1/6)*(k1 + 2*k2 + 2*k3 + k4);
end

end

