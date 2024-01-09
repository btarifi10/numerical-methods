function [x,y] = Heun(f,x0,y0,N,xf)
% HEUN'S METHOD TO SOLVE ORDINARY DIFFERENTIAL EQUATIONS
% Basheq Tarifi (1696842)
%
% Also known as the Trapezoidal Method or Improved Euler Method for
% solving first order ordinary differential equations numerically.
% This function computes the solution to the ordinary differential equation
% given by dy/dx = f, where f is a function handle. x0 and y0 specify the
% initial conditions, and xf specifies the stopping point of x. N is the
% number of subintervals between x0 and xf. The outputs are x and y vectors
% of length N+1, which specify all the x values and the approximations to
% the solution y. Heun's method has an error of the order of h^2, where h
% is given by (xf-x0)/N.

h = (xf-x0)/N;
x = (x0:h:xf)';
y = zeros(size(x));
y(1) = y0;

for i = 1:N
    k1 = f(x(i),y(i));
    k2 = f(x(i+1),y(i)+h*k1);
    y(i+1) = y(i) + (h/2)*(k1 + k2);
end

end

