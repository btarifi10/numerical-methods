function errorMat = ivpmethods(f,x0,y0,h,xf,ytrue)
% ERROR MATRIX GENERATION FOR NUMERICAL IVP SOLVERS
% Basheq Tarifi (1696842)
%
% This function solves the ODE given by dy/dx using the functions Euler.m,
% Heun.m and rk4.m, and then computes the error between them and the true
% solution, ytrue. The output is an error matrix showing the absolute error
% at each point.

x = (x0:h:xf)';
n = length(x);
N = n-1;

[~,y_euler] = Euler(f,x0,y0,N,xf);
[~,y_heun] = Heun(f,x0,y0,N,xf);
[~,y_rk4] = rk4(f,x0,y0,N,xf);

errorMat = [x abs(ytrue-y_euler) ...
    abs(ytrue-y_heun) abs(ytrue-y_rk4)];
end

