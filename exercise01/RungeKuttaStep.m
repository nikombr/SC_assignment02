function [tnew, ynew, error] = RungeKuttaStep(fun,t,y,h,varargin)
% Values from Butcher's table
b = [1/6 2/3 1/6]*h;
d = [1/12 -1/6 1/12]*h;
c = [0 1/2 1]*h;
a21 = 1/2*h;
a31 = -1*h;
a32 = 2*h;
% Stage 1
t1 = t;
xi1 = y;
f1 = feval(fun,t1,xi1,varargin{:});
% Stage 2
t2 = t + c(2)*h;
xi2 = y + a21*f1;
f2 = feval(fun,t2,xi2,varargin{:});
% Stage 3
t3 = t + c(3)*h;
xi3 = y + a31*f1 + a32*f2;
f3 = feval(fun,t3,xi3,varargin{:});
% Do update
tnew = t + h;
ynew = y + b(1)*f1 + b(2)*f2 + b(3)*f3;
% Compute error with respect to order 2 method
error = norm(d(1)*f1 + d(2)*f2 + d(3)*f3);