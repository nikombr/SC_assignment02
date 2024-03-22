function [tnew, ynew, error] = RungeKuttaStep(fun,t,y,h,varargin)

b = [1/6 2/3 1/6]*h;
d = [1/12 -1/6 1/12]*h;
c = [0 1/2 1]*h;
a21 = 1/2*h;
a31 = -1*h;
a32 = 2*h;

t1 = t;
xi1 = y;
f1 = feval(fun,t1,xi1,varargin{:});

t2 = t + c(2)*h;
xi2 = y + a21*f1;
f2 = feval(fun,t2,xi2,varargin{:});

t3 = t + c(3)*h;
xi3 = y + a31*f1 + a32*f2;
f3 = feval(fun,t3,xi3,varargin{:});


tnew = t + h;
ynew = y + b(1)*f1 + b(2)*f2 + b(3)*f3;
error = norm(d(1)*f1 + d(2)*f2 + d(3)*f3);
