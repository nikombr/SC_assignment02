function y = computeDerivative(u,x,alpha,beta,N,xbar)
k = 1;
h = 2/N;
mid = N/2+1;
a = -alpha*h;               % First displacement point in stencil
b = beta*h;                 % Last displacement point in stencil
n = alpha + beta + 1;       % Grid size
x = x((mid-alpha):(mid+beta)); % Get grid points
x
c = fdcoeffV(k,xbar,x);     % Get coefficients
us = u((mid-alpha):(mid+beta));                  % Get values for stencil
us
y = us*c';                  % Get derivative in point xbar
