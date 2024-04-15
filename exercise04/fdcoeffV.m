function c = fdcoeffV(k,xbar,x)
% function from the course book by Leveque

n = length(x);
A = ones(n,n);

xrow = (x(:)-xbar)'; % displacements as a row vector.

for i=2:n
    A(i,:) = (xrow .^ (i-1)) ./ factorial(i-1);
end

b = zeros(n,1); % b is right hand side,
b(k+1) = 1;     % so k’th derivative term remains
c = A\b;        % solve system for coefficients
c = c';         % row vector