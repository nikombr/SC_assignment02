function U = AdvectionDiffusion(fun,N,M,k,epsilon,varargin)
% Initialize
x = linspace(-1,1,N+1);
t = linspace(0,k*M,M+1);
h = 2/N;
frac1 = epsilon * k/h^2;
frac2 = k/h;
% Compute conditions
U = fun(N,M,x,t,varargin{:});
% Iterate through time 
for j = 1:length(t)-1
    U(j+1,2:(end-1)) = U(j,2:(end-1)) + frac1 * (U(j,1:(end-2)) - 2 * U(j,2:(end-1)) + U(j,3:end)) ...
        - frac2*U(j,2:(end-1)).*(U(j,2:(end-1))-U(j,1:(end-2)));
end