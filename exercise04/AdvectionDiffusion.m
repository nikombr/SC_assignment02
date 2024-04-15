function U = AdvectionDiffusion(fun,N,M,k,epsilon)

tmax = k*M
k=tmax/M
% Initialize
x = linspace(-1,1,N+1);
t = linspace(0,tmax,M+1);
%k = tmax/M;
k
x
t
h = 2/N;
epsilon
frac1 = epsilon * k/h^2;
frac2 = k/h;
% Compute conditions
U = fun(N,M,x,t,epsilon);
% Iterate through time 
for j = 1:length(t)-1
    U(j+1,2:(end-1)) = U(j,2:(end-1)) + frac1 * (U(j,1:(end-2)) - 2 * U(j,2:(end-1)) + U(j,3:end)) ...
        - frac2*U(j,2:(end-1)).*(U(j,2:(end-1))-U(j,1:(end-2)));
end