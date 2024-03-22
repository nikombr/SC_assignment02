function U = ForwardTimeCentralSpace(N,M,k,epsilon)
% Initialize
U = zeros(M+1,N+1); 
x = linspace(-1,1,N+1);
t = linspace(0,k*M,M+1);
h = 2/N;
frac = epsilon * k/h^2;
calc = 1-2*frac;
% Compute conditions
U(:,1) = fun(-1,t);    % gL - left boundary condition
U(:,end) = fun(1,t);   % gR - right boundary condition
U(1,:) = fun(x,0);     % eta - initial condition
% Iterate through time 
for j = 1:length(t)-1
    U(j+1,2:(end-1)) = frac * (U(j,1:(end-2)) + U(j,3:end)) + calc * U(j,2:(end-1));
end