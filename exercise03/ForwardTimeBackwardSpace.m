function U = ForwardTimeBackwardSpace(N,M,k,a)
% Initialize
U = zeros(M+1,N+1); 
x = linspace(-1,1,N+1);
h = 2/N;
frac = a*k/h;

A = diag(ones(N,1),-1)-diag(ones(N+1,1));
A(1,N) = 1;
A = frac*A;

% Compute conditions
U(1,:) = fun(x,0);     % eta - initial condition

% Iterate through time 
for j = 1:M
    U(j+1,:) = U(j,:)' + A*U(j,:)';
end