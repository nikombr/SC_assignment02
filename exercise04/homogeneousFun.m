function homogeneousFun(N,M,x,t)


U = zeros(M+1,N+1); 

U(1,:) = -sin(pi*x);     % eta - initial condition