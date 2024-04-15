function U = homogeneousBoundaryFun(N,M,x,t,epsilon)
U = zeros(M+1,N+1); 
U(1,:) = -sin(pi*x);     % eta - initial condition
U(:,1) = 0;    % gL - left boundary condition
U(:,end) = 0;   % gR - right boundary condition