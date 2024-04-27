function U = homogeneousBoundaryFun(x,t,epsilon)
U = zeros(length(t),length(x)); 
U(1,:) = -sin(pi*x);     % eta - initial condition
%U(:,1) = 0;    % gL - left boundary condition
%U(:,end) = 0;   % gR - right boundary condition