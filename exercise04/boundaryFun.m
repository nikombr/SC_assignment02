function U = boundaryFun(x,t,epsilon)
U = zeros(length(t),length(x)); 
U(:,1) = fun(-1,t,epsilon);    % gL - left boundary condition
U(:,end) = fun(1,t,epsilon);   % gR - right boundary condition
U(1,:) = fun(x,0,epsilon);     % eta - initial condition