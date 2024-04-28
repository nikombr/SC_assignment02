function U = boundaryFun(x,t,epsilon)
U        = zeros(length(t),length(x));  % Allocate space
U(:,1)   = tanhFun(-1,t,epsilon);       % gL - left boundary condition
U(:,end) = tanhFun(1,t,epsilon);        % gR - right boundary condition
U(1,:)   = tanhFun(x,0,epsilon);        % eta - initial condition