function [U,x,t] = AdvectionDiffusion(fun,N,M,tmax,epsilon,type,a)

if strcmp(type,"uniform") || strcmp(type,"higher")  % Uniform mesh
    
    % Initialize
    x = linspace(-1,1,N+1);
    
    m = length(x);


elseif strcmp(type,"nonuniform") % Non-uniform mesh
        
    % Initialize
    xi = linspace(-1,1,N+1);
    
    x = (1-a)*xi.^3 + a*xi;

end

% Time steps
k = tmax/M;
t = linspace(0,tmax,M+1);

% Compute conditions
U = fun(x,t,epsilon);

if strcmp(type,"uniform")

    h = 2/N;
    frac1 = epsilon * k/h^2;
    frac2 = k/h;
    
    % Iterate through time 
    for j = 1:length(t)-1
        U(j+1,2:(end-1)) = U(j,2:(end-1)) + frac1 * (U(j,1:(end-2)) - 2 * U(j,2:(end-1)) + U(j,3:end)) ...
            - frac2*U(j,2:(end-1)).*(U(j,2:(end-1))-U(j,1:(end-2)));
    end


elseif strcmp(type,"higher") % Higher order method

    % Compute second derivative of spatial variable using central difference
    FTCS2 = zeros(N+1);
    FTCS2(2,1:3) = fdcoeffF(2, x(2), x(1:3));
    FTCS2(3,1:5) = fdcoeffF(2, x(3), x(1:5));
    for i=4:(N-2)
        FTCS2(i,i-3:i+3) = fdcoeffF(2, x(i), x((i-3):(i+3)));
    end
    FTCS2(N-1,N-3:N+1) = fdcoeffF(2, x(N-1), x((N-3):(N+1)));
    FTCS2(N,N-1:N+1) = fdcoeffF(2, x(N), x((N-1):(N+1)));
    FTCS2 = FTCS2(2:(end-1),:); % FTCS
    

    % Compute first derivative of time variable using forward difference
    FTBS = zeros(N+1);
    FTBS(2,1:2) = fdcoeffF(1, x(2), x((1):(2)));
    for i=3:(N-1)
        FTBS(i,i-2:i) = fdcoeffF(1, x(i), x((i-2):(i)));
    end
    FTBS(N,N-1:N) = fdcoeffF(1, x(N), x((N-1):(N)));
    FTBS = FTBS(2:(end-1),:);

    % Iterate through time 
    for j = 1:M
        U(j+1,2:(end-1)) = U(j,2:(end-1))' + k*epsilon*FTCS2*U(j,:)' - U(j,2:(end-1))'.*(k*FTBS*U(j,:)');
    end
    
else
    % Compute second derivative of spatial variable using central difference
    FTCS2 = zeros(N+1);
    for i=2:(N)
        FTCS2(i,i-1:i+1) = fdcoeffF(2, x(i), x((i-1):(i+1)));
    end
    FTCS2 = FTCS2(2:(end-1),:); % FTCS

    % Compute first derivative of time variable using forward difference
    FTBS = zeros(N+1);
    for i=2:(N)
        FTBS(i,i-1:i) = fdcoeffF(1, x(i), x((i-1):(i)));
    end
    FTBS = FTBS(2:(end-1),:);

    % Iterate through time 
    for j = 1:M
        U(j+1,2:(end-1)) = U(j,2:(end-1))' + k*epsilon*FTCS2*U(j,:)' - U(j,2:(end-1))'.*(k*FTBS*U(j,:)');
    end

end



end