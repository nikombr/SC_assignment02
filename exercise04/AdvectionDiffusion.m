function [U,x,t] = AdvectionDiffusion(fun,N,M,tmax,epsilon,mesh,Nfine)

if nargin < 6 || strcmp(mesh,"uniform") % Uniform mesh
    k = tmax/M;
    % Initialize
    x = linspace(-1,1,N+1);
    t = linspace(0,tmax,M+1);
    h = 2/N;
    frac1 = epsilon * k/h^2;
    frac2 = k/h;
    % Compute conditions
    U = fun(x,t,epsilon);
    % Iterate through time 
    % old
%     for j = 1:length(t)-1
%         U(j+1,2:(end-1)) = U(j,2:(end-1)) + frac1 * (U(j,1:(end-2)) - 2 * U(j,2:(end-1)) + U(j,3:end)) ...
%             - frac2*U(j,2:(end-1)).*(U(j,2:(end-1))-U(j,1:(end-2)));
%     end

    % maybe smart


else % Other stuff

    if strcmp(mesh,"nonuniform")
        
        


    end




%     t = linspace(0,tmax,M+1);
%     k = tmax/M;
%     x1 = linspace(-1,-0.1,N);
%     x2 = linspace(-0.1,0.1,Nfine);
%     x3 = linspace(0.1,1,N);
%     x = [x1(1:(end-1)) x2 x3(2:end)];
%     h = 0.9/N;
%     hfine = 0.2/Nfine;
%     frac1 = epsilon * k/h^2;
%     frac2 = k/h;
%     frac1fine = epsilon * k/hfine^2;
%     frac2fine = k/hfine;
%     % Compute conditions
%     U = fun(x,t,epsilon);
%     % Iterate through time 
%     for j = 1:length(t)-1
%         U(j+1,2:(N-1)) = U(j,2:(N-1)) + frac1 * (U(j,1:(N-2)) - 2 * U(j,2:(N-1)) + U(j,3:N)) ...
%             - frac2*U(j,2:(N-1)).*(U(j,2:(N-1))-U(j,1:(N-2)));
%         U(j+1,N) = U(j,N) + frac1 * (U(j,N-1) - U(j,N)) + frac1fine*(U(j,N+1) - U(j,N)) ...
%             - frac2*U(j,N).*(U(j,N)-U(j,N-1));
%         U(j+1,(N+1):(Nfine+N-1)) = U(j,(N+1):(Nfine+N-1)) + frac1fine * (U(j,N:(Nfine+N-2)) - 2 * U(j,(N+1):(Nfine+N-1)) + U(j,(N+2):(Nfine+N))) ...
%             - frac2fine*U(j,(N+1):(Nfine+N-1)).*(U(j,(N+1):(Nfine+N-1))-U(j,N:(Nfine+N-2)));
%         U(j+1,Nfine+N-1) = U(j,Nfine+N-1) + ...
%                 frac1fine * (U(j,Nfine+N-2) - U(j,Nfine+N-1)) + ...
%                 frac1*(U(j,Nfine+N)- U(j,Nfine+N-1)) ...
%                 - frac2fine*U(j,Nfine+N-1).*(U(j,Nfine+N-1)-U(j,Nfine+N-2));
%         U(j+1,(Nfine+N):(end-1)) = U(j,(Nfine+N):(end-1)) + frac1 * (U(j,(Nfine+N-1):(end-2)) - 2 * U(j,(Nfine+N):(end-1)) + U(j,(Nfine+N+1):end)) ...
%             - frac2*U(j,(Nfine+N):(end-1)).*(U(j,(Nfine+N):(end-1))-U(j,(Nfine+N-1):(end-2)));
%     end

end