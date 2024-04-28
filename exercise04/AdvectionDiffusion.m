function [U,x,t] = AdvectionDiffusion(fun,N,M,tmax,epsilon,mesh,Nfine)

if nargin < 6 || strcmp(mesh,"uniform") % Uniform mesh
    
    % Initialize
    x = linspace(-1,1,N+1);
    
    %h = 2/N;
    %frac1 = epsilon * k/h^2;
    %frac2 = k/h;
    % Compute conditions
    %U = fun(x,t,epsilon);
    % Iterate through time 
    % old
%     for j = 1:length(t)-1
%         U(j+1,2:(end-1)) = U(j,2:(end-1)) + frac1 * (U(j,1:(end-2)) - 2 * U(j,2:(end-1)) + U(j,3:end)) ...
%             - frac2*U(j,2:(end-1)).*(U(j,2:(end-1))-U(j,1:(end-2)));
%     end

    % maybe smart
    m = length(x);

    % Second derivative of the spatial variable
    vec = fdcoeffV(2, x(2), x((1):(3))).*ones(m,1);
    FTCS2 = spdiags(vec,-1:1,m,m);
    FTCS2 = FTCS2(2:(end-1),:);
    
    % First derivative of the spatial variable
    vec = fdcoeffV(1, x(2), x((1):(2))).*ones(m,1);
    FTBS = spdiags(vec,-1:0,m,m);
    FTBS = FTBS(2:(end-1),:);


elseif strcmp(mesh,"nonuniform") % Non-uniform mesh

        
    x1 = linspace(-1,-0.1,N);
    x1 = x1(1:(end-1));
    x2 = linspace(-0.1,0.1,Nfine);
    x3 = linspace(0.1,1,N);
    x3 = x3(2:end);
    %x   = linspace(-1,1,2*N+Nfine-2)
    x   = [x1 x2 x3];

    m1 = length(x1);
    m2 = length(x2);
    m3 = length(x3);
    m = m1 + m2 + m3;

    % Second derivative of the spatial variable
    FTCS2 = sparse(zeros(m-2,m));

    vec = fdcoeffV(2, x(2), x((1):(3))).*ones(m1+1,1);
    temp = spdiags(vec,-1:1,m1+1,m1+1);
    FTCS2(1:(m1-1),1:(m1+1)) = temp(2:(end-1),:);
    FTCS2((m1+m2):end,(m1+m2):end) = temp(2:(end-1),:);

    vec = fdcoeffV(2, x2(2), x2((1):(3))).*ones(m2,1);
    temp = spdiags(vec,-1:1,m2,m2);
    temp = temp(2:(end-1),:);
    FTCS2((m1+1):(m1+m2-2),(m1+1):(m1+m2)) = temp;

    FTCS2(m1,m1:(m1+2)) = fdcoeffV(2, x2(1), [x1(end) x2(1) x2(2)]);
    FTCS2(m1+m2-1,(m1+m2-1):(m1+m2+1)) = fdcoeffV(2, x2(end), [x2(end-1) x2(end) x3(1)]);
    
    % First derivative of the spatial variable
    FTBS = sparse(zeros(m-2,m));

    vec = fdcoeffV(1, x1(2), x1(1:2)).*ones(m1+2,1);
    temp = spdiags(vec,-1:0,m1+2,m1+2);
    temp = temp(2:(end-1),:);
    FTBS(1:m1,1:(m1+2)) = temp;
    FTBS((end-m3+2):end,(end-m3):end) = temp(2:end,2:end);

    vec = fdcoeffV(1, x2(2), x2(1:2)).*ones(m2,1);
    temp = spdiags(vec,-1:0,m2,m2);
    temp = temp(2:end,:);
    FTBS((m1+1):(m1+m2-1),(m1+1):(m1+m2)) = temp;

end

FTBS = full(FTBS);
FTCS2 = full(FTCS2);

% Time steps
k = tmax/M;
t = linspace(0,tmax,M+1);

% Compute conditions
U = fun(x,t,epsilon);

%m = length(x);


% 
% FTCS2 = zeros(m);
% 
% %temp = fdcoeffV(2, x(2), x(1:3));
% %A(1,1:2) = temp(2:3);
% for i=2:(m-1)
%     FTCS2(i,i-1:i+1) = fdcoeffV(2, x(i), x((i-1):(i+1)));
% end
% %temp = fdcoeffV(2, x(m-1), x((m-2):m));
% %A(m,(m-1):m) = temp(1:2);
% 
% FTCS2 = FTCS2(2:(end-1),:); %A(2:end-1,2:end-1); % FTCS
% FTCS2
% 
% FTBS = zeros(m);
% 
% %A(1,1:2) = fdcoeffV(1, x(2), x(1:2));
% for i=2:(m-1)
%     FTBS(i,i-1:i) = fdcoeffV(1, x(i), x((i-1):(i)));
% end
% 
% FTBS = FTBS(2:(end-1),:);
%FTBS = A %A(2:end-1,2:end-1); % FTBS


% Iterate through time 
for j = 1:M
    U(j+1,2:(end-1)) = U(j,2:(end-1))' + k*epsilon*FTCS2*U(j,:)' - U(j,2:(end-1))'.*(k*FTBS*U(j,:)');
end
% h = 2/N;
% frac1 = epsilon * k/h^2;
% frac2 = k/h;
% for j = 1:length(t)-1
%      U(j+1,2:(end-1)) = U(j,2:(end-1)) + frac1 * (U(j,1:(end-2)) - 2 * U(j,2:(end-1)) + U(j,3:end)) ...
%          - frac2*U(j,2:(end-1)).*(U(j,2:(end-1))-U(j,1:(end-2)));
%  end


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