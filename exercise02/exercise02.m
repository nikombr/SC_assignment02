
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% Test of stability condition

clear; close all; clc;

% Parameters
N = 200;
tmax = 0.5;
epsilon = 0.1;
k = ((2/N)^2/(2*epsilon)); % Largest value while preserving stability
M = ceil(tmax/k);

% Call function
U = ForwardTimeCentralSpace(N,M,tmax,epsilon);

% Plot solution
plotSolution(U,M,N,tmax,epsilon);

% Save plot of solution
save_name = sprintf('../plots/exercise02/solution_plot_N_%d_M_%d.png',N,M);
exportgraphics(gcf,save_name,'Resolution',300);


% Show case of good stability
knew = k*0.8
M = ceil(tmax/knew);

% Call function
U = ForwardTimeCentralSpace(N,M,tmax,epsilon);

% Plot solution
plotSolution(U,M,N,tmax,epsilon);

% Save plot of solution
save_name = sprintf('../plots/exercise02/solution_plot_N_%d_M_%d.png',N,M);
exportgraphics(gcf,save_name,'Resolution',300);

% Show case of non-stability
knew = k*1.1
M = ceil(tmax/knew);

% Call function
U = ForwardTimeCentralSpace(N,M,tmax,epsilon);

% Plot solution
plotSolution(U,M,N,tmax,epsilon);

% Save plot of solution
save_name = sprintf('../plots/exercise02/solution_plot_N_%d_M_%d.png',N,M);
exportgraphics(gcf,save_name,'Resolution',300);

%% Solution in case of different sizes of grids

clear; close all; clc;

% Parameters
N = 120;
tmax = 0.5;
epsilon = 0.1;
k = ((2/N)^2/(2*epsilon)); % Largest value while preserving stability
M = ceil(tmax/k);

% Call function
U = ForwardTimeCentralSpace(N,M,tmax,epsilon);

% Plot solution
plotSolution(U,M,N,tmax,epsilon);

% Save plot of solution
save_name = sprintf('../plots/exercise02/solution_plot_N_%d_M_%d.png',N,M);
exportgraphics(gcf,save_name,'Resolution',300);


% Parameters
N = 60;
tmax = 0.5;
epsilon = 0.1;
k = ((2/N)^2/(2*epsilon)); % Largest value while preserving stability
M = ceil(tmax/k);

% Call function
U = ForwardTimeCentralSpace(N,M,tmax,epsilon);

% Plot solution
plotSolution(U,M,N,tmax,epsilon);

% Save plot of solution
save_name = sprintf('../plots/exercise02/solution_plot_N_%d_M_%d.png',N,M);
exportgraphics(gcf,save_name,'Resolution',300);

% Parameters
N = 30;
tmax = 0.5;
epsilon = 0.1;
k = ((2/N)^2/(2*epsilon)); % Largest value while preserving stability
M = ceil(tmax/k);

% Call function
U = ForwardTimeCentralSpace(N,M,tmax,epsilon);

% Plot solution
plotSolution(U,M,N,tmax,epsilon);

% Save plot of solution
save_name = sprintf('../plots/exercise02/solution_plot_N_%d_M_%d.png',N,M);
exportgraphics(gcf,save_name,'Resolution',300);
%% Convergence test

clear; close all; clc;

% Parameters
M = 900;
tmax = 0.5;
epsilon = 0.1;
N = round(linspace(10,200,100));

t = linspace(0,tmax,M+1);
h = 2./N;

% Allocate space
error = zeros(length(N),1);

for i = 1:length(N)
    
    x = linspace(-1,1,N(i)+1);
    [X, T] = meshgrid(x,t);
    utrue = fun(X,T);
    
    % Call function
    U = ForwardTimeCentralSpace(N(i),M,tmax,epsilon);

    error(i) = max(max(U-utrue));

    k  = tmax/M;
    
    if k > h(i)^2/(2*epsilon)
        disp("Stability not satisfied!")
    end

end

figure('Renderer', 'painters', 'Position', [400 400 600 300]);

loglog(h,error,'LineWidth',1.5,'DisplayName','$||U-u(x,t)||_{\infty}$')
hold on
loglog(h,h.^2,'--','LineWidth',1.5,'DisplayName','$\mathcal{O}(h^2)$')
grid on

legend('fontsize',15,'location','northwest')
xlabel('$h$','fontsize',15)
ylabel('Error','fontsize',15)

% Save plot
save_name = sprintf('../plots/exercise02/convergence_spatial_M_%d.png',M);
exportgraphics(gcf,save_name,'Resolution',300);

%% hmm spørg Allan

clear; close all; clc;

% Parameters
N = 200;
h = 2/N;
tmax = 0.5;
epsilon = 0.1;
k = linspace(10^(-4),((2/N)^2/(2*epsilon)),100);
k = k(2:end); 
M = ceil(tmax/k(1));

% True solution
x = linspace(-1,1,N+1);

error = zeros(length(k),1);

for i = 1:length(k)
    i
    M = ceil(tmax/k(i));
    t = linspace(0,tmax,M+1);
    [X, T] = meshgrid(x,t);
    utrue = fun(X,T);

    k(i) = tmax/M;
    
    % Call function
    U = ForwardTimeCentralSpace(N,M,tmax,epsilon);

    error(i) = max(max(U-utrue));

    if k(i) > h^2/(2*epsilon)
        disp("Stability not satisfied!")
    end

end



loglog(k,error)
hold on
loglog(k,k)