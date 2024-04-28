
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

clear; close all; clc;

N = 200;
tmax = 1;
a = 0.5;
k = ((2/N)/(a)); % Largest value while preserving stability
M = ceil(tmax/k);
M = M*100

U = ForwardTimeBackwardSpace(N,M,tmax,a);

plotSolution(U,M,N,tmax,a);

%% Test of stability condition

clear; close all; clc;

% Parameters
N = 1000;
tmax = 2;
a = 0.5;
h = 2/N;
k = h/a/1.1; % Largest value while preserving stability
M = ceil(tmax/k);

% Call function
U = ForwardTimeBackwardSpace(N,M,tmax,a);

% Plot solution
plotSolution(U,M,N,tmax,a);

% Save plot of solution
save_name = sprintf('../plots/exercise03/solution_plot_N_%d_M_%d.png',N,M);
exportgraphics(gcf,save_name,'Resolution',300);


% Show case of good stability
knew = k*0.8
M = ceil(tmax/knew);

% Call function
U = ForwardTimeBackwardSpace(N,M,tmax,a);

% Plot solution
plotSolution(U,M,N,tmax,a);

% Save plot of solution
save_name = sprintf('../plots/exercise03/solution_plot_N_%d_M_%d.png',N,M);
exportgraphics(gcf,save_name,'Resolution',300);

% Show case of non-stability
knew = k*1.7
M = ceil(tmax/knew);

% Call function
U = ForwardTimeBackwardSpace(N,M,tmax,a);

% Plot solution
plotSolution(U,M,N,tmax,a);

% Save plot of solution
save_name = sprintf('../plots/exercise03/solution_plot_N_%d_M_%d.png',N,M);
exportgraphics(gcf,save_name,'Resolution',300);

%% Solution in case of different sizes of grids

clear; close all; clc;

% Parameters
N = 120;
tmax = 2;
a = 0.5;
h = 2/N;
k = h/a/1.1; % Largest value while preserving stability
M = ceil(tmax/k);

% Call function
U = ForwardTimeBackwardSpace(N,M,tmax,a);

% Plot solution
plotSolution(U,M,N,tmax,a);

% Save plot of solution
save_name = sprintf('../plots/exercise03/solution_plot_N_%d_M_%d.png',N,M);
exportgraphics(gcf,save_name,'Resolution',300);


% Parameters
N = 60;
k = h/a/1.1; % Largest value while preserving stability
M = ceil(tmax/k);

% Call function
U = ForwardTimeBackwardSpace(N,M,tmax,a);

% Plot solution
plotSolution(U,M,N,tmax,a);

% Save plot of solution
save_name = sprintf('../plots/exercise03/solution_plot_N_%d_M_%d.png',N,M);
exportgraphics(gcf,save_name,'Resolution',300);

% Parameters
N = 30;
k = h/a/1.1; % Largest value while preserving stability
M = ceil(tmax/k);

% Call function
U = ForwardTimeBackwardSpace(N,M,tmax,a);

% Plot solution
plotSolution(U,M,N,tmax,a);

% Save plot of solution
save_name = sprintf('../plots/exercise03/solution_plot_N_%d_M_%d.png',N,M);
exportgraphics(gcf,save_name,'Resolution',300);

%% Convergence test

clear; close all; clc;

% Parameters
M = 1200;
tmax = 2;
a = 0.5;
N = round(linspace(10,1200,100));

t = linspace(0,tmax,M+1);
h = 2./N;

% Allocate space
error = zeros(length(N),1);

for i = 1:length(N)
    i
    x = linspace(-1,1,N(i)+1);
    [X, T] = meshgrid(x,t);
    utrue = fun(X,T);
    
    % Call function
    U = ForwardTimeBackwardSpace(N(i),M,tmax,a);

    error(i) = max(max(U-utrue));

    k  = tmax/M;
    
    if k > h(i)/a
        disp("Stability not satisfied!")
    end

end

figure('Renderer', 'painters', 'Position', [400 400 600 300]);

loglog(h,error,'LineWidth',1.5,'DisplayName','$||U-u(x,t)||_{\infty}$')
hold on
loglog(h,h,'--','LineWidth',1.5,'DisplayName','$\mathcal{O}(h)$')
grid on

legend('fontsize',15,'location','northwest')
xlabel('$h$','fontsize',15)
ylabel('Error','fontsize',15)
xlim([min(h),max(h)])

% Save plot
save_name = sprintf('../plots/exercise03/convergence_spatial_M_%d.png',M);
exportgraphics(gcf,save_name,'Resolution',300);


%% hmm spørg Allan

clear; close all; clc;

% Parameters
N = 200;
h = 2/N;
tmax = 2;
a = 0.5;
k = linspace(10^(-4),h/a,100);
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
    U = ForwardTimeBackwardSpace(N,M,tmax,a);

    error(i) = max(max(U-utrue));

    if k(i) > h/a
        disp("Stability not satisfied!")
    end

end



loglog(k,error)
hold on
loglog(k,k)