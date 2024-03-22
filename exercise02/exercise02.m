
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

%% Test of method

clear; close all; clc;

% Parameters
N = 200;
tmax = 0.5;
epsilon = 0.1;
k = ((2/N)^2/(2*epsilon)); % Largest value while preserving stability
M = ceil(tmax/k);

% Call function
U = ForwardTimeCentralSpace(N,M,k,epsilon);

% Plot solution
plotSolution(U,M,N,k);

% Save plot of solution
save_name = sprintf('../plots/exercise02/solution_plot_N_%d_M_%d.png',N,M);
exportgraphics(gcf,save_name,'Resolution',300);

%% Convergence test

clear; close all; clc;

% Parameters
N = 100;
tmax = 0.5;
epsilon = 0.1;
k = linspace(0,((2/N)^2/(2*epsilon)),100);
k = k(2:end); % Remove 0
M = ceil(tmax/k(1));

% True solution
x = linspace(-1,1,N+1);
t = linspace(0,k(1)*M,M+1);
[X, T] = meshgrid(x,t);
utrue = fun(X,T);

error = zeros(length(k),1);

for i = 1:length(k)
    M = ceil(tmax/k(i));
    t = linspace(0,k(i)*M,M+1);
    [X, T] = meshgrid(x,t);
    utrue = fun(X,T);
    
    % Call function
    U = ForwardTimeCentralSpace(N,M,k(i),epsilon);

    error(i) = max(max(U-utrue));

end



plot(error)