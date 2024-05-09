
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

clear; close all; clc;

delta = 0.06;

y0 = delta;
t0 = 0;
tend = 2/delta;
h0 = 1/tend;
reps = 10^(-4);
aeps = 10^(-8);

[T, Y, H, E, numStep] = RungeKutta(@fun,y0,t0,tend,h0,reps,aeps);

figure;
plot(T,Y,'.-')

figure;
plot(T,H)

figure;
plot(T(2:end),E)

%% Compute solution for varying deltas

clear; close all; clc;
reps = 10^(-6);
aeps = 10^(-8);
opts = odeset('RelTol',reps,'AbsTol',aeps);
b = 10;
a = -log(200)/10;
lis = linspace(0,b,9);

%deltas = linspace(0.02, 0.0001,9);
deltas = 1/50*exp(a*lis);
T_cell = cell(length(deltas),1);
Y_cell = cell(length(deltas),1);
T_cell_ode23 = cell(length(deltas),1);
Y_cell_ode23 = cell(length(deltas),1);
H_cell = cell(length(deltas),1);
E_cell = cell(length(deltas),1);
numSteps = [];

for i = 1:length(deltas)
    delta = deltas(i);
    y0 = delta;
    t0 = 0;
    tend = 2/delta;
    h0 = 1/tend;
    
    [T_cell{i}, Y_cell{i}, H_cell{i}, E_cell{i}, numStep] = RungeKutta(@fun,y0,t0,tend,h0,reps,aeps);
    numSteps = [numSteps numStep];

    [T_cell_ode23{i}, Y_cell_ode23{i}] = ode23(@fun,[t0 tend],delta,opts);

end

plot_varying_delta(deltas, T_cell, Y_cell, H_cell, E_cell, numSteps, T_cell_ode23, Y_cell_ode23,reps,aeps);

%% Compute solution for different tolerances

clear; close all; clc;

b = 10;
a = -log(200)/10;
lis = linspace(0,b,100);
deltas = 1/50*exp(a*lis);


%reps = [10^(-2) 10^(-3) 10^(-4) 10^(-5) 10^(-6)];
%aeps = [10^(-4) 10^(-5) 10^(-6) 10^(-7) 10^(-8)];
reps = [10^(-3) 10^(-4) 10^(-5) 10^(-6)];
aeps = [10^(-5) 10^(-6) 10^(-7) 10^(-8)];

numSteps = cell(length(reps),length(aeps));

%deltas = linspace(0.02, 0.0001,9);
for k = 1:length(deltas)
    delta = deltas(k);
    y0 = delta;
    t0 = 0;
    tend = 2/delta;
    h0 = 1/tend;
    for i = 1:length(reps)
        for j = 1:length(aeps)
            
            [T, Y, H, E, numStep] = RungeKutta(@fun,y0,t0,tend,h0,reps(i),aeps(j));
            numSteps{i,j} = [numSteps{i,j} numStep];
        end
    
    end
end

figure('Renderer', 'painters', 'Position', [400 400 1200 800]);
tiledlayout(length(reps),length(aeps),'TileSpacing','compact');

for i = 1:length(reps)
    for j = 1:length(aeps)
        nexttile;
        semilogx(deltas,numSteps{i,j});
        if i == 1
            title(sprintf('\\texttt{{aeps = %.4e}}',aeps(j)),'FontSize',12)
        end
        grid on
        if i == length(reps)
            xlabel('$\delta$','FontSize',12)
        end
        if j == 1
            ylabel({sprintf('\\texttt{{reps = %.4e}}',reps(i)),'Number of Steps'},'FontSize',12)
        else
            %ylabel('Number of Steps','FontSize',12)
        end
    end

end

exportgraphics(gcf,'../plots/exercise01/tolerance_test.png','Resolution',300);

%% Convergence plot


clear; close all; clc;

b = 10;
a = -log(200)/10;
lis = linspace(0,b,100);
deltas = 1/50*exp(a*lis);


%reps = [10^(-2) 10^(-3) 10^(-4) 10^(-5) 10^(-6)];
%aeps = [10^(-4) 10^(-5) 10^(-6) 10^(-7) 10^(-8)];
reps = [10^(-3) 10^(-4) 10^(-5) 10^(-6)];
aeps = [10^(-5) 10^(-6) 10^(-7) 10^(-8)];
expo = linspace(-10,-3,10);
reps = 10.^expo;
aeps = 10.^expo;

error = cell(length(reps),length(aeps));
error_estimation = cell(length(reps),length(aeps));


%deltas = linspace(0.02, 0.0001,9);
for k = 1:length(deltas)
    delta = deltas(k);
    y0 = delta;
    t0 = 0;
    tend = 2/delta;
    h0 = 1/tend;
    for i = 1:length(reps)
        for j = 1:length(aeps)
            
            [T, Y, H, E, numStep] = RungeKutta(@fun,y0,t0,tend,h0,reps(i),aeps(j));
            error_estimation(i,j) = norm(E,'inf');
            error = 
        end
    
    end
end