
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

%%

clear; close all; clc;
reps = 10^(-4);
aeps = 10^(-8);
opts = odeset('RelTol',reps);
b = 10;
a = -log(200)/10;
lis = linspace(0,10,9);

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

plot_varying_delta(deltas, T_cell, Y_cell, H_cell, E_cell, numSteps, T_cell_ode23, Y_cell_ode23);

lis = linspace(0,10,100);
numSteps = [];

%deltas = linspace(0.02, 0.0001,9);
deltas = 1/50*exp(a*lis);
for i = 1:length(deltas)
    delta = deltas(i);
    y0 = delta;
    t0 = 0;
    tend = 2/delta;
    h0 = 1/tend;
    
    [T_cell{i}, Y_cell{i}, H_cell{i}, E_cell{i}, numStep] = RungeKutta(@fun,y0,t0,tend,h0,reps,aeps);
    numSteps = [numSteps numStep];

end

figure;
semilogx(deltas,numSteps)

exportgraphics(gcf,'../plots/exercise01/number_of_steps.png','Resolution',300);
