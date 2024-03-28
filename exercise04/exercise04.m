
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

clear; close all; clc;

N = 200;
tmax = 0.5;
epsilon = 0.1;
k = ((2/N)/(epsilon)); % Largest value while preserving stability
k = 0.001;
k = 0.00001;
M = ceil(tmax/k);
M = 20;

U = AdvectionDiffusion(@boundaryFun,N,M,k,epsilon,epsilon);

plotSolution(@fun,U,M,N,k,epsilon);

%%


clear; close all; clc;

N = 200;
tmax = 0.5;
epsilon = 0.01/pi;
k = ((2/N)/(epsilon)); % Largest value while preserving stability
k = 0.001;
k = 0.00001;
M = ceil(tmax/k);
M = 20;

U = AdvectionDiffusion(@boundaryFun,N,M,k,epsilon);