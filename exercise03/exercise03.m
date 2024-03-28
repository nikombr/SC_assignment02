
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

clear; close all; clc;

N = 2000;
tmax = 0.5;
a = 0.1;
k = ((2/N)/(a)); % Largest value while preserving stability
k = 0.001;
k = 0.00001;
M = ceil(tmax/k);
M = 200;

U = ForwardTimeBackwardSpace(N,M,k,a);

plotSolution(U,M,N,k);