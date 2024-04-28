
%% Parabolic

clear; close all; clc;

load 'results/parabolic.mat'


figure('Renderer', 'painters', 'Position', [400 400 800 300]);
tiledlayout(1,2,'TileSpacing','compact');

nexttile
plot(loss)
ylabel('loss','FontSize',15)
xlabel('epoch','FontSize',15)
axis square
nexttile
imagesc(x,t,U)
axis square


%% Hyperbolic

clear; close all; clc;

load 'results/hyperbolic.mat'


%% Advection

clear; close all; clc;

load 'results/advection.mat'
