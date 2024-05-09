
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
M = 4000;
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

figure('Renderer', 'painters', 'Position', [400 400 1000 300]);
tiledlayout(1,2,'TileSpacing','compact');
nexttile;

loglog(h,error,'LineWidth',1.5,'DisplayName','$||U-u(x,t)||_{\infty}$')
hold on
loglog(h,6*h,'--','LineWidth',1.5,'DisplayName','$\mathcal{O}(h)$')
grid on

legend('fontsize',15,'location','northwest')
xlabel('$h$','fontsize',15)
ylabel('Error','fontsize',15)
xlim([min(h),max(h)])

title(sprintf('$k=%.2f\\cdot 10^{-4}$',k*10^(4)),'FontSize',15)

% Parameters
tmax = 2;
a = 0.5;
N = round(linspace(10,1200,100));

h = 2./N;
k = h;
M = ceil(tmax./k);

% Allocate space
error = zeros(length(N),1);

for i = 1:length(N)
    i
    x = linspace(-1,1,N(i)+1);
    t = linspace(0,tmax,M(i)+1);
    [X, T] = meshgrid(x,t);
    utrue = fun(X,T);
    
    % Call function
    U = ForwardTimeBackwardSpace(N(i),M(i),tmax,a);

    error(i) = max(max(U-utrue));
    
    if k(i) > h(i)/a
        disp("Stability not satisfied!")
    end

end

nexttile;

loglog(h,error,'LineWidth',1.5,'DisplayName','$||U-u(x,t)||_{\infty}$')
hold on
loglog(h,4*h,'--','LineWidth',1.5,'DisplayName','$\mathcal{O}(h)$')
grid on

legend('fontsize',15,'location','northwest')
xlabel('$h$','fontsize',15)
ylabel('Error','fontsize',15)
xlim([min(h),max(h)])

title('$k=h$','FontSize',15)

% Save plot
save_name = '../plots/exercise03/convergence_combined.png';
exportgraphics(gcf,save_name,'Resolution',300);


%% Dissipation and Dispersion

clear; clc; close all;

dx = 1/100;
xi = 2*pi;
cr = 0.8;
a = 0.5;
dt = cr*dx/a;
x = -1:1/100:1;
T = 1/a;
tmax = 40*T;
n = 1:(tmax/dt);
t = n*dt;
M = length(t)-1;
N = length(x)-1;

theta = xi*dx;

absg = ((1+cr*(cos(theta)-1))^2+cr^2*sin(theta)^2)^(1/2);
img = -cr*sin(theta);
reg = 1 - cr + cr*cos(theta);
y = atan2(img,reg); % phase angle
dispersion = y/(-cr*theta);

dissipation = absg.^n;

figure('Renderer', 'painters', 'Position', [400 400 1200 350]);
tiledlayout(1,2,'TileSpacing','compact');
nexttile;
plot(t/T,dissipation,'LineWidth',1.5)
grid on
xlabel('$t/T$','fontsize',15)
ylabel('$|g(\xi)|^n$','fontsize',15)

[X, Tt] = meshgrid(x,t);
utrue = fun(X,Tt);

U = ForwardTimeBackwardSpace(N,M,tmax,a);

nexttile;
plot(x,U(end,:),'LineWidth',1.5,'DisplayName','Numerical Solution')
hold on
plot(x,utrue(end,:),'LineWidth',1.5,'DisplayName','Exact Solution')
plot(x,x*0+dissipation(end),'LineWidth',1.5,'DisplayName','Theoretical Dissipation')
grid on
legend('NumColumns',3,'Location','southoutside','FontSize',15,'Box','off')
xlabel('$x$','fontsize',15)
ylabel('$u(40T,x)$','fontsize',15)

save_name = '../plots/exercise03/dissipation.png';
exportgraphics(gcf,save_name,'Resolution',300);
%%
m = 126;
figure
plot(utrue(1,:))
hold on
for i = 2:m
    plot(utrue(m,:))
end

t(m) - t(1);

%%
clear; clc; close all;

N = 100;
tmax = 20;
a = 0.5;
k = ((2/N)/(a)); % Largest value while preserving stability
M = ceil(tmax/k);
M = ceil(M*1.1);

U = ForwardTimeBackwardSpace(N,M,tmax,a);

plotSolution(U,M,N,tmax,a);

% Save plot of solution
save_name = '../plots/exercise03/dissipation_dispersion_2d.png';
exportgraphics(gcf,save_name,'Resolution',300);

x = linspace(-1,1,N+1);
k = tmax/M;
t = linspace(0,tmax,M+1);
[X, T] = meshgrid(x,t);
utrue = fun(X,T);

figure('Renderer', 'painters', 'Position', [400 400 1000 500]);
tiledlayout(2,3,'TileSpacing','compact');

for i=1:6
    idx = 80*i;
    nexttile
    plot(x,utrue(idx,:),'r-','LineWidth',1.5,'Displayname','Exact Solution')
    hold on
    plot(x,U(idx,:),'b-','LineWidth',1.5,'Displayname','Numerical Solution')
    grid on
    title(sprintf("$t=%.3f$",t(idx)))

end

l = legend('box','off','fontsize',15,'NumColumns',2);
l.Layout.Tile = 'south';

% Save plot of solution
save_name = '../plots/exercise03/dissipation_dispersion_1d.png';
exportgraphics(gcf,save_name,'Resolution',300);


phaseAmplitudeError(utrue,U,x,t)