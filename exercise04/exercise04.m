
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
%% Testing
clear; close all; clc;

N = 5;
tmax = 0.1;   
epsilon = 0.1;
k = ((2/N)/(epsilon)); % Largest value while preserving stability
k = 0.001;
k = 0.00001;
%M = ceil(tmax/k);
M = 10;
kmax = ((2/N)/(epsilon))
k = tmax/M
U = AdvectionDiffusion(@boundaryFun,N,M,tmax,epsilon);

plotSolution(@tanhFun,U,M,N,tmax,true,epsilon);

%% Illustration

clear; close all; clc;

epsilon = 0.1;
tmax = 0.4;
Ns = [10 40 160];
Ms = [5 50 570]; % Chosen smallest to keep stability
k = tmax./Ms
kmax = ((2./Ns)/(epsilon))
figure('Renderer', 'painters', 'Position', [400 400 1000 900]);
tiledlayout(3,3,'TileSpacing','compact');
for i = 1:3
    N = Ns(i);
    M = Ms(i);
    U = AdvectionDiffusion(@boundaryFun,N,M,tmax,epsilon);
    plotSolution(@tanhFun,U,M,N,tmax,false,epsilon);
end

exportgraphics(gcf,'../plots/exercise04/illustration.png','Resolution',300);

%% Convergence

clear; clc; close all;

numEval = 60;
epsilon = 0.1;

Nmin = 50;
Nmax = 1500;
M = 15000;
tmax = 0.1;
Ns = round(linspace(Nmin,Nmax,numEval));
hs = 2./Ns;
t = linspace(0,tmax,M+1);
error_spatial = zeros(length(Ns),1);

for i = 1:length(Ns)
    N = Ns(i)
    U = AdvectionDiffusion(@boundaryFun,N,M,tmax,epsilon);
    x = linspace(-1,1,N+1);
    [X, T] = meshgrid(x,t);
    utrue = tanhFun(X,T,epsilon);
    error_spatial(i) = max(max(abs(U - utrue)));
    clear U
end


N = 300;
Mmax = 10000;
Mmin = 300;
Ms = round(linspace(Mmin,Mmax,numEval));
ks = tmax./Ms;
x = linspace(-1,1,N+1);
error_time = zeros(length(Ns),1);

for i = 1:length(Ms)
    M = Ms(i)
    U = AdvectionDiffusion(@boundaryFun,N,M,tmax,epsilon);
    t = linspace(0,tmax,M+1);
    [X, T] = meshgrid(x,t);
    utrue = fun(X,T,epsilon);
    error_time(i) = max(max(abs(U - utrue)));
    clear U
end


figure('Renderer', 'painters', 'Position', [400 400 1000 300]);
tiledlayout(1,2,'TileSpacing','compact');
nexttile;
loglog(hs,error_spatial,'DisplayName',"$||U-u(x,t)||_\infty$")
hold on
loglog(hs,hs,'DisplayName','$h$')
grid on
legend('FontSize',15,'Location','northwest')
xlabel('$h$','FontSize',15)
ylabel('Error','FontSize',15)
nexttile;
loglog(ks,error_time,'DisplayName',"$||U-u(x,t)||_\infty$")
hold on
%loglog(ks,ks,'DisplayName','$k$')
grid on
legend('FontSize',15)
xlabel('$k$','FontSize',15)
ylabel('Error','FontSize',15)
exportgraphics(gcf,'../plots/exercise04/convergence.png','Resolution',300);

%% Homogenous boundary, convergence of spatial variable

clear; close all; clc;

numEval = 20;

epsilon = 0.01/pi;
Nmax = 12000; % Maximal N to keep stability
Nmin = 600; % Minimum N to keep stability
M = 120000; % Fixed number of time steps to have stability at Nmax
dudx = -152.00516; % True derivative

tstar = 1.6037/pi; % Time of evaluation of derivative
tmax = tstar;
xbar = 0;
Ns = linspace(Nmin,Nmax,numEval);
Ms = Ns*10

error = zeros(length(Ns),1);

for i = 1:length(Ns)
    N = Ns(i)
    M = Ms(i)
    x = linspace(-1,1,N+1);
    t = linspace(0,tmax,M+1);
    U = AdvectionDiffusion(@homogeneousBoundaryFun,N,M,tmax,epsilon);
    u = U(end,:);
    y = computeDerivative(u,x,5,5,N,xbar);
    error(i) = abs(dudx - y);

    clear U
end
%%
hs = 2./Ns;
figure('Renderer', 'painters', 'Position', [400 400 1000 300]);
loglog(hs,error,'.-','DisplayName',"$||U-u(x,t)||_\infty$")
hold on
loglog(hs,1e8*hs.^2,'DisplayName','$\mathcal{O}(h^2)$')
grid on
legend('FontSize',15,'Location','northwest')
xlabel('$h$','FontSize',15)
ylabel('Error','FontSize',15)

exportgraphics(gcf,'../plots/exercise04/derivative_estimate_convergence.png','Resolution',300);
%% Illustration, varying N
clear; close all; clc;

epsilon = 0.01/pi;
Nmax = 9000; % Maximal N to keep stability
Nmin = 600; % Minimum N to keep stability
M = 90000; % Fixed number of time steps to have stability at Nmax
dudx = -152.00516; % True derivative

tstar = 1.6037/pi; % Time of evaluation of derivative
tmax = tstar;
Ns = linspace(Nmin,Nmax,3);
Us = cell(3,1);

for i = 1:length(Ns)
    N = Ns(i)
    Us{i} = AdvectionDiffusion(@homogeneousBoundaryFun,N,M,tmax,epsilon); 
end

figure('Renderer', 'painters', 'Position', [400 400 1000 600],'visible','off');
tiledlayout(3,3,'TileSpacing','compact');

for i = 1:length(Ns)
    N = Ns(i);
    U = Us{i};
    plotDerivativeSolution(U,N,M,tmax)
end

exportgraphics(gcf,'../plots/exercise04/derivative_estimate_illustration.png','Resolution',300);

%% Illustration, varying M

clear; close all; clc;

epsilon = 0.01/pi;
N = 2000; 
Mmin = 4000; 
Mmax = 90000; 
dudx = -152.00516; % True derivative

tstar = 1.6037/pi; % Time of evaluation of derivative
tmax = tstar;
Ms = linspace(Mmin,Mmax,3);
Us = cell(3,1);

for i = 1:length(Ms)
    M = Ms(i)
    Us{i} = AdvectionDiffusion(@homogeneousBoundaryFun,N,M,tmax,epsilon); 
end

figure('Renderer', 'painters', 'Position', [400 400 1000 600],'visible','off');
tiledlayout(3,3,'TileSpacing','compact');

for i = 1:length(Ms)
    M = Ms(i);
    U = Us{i};
    plotDerivativeSolution(U,N,M,tmax)
end

exportgraphics(gcf,'../plots/exercise04/derivative_estimate_illustration_time.png','Resolution',300);

%% Non-uniform

clear; close all; clc;
epsilon = 0.01/pi;
tstar = 1.6037/pi; % Time of evaluation of derivative
tmax = tstar;
N = 100;
Nfine = 1000;
M = 4000;

N = 600;
Nfine = 600;
M = 90000;

[U,x,t] = AdvectionDiffusion(@homogeneousBoundaryFun,N,M,tmax,epsilon,'nonuniform',Nfine);
figure
plot(x,x,'o-')
figure;
imagesc(x,t,U)
axis square
title(sprintf("$N = %d$, $M=%d$",N,M),'fontsize',13)
colorbar
clim([-1,1])
fontsize=13;
xlabel('$x$','FontSize',fontsize)
ylabel('$t$','FontSize',fontsize)
%% tester lige
clear; clc; close all;

m = 10

x = linspace(1,m+2,m+2)

A = zeros(m+2);

for i=2:m+1
    A(i,i-1:i+1) = fdcoeffV(2, x(i), x((i-1):(i+1)));
end

FTCS2 = A(2:end-1,2:end-1); % Exercise 2

for i=2:m+1
    A(i,i-1:i+1) = 2*fdcoeffV(1, x(i), x((i-1):(i+1)));
end

FTCS = A(2:end-1,2:end-1) % FTCS

for i=2:m+1
    A(i,i-1:i) = 2*fdcoeffV(1, x(i), x((i-1):(i)));
end

FTBS = A(2:end-1,2:end-1) % FTBS


%% Matlab solution

clear; close all; clc;
Ns = [100, 1000, 10000];
Ms = [100, 1000, 10000];
epsilon = 0.01/pi;
tmax = 1.6037/pi;
figure('Renderer', 'painters', 'Position', [400 400 1000 900]);
tiledlayout(3,3,'TileSpacing','compact');
for i = 1:3
    M = Ms(i)
    N = Ns(i)
    x = linspace(-1,1,N+1);
    t = linspace(0,tmax,M+1);
    m = 0;
    sol = pdepe(m,@heatpde,@heatic,@heatbc,x,t);
    h = 2/N;
    plotDerivativeSolution(sol,N,M,tmax)
end
exportgraphics(gcf,'../plots/exercise04/matlab_solution.png','Resolution',300);

function [c,f,s] = heatpde(x,t,u,dudx)  
  epsilon = 0.01/pi;
  c = 1;
  f = epsilon*dudx;
  s = -dudx*u;
end

function u0 = heatic(x)
u0 = -sin(pi*x);
end

function [pl,ql,pr,qr] = heatbc(xl,ul,xr,ur,t)
pl = ul;
ql = 0;
pr = ur;
qr = 0;
end

%%
