
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
U = AdvectionDiffusion(@boundaryFun,N,M,tmax,epsilon,"uniform");

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
    U = AdvectionDiffusion(@boundaryFun,N,M,tmax,epsilon,"uniform");
    plotSolution(@tanhFun,U,M,N,tmax,false,epsilon);
end

exportgraphics(gcf,'../plots/exercise04/illustration.png','Resolution',300);

%% Convergence

clear; clc; close all;

numEval = 60;
epsilon = 0.1;

Nmin = 50;
Nmax = 1000;
M = 10000;
tmax = 0.1;
N = round(linspace(Nmin,Nmax,numEval));
h = 2./N;
t = linspace(0,tmax,M+1);
error_spatial = zeros(length(N),1);

for i = 1:length(N)
    i
    U = AdvectionDiffusion(@boundaryFun,N(i),M,tmax,epsilon,'uniform');
    x = linspace(-1,1,N(i)+1);
    [X, T] = meshgrid(x,t);
    utrue = tanhFun(X,T,epsilon);
    error_spatial(i) = max(max(abs(U - utrue)));
    clear U
end
%%
Nmin = 50;
Nmax = 500;
tmax = 0.1;
N = round(linspace(Nmin,Nmax,numEval));
h = 2./N;
k = h.^2;
M = ceil(tmax./k);

error = zeros(length(N),1);

for i = 1:length(N)
    i
    U = AdvectionDiffusion(@boundaryFun,N(i),M(i),tmax,epsilon,"uniform");
    t = linspace(0,tmax,M(i)+1);
    x = linspace(-1,1,N(i)+1);
    [X, T] = meshgrid(x,t);
    utrue = tanhFun(X,T,epsilon);
    error(i) = max(max(abs(U - utrue)));
    clear U
end
%%
figure('Renderer', 'painters', 'Position', [400 400 1000 300]);
tiledlayout(1,2,'TileSpacing','compact');
nexttile;
loglog(h,error_spatial,'LineWidth',1.5,'DisplayName',"$||U-u(x,t)||_\infty$")
hold on
loglog(h,h,'--','LineWidth',1.5,'DisplayName','$\mathcal{O}(h)$')
grid on
legend('FontSize',15,'Location','northwest')
xlabel('$h$','FontSize',15)
ylabel('Error','FontSize',15)
title(sprintf('$k=%.2f\\cdot 10^{-5}$',tmax/10000*10^(5)),'FontSize',15)
nexttile;
loglog(h,error,'LineWidth',1.5,'DisplayName',"$||U-u(x,t)||_\infty$")
hold on
loglog(h,h+h.^2,'--','LineWidth',1.5,'DisplayName','$\mathcal{O}(h+h^2)=\mathcal{O}(h+k)$')
grid on
legend('FontSize',15,'Location','northwest')
title('$k=h^2$','FontSize',15)

xlabel('$h$','FontSize',15)
ylabel('Error','FontSize',15)
exportgraphics(gcf,'../plots/exercise04/convergence.png','Resolution',300);

%% Homogenous boundary, convergence

clear; close all; clc;

numEval = 20;

epsilon = 0.01/pi;
Nmax = 13000; % Maximal N to keep stability
Nmin = 600; % Minimum N to keep stability
M = 120000; % Fixed number of time steps to have stability at Nmax
dudx = -152.00516; % True derivative

tstar = 1.6037/pi; % Time of evaluation of derivative
tmax = tstar;
xbar = 0;
Ns = ceil(linspace(Nmin,Nmax,numEval));
h = 2./Ns;
k = h/10;
Ms = ceil(tmax./k);
Ms = 10*Ns

error = zeros(length(Ns),1);

for i = 1:length(Ns)
    N = Ns(i)
    M = Ms(i)
    x = linspace(-1,1,N+1);
    t = linspace(0,tmax,M+1);
    U = AdvectionDiffusion(@homogeneousBoundaryFun,N,M,tmax,epsilon,"uniform");
    u = U(end,:);
    y = computeDerivative(u,x,5,5,N,xbar);
    error(i) = abs(dudx - y);
    
    clear U
end
%%
figure('Renderer', 'painters', 'Position', [400 400 1000 300]);
loglog(h(1:end-2),error(1:end-2),'.-','DisplayName',"$||U-u(x,t)||_\infty$",'linewidth',1.5,'markersize',20)
hold on
%loglog(h,h.^2,'DisplayName','$\mathcal{O}(h^2)$')
grid on
legend('FontSize',15,'Location','northwest')
xlabel('$h$','FontSize',15)
ylabel('Error','FontSize',15)

exportgraphics(gcf,'../plots/exercise04/derivative_estimate_convergence2.png','Resolution',300);
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
    Us{i} = AdvectionDiffusion(@homogeneousBoundaryFun,N,M,tmax,epsilon,"uniform"); 
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
    Us{i} = AdvectionDiffusion(@homogeneousBoundaryFun,N,M,tmax,epsilon,"uniform"); 
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

N = 300;
M = 20000;
a = [0.9 0.7 0.5 0.3];
xbar = 0;

figure('Renderer', 'painters', 'Position', [400 400 1000 1000]);
tiledlayout(4,3,'TileSpacing','compact');

for i = 1:4

    [U,x,t] = AdvectionDiffusion(@homogeneousBoundaryFun,N,M,tmax,epsilon,'nonuniform',a(i));

    nexttile;
    plot(x(2:end),diff(x),'-')

    nexttile;
    imagesc(x,t,U)
    axis square
    
    title(sprintf("$N = %d$, $M=%d$",N,M),'fontsize',13)
    colorbar
    clim([-1,1])
    fontsize=13;
    xlabel('$x$','FontSize',fontsize)
    ylabel('$t$','FontSize',fontsize)

    nexttile;
    u = U(end,:);
    plot(x,u)
    grid on

    y = computeDerivative(u,x,5,5,N,xbar)

end

%% Higher order

clear; close all; clc;
epsilon = 0.01/pi;
tstar = 1.6037/pi; % Time of evaluation of derivative
tmax = tstar;

N = 1100; % Det her virker helt vildt godt!!
M = 3000;
xbar = 0;

figure('Renderer', 'painters', 'Position', [400 400 1000 300]);
tiledlayout(1,2,'TileSpacing','compact');

[U,x,t] = AdvectionDiffusion(@homogeneousBoundaryFun,N,M,tmax,epsilon,'higher');


nexttile;
imagesc(x,t,U)
axis square

title(sprintf("$N = %d$, $M=%d$",N,M),'fontsize',13)
colorbar
clim([-1,1])
fontsize=13;
xlabel('$x$','FontSize',fontsize)
ylabel('$t$','FontSize',fontsize)

nexttile;
u = U(end,:);
plot(x,u)
grid on

y = computeDerivative(u,x,5,5,N,xbar)



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
