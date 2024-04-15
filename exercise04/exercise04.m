
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');

clear; close all; clc;

N = 100;
tmax = 1;   
epsilon = 0.1;
k = ((2/N)/(epsilon)); % Largest value while preserving stability
k = 0.001;
k = 0.00001;
%M = ceil(tmax/k);
M = 200;

U = AdvectionDiffusion(@boundaryFun,N,M,k,epsilon);

plotSolution(@fun,U,M,N,k,epsilon);

%%

clear; close all; clc;

N = 11000;
h = 1/(N+1);
tstar = 1.6037/pi;
tmax = tstar;
epsilon = 0.01/pi;
%epsilon = 0.01
%k = ((2/N)/(epsilon)); % Largest value while preserving stability

%k = min(h^2/(2*epsilon),2*h)
%k = 0.001;
%k = 0.001;
%k = 0.00001;
%M = ceil(tmax/k);
%M = 20;

M = 110000;
k = tmax/M;
hm = min(h^2/(2*epsilon),2*h)

M = ceil(tmax/k);

U = AdvectionDiffusion(@homogeneousBoundaryFun,N,M,k,epsilon);


x = linspace(-1,1,N+1);
t = linspace(0,M*k,M+1);
imagesc(x,t,U)
hold on
%plot(x,x*0+tstar,'k--','LineWidth',1.5)
%plot(t*0,t,'k--','LineWidth',1.5)
%plot(0,tstar,'r.','markersize',30)
colorbar


figure;
%u = interp2(x,t,U,x,tstar);
u = U(end,:)
plot(x,u,'.-')
grid on
xbar = 0

nabo = 40
y = computeDerivative(u,x,5,5,N,xbar)

%%
N = 2000;
epsilon = 0.01/pi;
x = linspace(-1,1,N+1);
t = linspace(0,M*k,M+1);
xbar = 0;
m = 0;
sol = pdepe(m,@heatpde,@heatic,@heatbc,x,t);
h = 2/N;

figure
imagesc(sol-U)

usol = sol(end,:);
figure('Renderer', 'painters', 'Position', [400 400 1000 300]);
plot(x,usol,'.-')
hold on 
plot(x,u,'.-')


figure('Renderer', 'painters', 'Position', [400 400 1000 300]);
tiledlayout(1,3,'TileSpacing','compact');
nexttile;
imagesc(x,t,sol)

nexttile;
plot(x,usol,'.-')
hold on 
plot(x,u,'.-')
nexttile;
plot(x(2:end),diff(u)/h,'.-')

y = computeDerivative(u,x,3,3,N,xbar)

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
