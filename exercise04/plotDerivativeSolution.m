function plotDerivativeSolution(U,N,M,tmax)
fontsize = 13;
xbar = 0;
x = linspace(-1,1,N+1);
t = linspace(0,tmax,M+1);
y = zeros(5,1);
nexttile;
imagesc(x,t,U)
axis square
title(sprintf("$N = %d$, $M=%d$",N,M),'fontsize',13)
colorbar
clim([-1,1])
xlabel('$x$','FontSize',fontsize)
ylabel('$t$','FontSize',fontsize)
nexttile;
u = U(end,:);
plot(x,u,'.-')
xlabel('$x$','FontSize',fontsize)
ylabel('$u(t^*,x)$','FontSize',fontsize)
grid on
nabo = [2 5 10];
for j = 1:3
    y(j) = computeDerivative(u,x,nabo(j),nabo(j),N,xbar);
end
title(sprintf("$\\partial_xu|_{{x=0}}=\\{{%.2f,%.2f,%.2f\\}}$",y(1),y(2),y(3)),'fontsize',13)
nexttile;
h = 2/N;
plot(x(2:end),diff(u)/h,'.-')
xlim([-0.05,0.05])
grid on
xlabel('$x$','FontSize',fontsize)
ylabel('$\partial_xu(t^*,x)$','FontSize',fontsize)