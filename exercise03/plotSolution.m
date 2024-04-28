function plotSolution(U,M,N,tmax,a)

k = tmax/M;
h = 2/N;
fontsize = 13;
x = linspace(-1,1,N+1);
t = linspace(0,tmax,M+1);
[X, T] = meshgrid(x,t);
utrue = fun(X,T);


figure('Renderer', 'painters', 'Position', [400 400 1000 300]);
tiledlayout(1,3,'TileSpacing','compact');
nexttile;
imagesc(x,t,utrue)
colorbar
clim([-1.5,3])
axis square
title('\textbf{True Solution}','FontSize',fontsize+2)
xlabel('$x$','FontSize',fontsize)
ylabel('$t$','FontSize',fontsize)

nexttile;
imagesc(x,t,U)
colorbar
clim([-1.5,3])
axis square
xlabel('$x$','FontSize',fontsize)
ylabel('$t$','FontSize',fontsize)

title('\textbf{Computed Solution}','FontSize',fontsize+2)

nexttile;
imagesc(x,t,abs(U-utrue))
axis square
colorbar
title('\textbf{Absolute Error}','FontSize',fontsize+2)
xlabel('$x$','FontSize',fontsize)
ylabel('$t$','FontSize',fontsize)

if k <= h/a
    sgtitle(sprintf('$N=%d$, $M=%d$, $k=%.5f\\leq %.5f=h/a$',N,M,k,h/a),'FontSize',fontsize+4,'interpreter','latex')
else
    sgtitle(sprintf('$N=%d$, $M=%d$, $k=%.5f\\not\\leq %.5f=h/a$',N,M,k,h/a),'FontSize',fontsize+4,'interpreter','latex')
end