function plotSolution(U,M,N,k)

fontsize = 13;
x = linspace(-1,1,N+1);
t = linspace(0,k*M,M+1);
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

sgtitle(sprintf('$N=%d$, $M=%d$, $k=%.4f$',N,M,k),'FontSize',fontsize+3,'interpreter','latex')