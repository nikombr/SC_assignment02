
%% Parabolic

clear; close all; clc;

h = [20 40 80];
H = cell(3,1);

figure('Renderer', 'painters', 'Position', [400 400 1000 1000]);
tiledlayout(3,3,'TileSpacing','compact');

for i = 1:3
    H{1} = [h(i) 1];
    H{2} = [h(i) h(i) h(i) 1];
    H{3} = [h(i) h(i) h(i) h(i) h(i) 1];
    for j = 1:3
        load(sprintf('results/parabolic_%d_%d.mat',i-1,j-1))
            nexttile
            [X, T] = meshgrid(x,t);
            utrue = funParabolic(X,T);
            imagesc(x,t,U)
            axis square
            xlabel('$x$','FontSize',15)
            ylabel('$t$','FontSize',15)
            
            clim([-1.5,3])

            str = "Layers = \{2";
            for Hh = H{j}
                str = sprintf("%s, %d",str,Hh);
            end
            str = sprintf("%s\\}",str);

            disp(sprintf("%s, Loss = %.4f, Global Error = %.4f",str, loss(end),max(max(abs(U-utrue)))))

            title(str,'fontsize',15)

    end 
end

c = colorbar;
c.Layout.Tile = 'east';

exportgraphics(gcf,'../plots/exercise05/parabolic_hypertest_solution.png','Resolution',300);

figure('Renderer', 'painters', 'Position', [400 400 1000 1000]);
tiledlayout(3,3,'TileSpacing','compact');

for i = 1:3
    H{1} = [h(i) 1];
    H{2} = [h(i) h(i) h(i) 1];
    H{3} = [h(i) h(i) h(i) h(i) h(i) 1];
    for j = 1:3
        load(sprintf('results/parabolic_%d_%d.mat',i-1,j-1))
            nexttile
            plot(loss)
            ylabel('Loss','FontSize',15)
            xlabel('Epoch','FontSize',15)
            axis square
            grid on
            str = "Layers = \{2";
            for Hh = H{j}
                str = sprintf("%s, %d",str,Hh);
            end
            str = sprintf("%s\\}",str);
            title(str,'fontsize',15)

    end 
end

exportgraphics(gcf,'../plots/exercise05/parabolic_hypertest_loss.png','Resolution',300);

%% Parabolic convergence

clear; clc; close all;
tmax = 0.5;
M = 10000;
k = tmax/M;
N = round(linspace(10,200,10));
error_spatial = zeros(length(N),1);
error = zeros(length(N),1);
h = 2./N;

for i = 1:length(N)

    load(sprintf('results/parabolic_convergence_spatial_N_%d.mat',N(i)));

    [X, T] = meshgrid(x,t);
    utrue = funParabolic(X,T);

    error_spatial(i) = max(max(abs(U-utrue)));


end

for i = 1:length(N)

    load(sprintf('results/parabolic_convergence_combined_N_%d.mat',N(i)));

    [X, T] = meshgrid(x,t);
    utrue = funParabolic(X,T);

    error(i) = max(max(abs(U-utrue)));


end



figure('Renderer', 'painters', 'Position', [400 400 1000 300]);
tiledlayout(1,2,'TileSpacing','compact');
nexttile;
loglog(h,error_spatial,'LineWidth',1.5,'DisplayName',"$||U-u(x,t)||_\infty$")
hold on
loglog(h,h.^2,'--','LineWidth',1.5,'DisplayName','$\mathcal{O}(h^2)$')
grid on
legend('FontSize',15,'Location','southeast')
xlabel('$h$','FontSize',15)
ylabel('Error','FontSize',15)
title(sprintf('$k=%.2f\\cdot 10^{-5}$',tmax/10000*10^(5)),'FontSize',15)
nexttile;
loglog(h,error,'LineWidth',1.5,'DisplayName',"$||U-u(x,t)||_\infty$")
hold on
loglog(h,h.^2,'--','LineWidth',1.5,'DisplayName','$\mathcal{O}(h^2)=\mathcal{O}(k)$')
grid on
legend('FontSize',15,'Location','southeast')
title('$k=h^2$','FontSize',15)

xlabel('$h$','FontSize',15)
ylabel('Error','FontSize',15)
exportgraphics(gcf,'../plots/exercise05/parabolic_convergence.png','Resolution',300);


%% Hyperbolic

clear; close all; clc;

load 'results/hyperbolic.mat'

[X, T] = meshgrid(x,t);
utrue = funHyperbolic(X,T);

[X, T] = meshgrid(x,t_pred);
utrue_pred = funHyperbolic(X,T);

figure('Renderer', 'painters', 'Position', [400 400 800 300]);

plot(loss)
ylabel('Loss','FontSize',15)
xlabel('Epoch','FontSize',15)
axis square
grid on

figure('Renderer', 'painters', 'Position', [400 400 1000 600]);
tiledlayout(2,3,'TileSpacing','compact');

nexttile
imagesc(x,t,utrue)
axis square
xlabel('$x$','FontSize',15)
ylabel({"Training Window",'$t$'},'FontSize',15)
colorbar
clim([-1,1])

title("\textbf{True Solution}",'FontSize',15)

nexttile
imagesc(x,t,U)
axis square
xlabel('$x$','FontSize',15)
ylabel('$t$','FontSize',15)
colorbar
clim([-1,1])

title("\textbf{PINN Solution}",'FontSize',15)


nexttile
imagesc(x,t,abs(U-utrue))
axis square
xlabel('$x$','FontSize',15)
ylabel('$t$','FontSize',15)
colorbar
title("\textbf{Error}",'FontSize',15)


nexttile
imagesc(x,t_pred,utrue_pred)
axis square
xlabel('$x$','FontSize',15)
ylabel({"Testing Window",'$t$'},'FontSize',15)
colorbar
clim([-1,1])

nexttile
imagesc(x,t_pred,U_pred)
axis square

xlabel('$x$','FontSize',15)
ylabel('$t$','FontSize',15)
colorbar
clim([-1,1])

nexttile
imagesc(x,t_pred,abs(U_pred-utrue_pred))
axis square
xlabel('$x$','FontSize',15)
ylabel('$t$','FontSize',15)
colorbar

exportgraphics(gcf,'../plots/exercise05/hyperbolic.png','Resolution',300);


figure('Renderer', 'painters', 'Position', [400 400 500 250]);
plot(loss)
ylabel('Loss','FontSize',15)
xlabel('Epoch','FontSize',15)
grid on

exportgraphics(gcf,'../plots/exercise05/hyperbolic_loss.png','Resolution',300);

%% Advection

clear; close all; clc;

load 'results/advection.mat'


figure('Renderer', 'painters', 'Position', [400 400 1200 300]);
tiledlayout(1,3,'TileSpacing','compact');

nexttile
plot(loss)
ylabel('Loss','FontSize',15)
xlabel('Epoch','FontSize',15)
grid on

nexttile
imagesc(x,t,U)
xlabel('$x$','FontSize',15)
ylabel('$t$','FontSize',15)
colorbar
clim([-1,1])

nexttile
u = U(end,:);
plot(x,u)
N = length(x);
xbar =0;
grid on
y = computeDerivative(u,x,5,5,N,xbar)
title(sprintf('$\\partial_x u|_{x=0}\\approx %.3f$',y),'FontSize',15)

exportgraphics(gcf,'../plots/exercise05/advection.png','Resolution',300);
