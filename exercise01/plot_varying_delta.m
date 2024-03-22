function plot_varying_delta(deltas, T_cell, Y_cell, H_cell, E_cell, numSteps,T_cell_ode23, Y_cell_ode23)

figure('Renderer', 'painters', 'Position', [400 400 1200 700]);
tiledlayout(3,3,'TileSpacing','compact');
fontsize = 13;
for i = 1:length(deltas)
    delta = deltas(i);
    nexttile;
    plot(T_cell{i},Y_cell{i},'.-')
    hold on
    plot(T_cell_ode23{i},Y_cell_ode23{i},'--')
    title(sprintf('$\\delta = %.5f$',delta),'FontSize',fontsize+2)
    ylim([-0.05,1.05])
    grid on
    xlabel('$t$','FontSize',fontsize)
    ylabel('$y(t)$','FontSize',fontsize)
end
exportgraphics(gcf,'../plots/exercise01/solutions_varying_delta.png','Resolution',300);

figure('Renderer', 'painters', 'Position', [400 400 1200 700]);
tiledlayout(3,3,'TileSpacing','compact');
fontsize = 13;
for i = 1:length(deltas)
    delta = deltas(i);
    nexttile;
    T = T_cell{i};
    E = E_cell{i};
    plot(T(2:end),E,'.-')
    title(sprintf('$\\delta = %.5f$',delta),'FontSize',fontsize+2)
    %ylim([-0.05,1.05])
    grid on
    xlabel('$t$','FontSize',fontsize)
    ylabel('$y(t)$','FontSize',fontsize)
end