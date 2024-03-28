function plot_varying_delta(deltas, T_cell, Y_cell, H_cell, E_cell, numSteps,T_cell_ode23, Y_cell_ode23,reps,aeps)

figure('Renderer', 'painters', 'Position', [400 400 1200 700]);
tiledlayout(3,3,'TileSpacing','compact');
fontsize = 13;
for i = 1:length(deltas)
    delta = deltas(i);
    nexttile;
    plot(T_cell{i}*delta,Y_cell{i},'.-')
    hold on
    plot(T_cell_ode23{i}*delta,Y_cell_ode23{i},'--')
    title(sprintf('$\\delta = %.5f$',delta),'FontSize',fontsize+2)
    ylim([-0.05,1.05])
    xlim([0,2])
    grid on
    xlabel('$t\cdot \delta$','FontSize',fontsize)
    ylabel('$y(t)$','FontSize',fontsize)
end

sgtitle(sprintf('\\texttt{{reps = %.4e, aeps = %.4e}}',reps,aeps),'fontsize',16,'interpreter','latex')

exportgraphics(gcf,sprintf('../plots/exercise01/solutions_varying_delta_reps_%d_aeps_%d.png',-log10(reps),-log10(aeps)),'Resolution',300);

figure('Renderer', 'painters', 'Position', [400 400 1200 700]);
tiledlayout(3,3,'TileSpacing','compact');

for i = 1:length(deltas)
    delta = deltas(i);
    nexttile;
    T = T_cell{i};
    E = E_cell{i};
    plot(T(2:end)*delta,E,'.-')
    title(sprintf('$\\delta = %.5f$',delta),'FontSize',fontsize+2)
    xlim([0,2])
    grid on
    xlabel('$t\cdot \delta$','FontSize',fontsize)
    ylabel('Error','FontSize',fontsize)
end
sgtitle(sprintf('\\texttt{{reps = %.4e, aeps = %.4e}}',reps,aeps),'fontsize',16,'interpreter','latex')

exportgraphics(gcf,sprintf('../plots/exercise01/error_varying_delta_%d_aeps_%d.png',-log10(reps),-log10(aeps)),'Resolution',300);


figure('Renderer', 'painters', 'Position', [400 400 1200 700]);
tiledlayout(3,3,'TileSpacing','compact');

for i = 1:length(deltas)
    delta = deltas(i);
    nexttile;
    T = T_cell{i};
    H = H_cell{i};
    plot(T*delta,H,'.-')
    title(sprintf('$\\delta = %.5f$',delta),'FontSize',fontsize+2)
    xlim([0,2])
    grid on
    xlabel('$t\cdot \delta$','FontSize',fontsize)
    ylabel('Step Size','FontSize',fontsize)
end
sgtitle(sprintf('\\texttt{{reps = %.4e, aeps = %.4e}}',reps,aeps),'fontsize',16,'interpreter','latex')

exportgraphics(gcf,sprintf('../plots/exercise01/step_size_varying_delta_%d_aeps_%d.png',-log10(reps),-log10(aeps)),'Resolution',300);