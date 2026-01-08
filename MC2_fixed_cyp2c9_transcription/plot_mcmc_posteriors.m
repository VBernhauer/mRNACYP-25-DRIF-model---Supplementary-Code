function [] = plot_mcmc_posteriors()

    %%% command:
    %%% plot_mcmc_posteriors()
   
    clc;
    close all
    set(0,'DefaultFigureVisible','on');

    if ~exist('./figures', 'dir')
        mkdir('./figures')
    end

    chains = [];
    for ii = 1:5
        jjchains = load(strcat('./chains/chains_',num2str(ii),'.mat'));
        jjchains = jjchains.chains(:,:);
        chains = [chains;jjchains];
    end

    %%% Aggregated posterior plots %%%
    fontsize = 12;
    leg = {'$k_{pxr,drif,cyp}$',...
           '$k_{{mRNA}^{fold}_{cyp3a4}}$',...
           '$k_{{mRNA}^{fold}_{cyp2b6}}$'};
    fcolor = {[0 0 0],[1 0 0],[0 0 1]};
    figure(1);
    tiledplot = tiledlayout(1,1,'TileSpacing','Compact');
    set(gcf, 'Position',  [300, 100, 600, 400]);
    ax(1) = nexttile();
    set(ax(1),...
        'box','on',...
        'FontSize',10);
    set(gca,'TickLength',[0.025, 0.01])
    hold on;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for nn = 1:3
        histogram(log10(chains(:,nn)),...
              'FaceColor', 'none',...
              'EdgeColor', fcolor{nn},...
              'FaceAlpha', 1,...
              'EdgeAlpha', 1,...
              'Normalization', 'probability',...
              'DisplayStyle', 'stairs',...
              'LineWidth',1.5);
        hold on;
    end
    l=legend(leg);
    set(l,'Location','northwest');
    set(l,'FontSize',fontsize);
    set(l,'interpreter','latex');
       
    xlabel(tiledplot,'Log_{10} value','FontSize',fontsize);
    ylabel(tiledplot,'Normalized frequency','FontSize',fontsize);
    tiledplot.TileSpacing = 'compact';
    tiledplot.Padding = 'compact';

    savefig(strcat('figures/posteriors_fixed_cyp2c9_transcription.fig'));
    exportgraphics(gcf,'figures/posteriors_fixed_cyp2c9_transcription.png');
    exportgraphics(gcf,'figures/posteriors_fixed_cyp2c9_transcription.eps','ContentType','vector');               

    end