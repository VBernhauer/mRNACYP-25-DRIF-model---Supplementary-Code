function [] = plot_practical_identifiability()

    clc;
    clear all;
    close all;

    if ~exist('./figures', 'dir')
        mkdir('./figures')
    end

    k_pxrmax      = load('./k_pxrmax.mat');
    k_mRNAcyp3a4  = load('./k_mRNAcyp3a4.mat');
    k_mRNAcyp2c9  = load('./k_mRNAcyp2c9.mat');
    k_mRNAcyp2b6  = load('./k_mRNAcyp2b6.mat');

    data = {k_pxrmax, k_mRNAcyp3a4, k_mRNAcyp2c9, k_mRNAcyp2b6};

    MLvalues = load('./MLvalues.mat');

    %%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    figlabs = {'(A)','(B)','(C)','(D)'};
    xLabels =  {'$k_\mathrm{{pxr,max}}\, \mathrm{(\mu M^{-1}h^{-1})}$',...
                '$k_\mathrm{{mRNA_{cyp3a4}^{fold}}}\, \mathrm{(h^{-1})}$',...
                '$k_\mathrm{{mRNA_{cyp2c9}^{fold}}}\, \mathrm{(h^{-1})}$',...
                '$k_\mathrm{{mRNA_{cyp2b6}^{fold}}}\, \mathrm{(h^{-1})}$'};
    yLims = {[24 28],...
             [10 120],...
             [20 60],...
             [20 40]};
    fontsize    = 14;
    tiledplot = tiledlayout(2,2,'TileSpacing','Compact');
    set(gcf, 'Position',  [300, 100, 800, 600]);
    for aa = 1:4
        fval = data{aa};
        xXaxis = fval.funval(:,1);
        yYaxis = fval.funval(:,2);
        ax(aa) = nexttile(aa);            
        set(ax(aa),...
            'box','on',...
            'YLim',yLims{aa},...
            'XTickLabelRotation',0,...
            'FontSize',fontsize);
        set(gca,'TickLength',[0.025, 0.01])
        hold on
        plot(ax(aa),xXaxis,yYaxis,'Color',[0 0 0 1],'LineWidth',1,'LineStyle','-')
        hold on
        plot(ax(aa),xXaxis,(MLvalues.MLvalues(1,1)+0.5*chi2inv(0.95,1))*ones(1,length(xXaxis)),'Color',[0 0 0 1],'LineWidth',1,'LineStyle','--')
        hold on
        plot(ax(aa),MLvalues.MLvalues(1,aa+1),MLvalues.MLvalues(1,1),'Marker','o',...
                'MarkerSize',8,...
                'MarkerFaceColor',[1 0 0],...
                'MarkerEdgeColor',[1 0 0],...
                'Color',[1 0 0],...
                'LineStyle','none')
        hold on
        xlabel(ax(aa),xLabels{aa},'FontSize',14,'Interpreter','latex');
        ylabel(ax(aa),'$-$ Profile likelihood','FontSize',14,'Interpreter','latex');
        text(0.075,0.9,figlabs{aa},...
                'Units','Normalized',...
                'HorizontalAlignment','left',...
                'FontSize',fontsize,...
                'FontWeight','Normal'); 
        
    end
    tiledplot.TileSpacing = 'compact';
    tiledplot.Padding = 'compact';

    savefig(strcat('./figures/practical_identifiability.fig'));
    exportgraphics(gcf,strcat('./figures/practical_identifiability.png'));
    exportgraphics(gcf,strcat('./figures/practical_identifiability.eps'),'ContentType','vector');


end