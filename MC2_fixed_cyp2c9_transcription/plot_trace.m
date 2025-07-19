function [] = plot_trace()

    %%% command:
    %%% plot_trace()
   
    clc;
    close all
    set(0,'DefaultFigureVisible','on');
    
    if ~exist('./figures', 'dir')
        mkdir('./figures')
    end
    
    a{1} = '$\mathrm{Log}_{10}$ $k_\mathrm{pxr,max}$';
    a{2} = '$\mathrm{Log}_{10}$ $k_{{mRNA}_\mathrm{cyp3a4}^\mathrm{fold}}$';
    a{3} = '$\mathrm{Log}_{10}$ $k_{{mRNA}_\mathrm{cyp2b6}^\mathrm{fold}}$';

    labs = {a{1},...
            a{2},...
            a{3}};

    col = {[1 0 0]; [0 1 0]; [0 0 1]; [1 0 1]; [1 1 0]};
        
    tiledplot = tiledlayout(1,3,'TileSpacing','Compact');
    set(gcf, 'Position',  [500, 100, 1200, 300]);
    for jj = 1:5
        chains = load(strcat('./chains/chains_',num2str(jj),'.mat'));
        chains = chains.chains(:,:);
        mm = 0;
        for aa = 1:length(labs)
            mm = mm+1;
            ax(mm) = nexttile(mm);
            plot(ax(mm),1:1:length(chains(:,aa)),log10(chains(:,aa)),'Color',[col{jj} 1],'LineWidth',0.1,'LineStyle','-');
            hold on;
            ylabel(ax(mm),labs{aa},'FontSize',12,'Interpreter','latex');
            xlabel(ax(mm),'Iteration','FontSize',12);
            set(ax(mm),'xlim',[1 length(chains(:,aa))]);
        end
    end   
    tiledplot.TileSpacing = 'compact';
    tiledplot.Padding = 'compact';

    savefig(strcat('./figures/trace_fixed_cyp2c9_transcription.fig'));
    exportgraphics(gcf,strcat('./figures/trace_fixed_cyp2c9_transcription.png'));
    exportgraphics(gcf,strcat('./figures/trace_fixed_cyp2c9_transcription.eps'),'ContentType','vector');
    
end