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
    a{3} = '$\mathrm{Log}_{10}$ $k_{{mRNA}_\mathrm{cyp2c9}^\mathrm{fold}}$';
    a{4} = '$\mathrm{Log}_{10}$ $k_{{mRNA}_\mathrm{cyp2b6}^\mathrm{fold}}$';

    labs = {a{1},...
            a{2},...
            a{3},...
            a{4}};

    col = {[1 0 0]; [0 1 0]; [0 0 1]; [1 0 1]; [1 1 0]};
        
    tiledplot = tiledlayout(2,2,'TileSpacing','Compact');
    set(gcf, 'Position',  [500, 100, 800, 500]);
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

    savefig('figures/trace_varying_transcription.fig');
    exportgraphics(gcf,'figures/trace_varying_transcription.png');
    exportgraphics(gcf,'../../LaTeX/figures/trace_varying_transcription.eps','ContentType','vector');
    
end