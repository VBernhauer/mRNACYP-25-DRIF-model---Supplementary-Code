function [] = plot_data()

    %%% command:
    %%% plot_data()
   
    clc;
    close all;
    clear all;
    set(0,'DefaultFigureVisible','on');
        
    % load data and parameters
    metadata = load('datamat.mat');
    
    %%% CYPs %%%
        
    data    = {metadata.dataCYP3A4_des;...
               metadata.dataCYP2C9_des;...
               metadata.dataCYP2B6_des};
    time    = {metadata.time_des;...
               metadata.time_des;...
               metadata.time_des};
    timeXaxis = {[0,12,24,48,72,120];...
                 [0,12,24,48,72,120];...
                 [0,12,24,48,72,120]
                 };
    YLims = {[-0.35 4.35];...
             [-0.35 4.35];...
             [-0.35 4.35]};    
    YTick = {0:1:4;...
             0:1:4;...
             0:1:4};  
    ylabels = {' mRNA fold change';...
               ' mRNA fold change';...
               ' mRNA fold change'
               };
    Text  = {'CYP3A4';'CYP2C9';'CYP2B6'};
    
    if ~exist('./figures', 'dir')
        mkdir('./figures')
    end

    %%%% Each donor in different color %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    facecol     = {'red','blue','green'};
    names       = {'Donor 1','Donor 2','Donor 3'};
    markers     = {'o','diamond','square'};
    figlabs     = {'(A)','(B)','(C)'};
    markersize  = 50;
    markerfa    = {0.5,0.5,0.5};
    mcolors     = [0 0 0];
    fontsize    = 14;
    
    figure(1);
    clc;
    tiledplot = tiledlayout(1,length(data));
    set(gcf, 'Position',  [500, 300, 1200, 350]);
    clc;
    mm = 0;
    for aa = 1:length(data)
        mm = mm+1;
        ax(mm) = nexttile(mm);
        set(ax(mm),...
            'box','on',...
            'XLim',[-0.05*max(timeXaxis{mm}) 1.05*max(timeXaxis{mm})],...
            'XTick',timeXaxis{mm},...
            'XTickLabel',timeXaxis{mm},...
            'XTickLabelRotation',0,...
            'YLim',YLims{mm},...
            'YTick',YTick{mm},...
            'FontSize',fontsize);
        set(gca,'TickLength',[0.025, 0.01])
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        hold on;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for kk = 1:3
            scatter(time{mm},data{mm}(kk,:),markersize,... 
                'Marker',markers{kk}, ...
                'MarkerEdgeColor',mcolors,...
                'MarkerFaceColor',facecol{kk},...
                'MarkerFaceAlpha',markerfa{kk});
        end
        hold on;
        plot(ax(mm),0:10:120,ones(length(0:10:120)),'Color',[0 0 0],'LineWidth',1,'LineStyle','--')
        hold on;
        xlabel(ax(aa),'Time (hours)','FontSize',14);
        title(strcat('\rm',Text{mm},ylabels{mm}),'FontSize',fontsize);
        text(0.1,0.9,figlabs{mm},...
                'Units','Normalized',...
                'HorizontalAlignment','center',...
                'FontSize',fontsize,...
                'FontWeight','Normal');                    
    end
    leg = legend(ax(1),names,'Location','SouthOutside','FontSize',fontsize,'Orientation','Horizontal');
    leg.Layout.Tile = 'North';
    
    tiledplot.TileSpacing = 'compact';
    tiledplot.Padding = 'compact';

    savefig(strcat('figures/data.fig'));
    exportgraphics(gcf,'figures/data.png');
    exportgraphics(gcf,'../LaTeX/figures/data.eps','ContentType','vector');

    end