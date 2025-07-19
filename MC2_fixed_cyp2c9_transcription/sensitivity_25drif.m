function [] = sensitivity_25drif()

    %%% command:
    %%% sensitivity_25drif()  
   
    clc;
    close all
    set(0,'DefaultFigureVisible','on');
        
    metadata = load('../datamat.mat');

    %%% estimated parameters
    pars = readtable('maxLikValues.txt');
    pars = pars.Var2(2:end);

    %%% fixed parameters
    k_r                 = 0.049;
    k_pxrdeg            = 0.011;
    k_mRNAcyp3a4deg     = 0.044;
    k_mRNAcyp2c9fold    = 0.040;
    k_mRNAcyp2c9deg     = 0.036;
    k_mRNAcyp2b6deg     = 0.034;

    Dose = [1,5,10,20,50,100,200];
      
    T = metadata.t_des;
    tmax = metadata.tmax_des;
    timeXaxis = {[0,12,24,48,72,120];...
                 [0,12,24,48,72,120];...
                 [0,12,24,48,72,120];...
                 [0,12,24,48,72,120]
                 };
    YLims = {[-0.1 1.1];...
             [-0.5 5.5];...
             [-0.5 5.5];...
             [-0.5 5.5]};    
    YTick = {0:0.2:1;...
             0:1:5;...
             0:1:5;...
             0:1:5}; 
    ylabels = {' activation';...
               ' mRNA fold change';...
               ' mRNA fold change';...
               ' mRNA fold change'};
    Text  = {'PXR';...
             'CYP3A4';...
             'CYP2C9';...
             'CYP2B6'};
    
    if ~exist('./figures', 'dir')
        mkdir('./figures')
    end
    
    
    %%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    colors = jet(length(Dose));
    fontsize = 14;
    for ii = 1:length(Dose)
        Legend{ii} = strcat(num2str(Dose(ii)),' µM');
    end
    figlabs = {'(B1)','(B2)','(B3)','(B4)'};
    tiledplot = tiledlayout(1,4,'TileSpacing','Compact');
    set(gcf, 'Position',  [300, 100, 1500, 300]);
    mm = 0;
    for aa = 1:4
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
            xlabel(ax(aa),'Time (hours)','FontSize',fontsize);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            hold on;
            for i = 1:length(Dose)
                output = model(pars,Dose(i));
                plot(ax(mm),T,output{mm},'Color',colors(i,:),'LineWidth',1,'LineStyle','-')
                hold on;
            end
            if mm == 1
                plot(ax(mm),T,zeros(length(T)),'Color',[0 0 0],'LineWidth',1,'LineStyle','--')
                hold on;
            else
                plot(ax(mm),T,ones(length(T)),'Color',[0 0 0],'LineWidth',1,'LineStyle','--')
                hold on;
            end
            title(strcat('\rm',Text{mm},ylabels{mm}),'FontSize',fontsize);
            text(0.9,0.9,figlabs{mm},...
                'Units','Normalized',...
                'HorizontalAlignment','center',...
                'FontSize',fontsize,...
                'FontWeight','Normal');    
            hold on;
    end
    leg = legend(ax(4),Legend{:,:},'Location','EastOutside','FontSize',fontsize,'Orientation','Vertical','NumColumns',1); 
    title(leg,'25-DRIF (concentration)','FontWeight','Normal')
    tiledplot.TileSpacing = 'compact';
    tiledplot.Padding = 'compact';


    savefig(strcat('./figures/concentration_sensitivity.fig'));
    exportgraphics(gcf,'./figures/concentration_sensitivity.png');
    exportgraphics(gcf,'./figures/concentration_sensitivity.eps','ContentType','vector');
             

    
    %%% helper functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% output of the model %%%
    function [output] = model(pars,DRif) 
        solution = ode23s(@ode,[0 tmax],...
                               [0 1 1 1],...
                               [],...
                               pars,...
                               DRif);
        for kk = 1:4
            output{kk} = deval(solution,T,kk);
            kk = kk + 1;
        end
    end

 %%% ODE system - DESRIF %%%
    function [dxdt] = ode(t,x,pars,DRif)        
        dxdt = zeros(4,1);

        dxdt(1) = pars(1)*DRif*(1 - x(1))*exp(-k_r*t) - k_pxrdeg*x(1);                     % pxr
        dxdt(2) = pars(2)*x(1) + k_mRNAcyp3a4deg*(1 - x(2));                               % mRNA CYP3A4 
        dxdt(3) = k_mRNAcyp2c9fold*x(1) + k_mRNAcyp2c9deg*(1 - x(3));                      % mRNA CYP2C9
        dxdt(4) = pars(3)*x(1) + k_mRNAcyp2b6deg*(1 - x(4));                               % mRNA CYP2B6
    end
       
end