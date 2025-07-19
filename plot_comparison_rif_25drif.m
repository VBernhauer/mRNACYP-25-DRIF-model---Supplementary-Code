function [] = plot_comparison_rif_25drif()

    %%% command:
    %%% plot_comparison_rif_25drif()  
   
    clc;
    close all
    set(0,'DefaultFigureVisible','on');

    %%% fixed parameters 
    DRif                = 10;
    Rif                 = [1, 10];
    k_pxrmax            = 0.143;
    k_r                 = 0.049;
    k_pxrdeg            = 0.011;
    k_mRNAcyp3a4fold    = 0.083;
    k_mRNAcyp3a4deg     = 0.044;
    k_mRNAcyp2c9fold    = 0.040;
    k_mRNAcyp2c9deg     = 0.036;
    k_mRNAcyp2b6fold    = 0.139;
    k_mRNAcyp2b6deg     = 0.034;

    metadata = load('datamat.mat');

    %%% estimated parameters
    pars = readtable('./MC2_fixed_cyp2c9_transcription/maxLikValues.txt');
    pars = pars.Var2(2:end);
    
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
    colors = [[0 0 1];...
              [1 0 0];...
              [0 0 0]];
    fontsize    = 14;
    figlabs     = {'(A1)','(A2)','(A3)','(A4)'};
    Legend      = {'RIF (1 µM)', 'RIF (10 µM)', '25-DRIF (10 µM)'};
 
    tiledplot = tiledlayout(1,4,'TileSpacing','Compact');
    set(gcf, 'Position',  [300, 100, 1500, 300]);
    for mm = 1:4
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
            xlabel(ax(mm),'Time (hours)','FontSize',fontsize);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            hold on;
            for i = 1:3
                output = model(pars);
                plot(ax(mm),T,output{i}{mm},'Color',colors(i,:),'LineWidth',1,'LineStyle','-')
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
    title(leg,'Ligand (concentration)','FontWeight','Normal')
    tiledplot.TileSpacing   = 'compact';
    tiledplot.Padding       = 'compact';


    savefig(strcat('figures/comparison_rif_25drif.fig'));
    exportgraphics(gcf,'figures/comparison_rif_25drif.png');
    exportgraphics(gcf,'./figures/comparison_rif_25drif.eps','ContentType','vector');

    
    %%% helper functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% output of the model %%%
    function [output] = model(pars) 
        solution_rif_1uM = ode23s(@ode_rif,[0 tmax],...
                               [0 1 1 1],...
                               [],...
                               Rif(1));
        solution_rif_10uM = ode23s(@ode_rif,[0 tmax],...
                               [0 1 1 1],...
                               [],...
                               Rif(2));
        solution_drif_10uM = ode23s(@ode_drif,[0 tmax],...
                               [0 1 1 1],...
                               [],...
                               pars,...
                               DRif);
        for kk = 1:4
            output_rif_1uM{kk} = deval(solution_rif_1uM,T,kk);
            kk = kk + 1;
        end
        for kk = 1:4
            output_rif_10uM{kk} = deval(solution_rif_10uM,T,kk);
            kk = kk + 1;
        end
        for kk = 1:4
            output_drif{kk} = deval(solution_drif_10uM,T,kk);
            kk = kk + 1;
        end

        output = {output_rif_1uM, output_rif_10uM, output_drif};
    end

 %%% ODE system - DESRIF %%%
    function [dxdt] = ode_drif(t,x,pars,DRif)        
        dxdt = zeros(4,1);

        dxdt(1) = pars(1)*DRif*(1 - x(1))*exp(-k_r*t) - k_pxrdeg*x(1);                     % pxr
        dxdt(2) = pars(2)*x(1) + k_mRNAcyp3a4deg*(1 - x(2));                               % mRNA CYP3A4 
        dxdt(3) = k_mRNAcyp2c9fold*x(1) + k_mRNAcyp2c9deg*(1 - x(3));                      % mRNA CYP2C9
        dxdt(4) = pars(3)*x(1) + k_mRNAcyp2b6deg*(1 - x(4));                               % mRNA CYP2B6

    end

    function [dxdt] = ode_rif(t,x,Rif)        
        dxdt = zeros(4,1);

        dxdt(1) = k_pxrmax*Rif*(1 - x(1))*exp(-k_r*t) - k_pxrdeg*x(1);                     % pxr
        dxdt(2) = k_mRNAcyp3a4fold*x(1) + k_mRNAcyp3a4deg*(1 - x(2));                      % mRNA CYP3A4 
        dxdt(3) = k_mRNAcyp2c9fold*x(1) + k_mRNAcyp2c9deg*(1 - x(3));                      % mRNA CYP2C9
        dxdt(4) = k_mRNAcyp2b6fold*x(1) + k_mRNAcyp2b6deg*(1 - x(4));                      % mRNA CYP2B6
    end
       
end