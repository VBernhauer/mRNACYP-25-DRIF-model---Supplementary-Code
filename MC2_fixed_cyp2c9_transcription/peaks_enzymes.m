function [] = peaks_enzymes()

    %%% command:
    %%% peaks_enzymes()  

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

    tmax = metadata.tmax_des;
    timeXaxis = {[0,12,24,48,72,120];...
                 [0,12,24,48,72,120];...
                 [0,12,24,48,72,120];...
                 [0,12,24,48,72,120]
                 };
    YLims = {[-0.1 1.1];...
             [-0.35 4.35];...
             [-0.35 4.35];...
             [-0.35 4.35]};    
    YTick = {0:0.2:1;...
             0:1:4;...
             0:1:4;...
             0:1:4}; 
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
    figlabs = {'(A)','(B)','(C)','(D)'};
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
            output = peak(pars,Dose(i));
            disp(['Var: ', Text{mm}, ' Dose: ', num2str(Dose(i)), ' Time: ', num2str(output(mm,1)), ' Peak: ', num2str(output(mm,2))]);
            plot(ax(mm),output(mm,1),output(mm,2),...
                'Marker','o',...
                'MarkerSize',5,...
                'MarkerFaceColor',colors(i,:),...
                'MarkerEdgeColor',colors(i,:),...
                'Color',colors(i,:),...
                'LineStyle','none')
            hold on;
            title(strcat('\rm',Text{mm},ylabels{mm}),'FontSize',12);
            text(0.9,0.9,figlabs{mm},...
                'Units','Normalized',...
                'HorizontalAlignment','center',...
                'FontSize',fontsize,...
                'FontWeight','Normal');    
            hold on;
        end
        if mm == 1
            plot(ax(mm),0:1:120,zeros(length(0:1:120)),'Color',[0 0 0],'LineWidth',1,'LineStyle','--')
            hold on;
        else
            plot(ax(mm),0:1:120,ones(length(0:1:120)),'Color',[0 0 0],'LineWidth',1,'LineStyle','--')
            hold on;
        end
    end
    leg = legend(ax(4),Legend{:,:},'Location','EastOutside','FontSize',fontsize,'Orientation','Vertical'); 
    title(leg,'25-DRIF concentration','FontWeight','Normal')
    xlabel(tiledplot,'Time (hours)','FontSize',fontsize);    
    tiledplot.TileSpacing = 'compact';
    tiledplot.Padding = 'compact';


    %%% helper functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    %%% output of the model %%%
    function [output] = peak(pars,DRif) 
        solution = ode23s(@ode,[0 tmax],...
                               [0 1 1 1],...
                               [],...
                               pars,...
                               DRif);
        output = [];
        for kk = 1:4
            time = solution.x';
            sol  = solution.y';
            [solpeak,id] = max(sol(:,kk));
            output = [output; [time(id), solpeak]];
        end
    end

    %%% ODE system %%%
    function [dxdt] = ode(t,x,pars,DRif)        
        dxdt = zeros(4,1);

        dxdt(1) = pars(1)*DRif*(1 - x(1))*exp(-k_r*t) - k_pxrdeg*x(1);                     % pxr
        dxdt(2) = pars(2)*x(1) + k_mRNAcyp3a4deg*(1 - x(2));                               % mRNA CYP3A4 
        dxdt(3) = k_mRNAcyp2c9fold*x(1) + k_mRNAcyp2c9deg*(1 - x(3));                      % mRNA CYP2C9
        dxdt(4) = pars(3)*x(1) + k_mRNAcyp2b6deg*(1 - x(4));                               % mRNA CYP2B6
    end

end