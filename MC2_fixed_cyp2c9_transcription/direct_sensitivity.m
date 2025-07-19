function [] = direct_sensitivity()

    %%% command:
    %%% direct_sensitivity()  
   
    clc;
    close all
    set(0,'DefaultFigureVisible','on');
        
    metadata = load('../datamat.mat');

    %%% estimated parameters
    pars = readtable('maxLikValues.txt');
    pars = pars.Var2(2:end);
     
    %%% fixed parameters %%%
    DRif                = 10;
    k_r                 = 0.049;
    k_pxrdeg            = 0.011;
    k_mRNAcyp3a4deg     = 0.044;
    k_mRNAcyp2c9fold    = 0.040;
    k_mRNAcyp2c9deg     = 0.036;
    k_mRNAcyp2b6deg     = 0.034;

    parnames{1} = '$k_\mathrm{{pxr,max}}$';
    parnames{2} = '$k_\mathrm{{mRNA_{cyp3a4}^{fold}}}$';
    parnames{3} = '$k_\mathrm{{mRNA_{cyp2b6}^{fold}}}$';

    tmax = metadata.tmax_des;          
    timeXaxis = {[0,12,24,48,72,120];...
                 [];...
                 [];...
                 [0,12,24,48,72,120];...
                 [0,12,24,48,72,120];...
                 [0,12,24,48,72,120]}; 
    Text  = {'$p [pxr]_{p}/[pxr]$';...
             '';...
             '';...
             '$p {[mRNA_\mathrm{cyp3a4}^\mathrm{fold}]}_{p}/[mRNA_\mathrm{cyp3a4}^\mathrm{fold}]$';...
             '$p {[mRNA_\mathrm{cyp2c9}^\mathrm{fold}]}_{p}/[mRNA_\mathrm{cyp2c9}^\mathrm{fold}]$';...
             '$p {[mRNA_\mathrm{cyp2b6}^\mathrm{fold}]}_{p}/[mRNA_\mathrm{cyp2b6}^\mathrm{fold}]$'};
    
    if ~exist('./figures', 'dir')
        mkdir('./figures')
    end
    
    
    %%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    %colors = jet(length(pars));
    colors = [0 0 0; 1 0 0; 0 0 1]; % green for k_{pxr,max}, blue for k_{mRNA_{cyp3a4}^{fold}} , red for k_{mRNA_{cyp2b6}^{fold}}
    pos = {[40, 0.9, 0]; []; []; [70, 0.9, 0]; [70, 0.9, 0]; [70, 0.9, 0]};
    for ii = 1:length(pars)
        Legend{ii} = num2str(parnames{ii});
    end
    figlabs = {'(A)','','','(B)','(C)','(D)'};
    tiledplot = tiledlayout(2,3,'TileSpacing','Compact');
    set(gcf, 'Position',  [300, 100, 1200, 600]);
    mm = 0;
    nn = 0;
    for aa = 1:6
        mm = mm+1;
        ax(mm) = nexttile(mm);
        if ismember(mm,[2,3])
            set(ax(mm),'Visible','off');
        else
            nn = nn + 1;
            set(ax(mm),...
                'box','on',...
                'XLim',[-0.05*max(timeXaxis{mm}) 1.05*max(timeXaxis{mm})],...
                'XTick',timeXaxis{mm},...
                'XTickLabel',timeXaxis{mm},...
                'XTickLabelRotation',0,...
                'YLim',[-0.1 1.1],...
                'FontSize',14);
            set(gca,'TickLength',[0.025, 0.01])
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            hold on;
            for i = 1:length(parnames)
                [tt, output] = model(pars,i);
                plot(ax(mm),tt,output(:,nn),'Color',colors(i,:),'LineWidth',1.5,'LineStyle','-')
                hold on;
            end
            title(strcat('\rm',Text{mm}),'FontSize',14,'Interpreter','latex','Position', pos{mm});
            text(0.1,0.9,figlabs{mm},...
                'Units','Normalized',...
                'HorizontalAlignment','center',...
                'FontSize',14,...
                'FontWeight','Normal');    
            hold on;
            xlabel(ax(mm),'Time (hours)','FontSize',14);
        end
    end
    leg = legend(ax(1),Legend{:,:},...
          'Location','EastOutside',...
          'FontSize',14,...
          'Orientation','Vertical',...
          'Interpreter','latex');
    rect = [0.6, 0.7, .15, .2];
    set(leg, 'Position', rect)
    title(leg,'Parameter')  
    ylabel(tiledplot,'Normalized sensitivity coefficient','FontSize',14);
    tiledplot.TileSpacing = 'compact';
    tiledplot.Padding = 'compact';


    savefig('./figures/direct_sensitivity.fig');
    exportgraphics(gcf,'./figures/direct_sensitivity.png');
    exportgraphics(gcf,'./figures/direct_sensitivity.eps','ContentType','vector');

    
    %%% helper functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% output of the model %%%
    function [t,output] = model(pars,parsnum) 
        [t,s] = ode23s(@ode,[0 tmax],...
                       [0 1 1 1 0 0 0 0],...
                       [],...
                       pars,...
                       parsnum); 
        output_m = s(:,1:4);
        output_s = s(:,5:end);
        output = pars(parsnum) * output_s ./ output_m;
    end

    %%% Jacobian times sensitivities %%%
    function [dsdt] = ode(t,x,pars,parsnum) 

        dfdp_matrix = [(1-x(1))*DRif*exp(-k_r*t) 0 0;...
                        0 x(1) 0;...
                        0 0 0;...
                        0 0 x(1)];
        dfdp = [0;0;0;0;dfdp_matrix(:,parsnum)];
        
        dsdt = [pars(1)*DRif*(1 - x(1))*exp(-k_r*t) - k_pxrdeg*x(1);...                     
                pars(2)*x(1) + k_mRNAcyp3a4deg*(1 - x(2));...                               
                k_mRNAcyp2c9fold*x(1) + k_mRNAcyp2c9deg*(1 - x(3));...                      
                pars(3)*x(1) + k_mRNAcyp2b6deg*(1 - x(4));...                              
                -(pars(1)*DRif*exp(-k_r*t) + k_pxrdeg)*x(5);...
                pars(2)*x(5) - k_mRNAcyp3a4deg*x(6);...
                k_mRNAcyp2c9fold*x(5) - k_mRNAcyp2c9deg*x(7);...
                pars(3)*x(5) - k_mRNAcyp2b6deg*x(8)] + dfdp;       
    end

end