function [] = plot_mcmc_kinetics()

    %%% command:
    %%% plot_mcmc_kinetics()  
   
    clc;
    close all
    set(0,'DefaultFigureVisible','on');
    
    %%% fixed parameters %%%
    DRif            = 10;
    k_r             = 0.049;
    k_pxrdeg        = 0.011;
    k_mRNAcyp3a4deg = 0.044;
    k_mRNAcyp2c9deg = 0.036;
    k_mRNAcyp2b6deg = 0.034;
        
    vars = [2,3,4];
    metadata = load('../datamat.mat');
    
    data    = {{[]};...
               metadata.dataCYP3A4_des;...
               metadata.dataCYP2C9_des;...
               metadata.dataCYP2B6_des
               };
    stdev   = {{[]};...
               metadata.std_dataCYP3A4_des;...
               metadata.std_dataCYP2C9_des;...
               metadata.std_dataCYP2B6_des};
    time    = {{[]};...
               metadata.time_des;...
               metadata.time_des;...
               metadata.time_des};
    tmax    = metadata.tmax_des;
    T       = metadata.t_des;        
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
               ' mRNA fold change'
               };
    Text  = {'PXR';...
             'CYP3A4';...
             'CYP2C9';...
             'CYP2B6'
             };
    mcolors = {[0 0 0];...
               [0 0 0];...
               [0 0 0];...
               [0 0 0]};

    data = vertcat(data{:});
    stdev = vertcat(stdev{:});
    
    if ~exist('./figures', 'dir')
        mkdir('./figures')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    outputvec = [];
    for jj = 1:5
        chains = load(strcat('./chains/chains_',num2str(jj),'.mat'));
        chains = chains.chains(:,:);
        for nn = 1:size(chains,1)
            if mod(nn,5)==0
                outputvec = [outputvec; model(chains(nn,:))];
            end
        end
    end  
    %%% PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    figlabs = {'(C1)','(C2)','(C3)','(C4)'};
    fontsize    = 14;
    tiledplot = tiledlayout(1,4,'TileSpacing','Compact');
    set(gcf, 'Position',  [300, 100, 1500, 300]);
    for aa = 1:4
        ax(aa) = nexttile(aa);            
        set(ax(aa),...
            'box','on',...
            'XLim',[-0.05*max(timeXaxis{aa}) 1.05*max(timeXaxis{aa})],...
            'XTick',timeXaxis{aa},...
            'XTickLabel',timeXaxis{aa},...
            'XTickLabelRotation',0,...
            'YLim',YLims{aa},...
            'YTick',YTick{aa},...
            'FontSize',fontsize);
        set(gca,'TickLength',[0.025, 0.01])
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        hold on;
        %%% shaded area %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        out_aa = [];
        for kk = 1:size(outputvec(:,aa),1)
            solution = outputvec(kk,aa);
            out_aa = [out_aa; solution{1}];
        end
        outputvecmin = min(out_aa,[],1);
        outputvecmax = max(out_aa,[],1);
        patch(ax(aa),[T,fliplr(T)],[outputvecmin,fliplr(outputvecmax)],mcolors{aa},'FaceAlpha',0.15,'EdgeColor',mcolors{aa},'EdgeAlpha',0.25);
        plot(ax(aa),T,outputvecmin,'Color',[mcolors{aa} 0.1],'LineWidth',0.1,'LineStyle','-')
        plot(ax(aa),T,outputvecmax,'Color',[mcolors{aa} 0.1],'LineWidth',0.1,'LineStyle','-')
        %%% maximum likelihood %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        MLpars = readmatrix('maxLikValues.txt');
        outputmaxlik = model(MLpars(2:end,2));
        plot(ax(aa),T,outputmaxlik{aa},'Color',mcolors{aa},'LineWidth',1,'LineStyle','-')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        hold on;
        if aa == 1
            plot(ax(aa),T,zeros(length(T)),'Color',[0 0 0],'LineWidth',1,'LineStyle','--')
            hold on;
        else
            plot(ax(aa),T,ones(length(T)),'Color',[0 0 0],'LineWidth',1,'LineStyle','--')
            hold on;
        end
        if ismember(aa,vars)
            errorbar(time{aa},nanmean(data{aa},1),stdev{aa},... 
                'Marker','o', ...
                'MarkerSize',7,...
                'MarkerEdgeColor',[0 0 0],...
                'MarkerFaceColor',[0.75 0.75 0.75],...
                'LineStyle','none',...
                'LineWidth',1,...
                'Color',[0 0 0],...
                'CapSize',7);
        end
        title(strcat('\rm',Text{aa},ylabels{aa}),'FontSize',fontsize);
        xlabel(ax(aa),'Time (hours)','FontSize',14);
        text(0.1,0.9,figlabs{aa},...
                'Units','Normalized',...
                'HorizontalAlignment','center',...
                'FontSize',fontsize,...
                'FontWeight','Normal'); 
        
    end
    tiledplot.TileSpacing = 'compact';
    tiledplot.Padding = 'compact';


    savefig(strcat('figures/kinetics_varying_transcription.fig'));
    exportgraphics(gcf,strcat('figures/kinetics_varying_transcription.png'));
    exportgraphics(gcf,strcat('../../LaTeX/figures/kinetics_varying_transcription.eps'),'ContentType','vector');
             

    %%% helper functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% output of the model %%%
    function [output] = model(pars) 
        solution_des = ode23s(@odeDES,[0 tmax],...
                       [0 1 1 1],...
                       [],...
                       pars);
    
        for ii = 1:4
            output{ii} = deval(solution_des,T,ii);
        end
    end


    %%% ODE system - DESRIF %%%
    function [dxdt] = odeDES(t,x,pars)        
        dxdt = zeros(4,1);
        
        dxdt(1) = pars(1)*DRif*(1 - x(1))*exp(-k_r*t) - k_pxrdeg*x(1);                              % pxr
        dxdt(2) = pars(2)*x(1) + k_mRNAcyp3a4deg*(1 - x(2));                               % mRNA CYP3A4 
        dxdt(3) = pars(3)*x(1) + k_mRNAcyp2c9deg*(1 - x(3));                               % mRNA CYP2C9
        dxdt(4) = pars(4)*x(1) + k_mRNAcyp2b6deg*(1 - x(4));                               % mRNA CYP2B6
    end

end