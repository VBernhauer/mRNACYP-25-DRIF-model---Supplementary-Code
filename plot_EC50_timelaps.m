function [] = plot_EC50_timelaps()

    %%% command:
    %%% plot_EC50_timelaps()

    clc;
    close all
    set(0,'DefaultFigureVisible','on');

    name = 'cyp3a4';
    ylabelfig = 'CYP3A4 mRNA fold change';
        

    %%% fixed parameters 
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
    
    tmax        = metadata.tmax_des;
    T           = [0.01:0.01:1,1.1:0.1:12,13:1:120]; 
    T_specific  = [4,8,12,16,24,48,72,96,120];

    idx_T_specific_vec = [];
    for idx_T_specific = 1:length(T_specific)
        for idx_T = 1:length(T)
            if T(idx_T) == T_specific(idx_T_specific)
                idx_T_specific_vec = [idx_T_specific_vec,idx_T];
            end
        end
    end
    logC        = -3:0.01:2;
    C           = 10 .^ logC;

    if ~exist('./figures', 'dir')
        mkdir('./figures')
    end


    ax = figure();
    set(gcf, 'Position',  [300, 0, 1100, 900]);
    fontsize    = 12;
    ec_fontsize = 7;

    varmin = 1;
    varmax = 3;

    rifname = 'RIF, EC_{50}=';
    drifname = '25-DRIF, EC_{50}=';

    for i = 1:length(T_specific)

        data_rif_min = load(strcat('./EC50_minmax/EC50_rif_min_T_',num2str(T_specific(i)),'.mat'));
        data_rif_c_min = data_rif_min.('data_rif_min');

        data_rif_max = load(strcat('./EC50_minmax/EC50_rif_max_T_',num2str(T_specific(i)),'.mat'));
        data_rif_c_max = data_rif_max.('data_rif_max');

        data_drif_min = load(strcat('./EC50_minmax/EC50_drif_min_T_',num2str(T_specific(i)),'.mat'));
        data_drif_c_min = data_drif_min.('data_drif_min');

        data_drif_max = load(strcat('./EC50_minmax/EC50_drif_max_T_',num2str(T_specific(i)),'.mat'));
        data_drif_c_max = data_drif_max.('data_drif_max');

        EC50_drif   = [];
        EC50_rif    = [];

        CYPmin  = 1;
        CYPmax  = 3;
        N = 2;
        logEC50 = 0.5;
        hill_pars = [CYPmax, N, logEC50];
        opts = optimset('Display','off');

        var50_rif = [];
        for ii = 1:length(C)
            output_rif = model_rif(C(ii),T_specific(i));
            var50_rif  = [var50_rif, output_rif];
        end
        
        [EC50_rif_pars, resnorm] = lsqcurvefit(@hill_cyp,hill_pars,logC,var50_rif,[],[],opts);
        EC50_rif    = [EC50_rif; T_specific(i), EC50_rif_pars(end)];

        var50_drif = [];
        for ii = 1:length(C)
            output_drif = model_drif(pars, C(ii), T_specific(i));
            var50_drif  = [var50_drif, output_drif];
        end
        [EC50_drif_pars, resnorm] = lsqcurvefit(@hill_cyp,hill_pars,logC,var50_drif,[],[],opts);
        EC50_drif   = [EC50_drif; T_specific(i), EC50_drif_pars(end)];   

        interp_rif  = csapi(logC,hill_cyp(EC50_rif_pars,logC));
        interp_drif = csapi(logC,hill_cyp(EC50_drif_pars,logC));


        ax(i) = subplot(3,3,i);%nexttile(i);
        set(ax(i),...
            'box','on',...
            'XLim',[logC(1)-0.25 logC(end)+0.25],...
            'XTick',[-3,-2,-1,0,1,2],...
            'YLim',[varmin-0.1 varmax+0.1],...
            'FontSize',fontsize);

        xlabel('Log_{10} concentration (µM)','FontSize',fontsize);
        ylabel(ylabelfig,'FontSize',fontsize);
        set(gca,'TickLength',[0.015, 0.015])
        hold on;  

        patch([logC,fliplr(logC)],[data_rif_c_min,fliplr(data_rif_c_max)],[0 0 1],'FaceAlpha',0.15,'EdgeColor',[0 0 1],'EdgeAlpha',0.25);

        rif=plot(logC, var50_rif,...
              'Color',[0 0 1], ...
              'LineWidth',1.5);

        patch([logC,fliplr(logC)],[data_drif_c_min,fliplr(data_drif_c_max)],[0 1 1],'FaceAlpha',0.15,'EdgeColor',[0 1 1],'EdgeAlpha',0.25);

        drif=plot(logC, var50_drif,...
              'Color',[0 1 1], ...
              'LineWidth',1.5);

        line1=plot([logC(1) logC(1)], [varmin varmax],...
              'Color',[0 0 0], ...
              'LineWidth',1.5,...
              'LineStyle','-');
        line2=plot([logC(1) logC(end)], [varmax varmax],...
              'Color',[0 0 0], ...
              'LineWidth',1.5,...
              'LineStyle','-');
        line3=plot([logC(end) logC(end)], [varmin varmax],...
              'Color',[0 0 0], ...
              'LineWidth',1.5,...
              'LineStyle','-');
        line4=plot([logC(1) logC(end)], [varmin varmin],...
              'Color',[0 0 0], ...
              'LineWidth',1.5,...
              'LineStyle','-');

        ec50_rif_ver=plot([EC50_rif_pars(end) EC50_rif_pars(end)], [1 fnval(interp_rif,EC50_rif_pars(end))],...
              'Color',[0 0 1], ...
              'LineWidth',1.5,...
              'LineStyle',':');
        ec50_rif_hor=plot([logC(1) EC50_rif_pars(end)], [fnval(interp_rif,EC50_rif_pars(end)) fnval(interp_rif,EC50_rif_pars(end))],...
              'Color',[0 0 1], ...
              'LineWidth',1.5,...
              'LineStyle',':');

        ec50_drif_ver=plot([EC50_drif_pars(end) EC50_drif_pars(end)], [1 fnval(interp_drif,EC50_drif_pars(end))],...
              'Color',[0 1 1], ...
              'LineWidth',1.0,...
              'LineStyle',':');
        ec50_drif_hor=plot([logC(1) EC50_drif_pars(end)], [fnval(interp_drif,EC50_drif_pars(end)) fnval(interp_drif,EC50_drif_pars(end))],...
              'Color',[0 1 1], ...
              'LineWidth',1.5,...
              'LineStyle',':');
        
        point_ec50_rif=plot(EC50_rif_pars(end), fnval(interp_rif,EC50_rif_pars(end)),...
              'o',...
              'MarkerFaceColor',[0 0 1], ...
              'MarkerEdgeColor',[0 0 0],...
              'MarkerSize',8,...
              'LineWidth',0.1,...
              'DisplayName',[rifname,' ',num2str(10^EC50_rif_pars(end),'%4.2f'),' µM']);
        point_ec50_drif=plot(EC50_drif_pars(end), fnval(interp_drif,EC50_drif_pars(end)),...
              'o',...
              'MarkerFaceColor',[0 1 1], ...
              'MarkerSize',8,...
              'MarkerEdgeColor',[0 0 0],...
              'LineWidth',0.1,...
              'DisplayName',[drifname,' ',num2str(10^EC50_drif_pars(end),'%4.2f'),' µM']);
      
        leg = legend(ax(i), [point_ec50_rif,point_ec50_drif]);
        set(leg,'units','normalized');
        pos = get(leg,'Position');        
        pos(2) = pos(2) + 0.005;
        if i > 5
            pos(1) = pos(1) - 0.047;
        else
            pos(1) = pos(1) - 0.041;
        end
        
        set(leg,'Position',pos);
        set(leg,'FontSize',ec_fontsize); 
        
        title(strcat(num2str(T_specific(i)),{' h post-treatment'}),...
            'FontWeight','Normal',...
            'FontSize',fontsize);

        disp(strcat({'Image '},num2str(i),{' done.'}));
    end

    savefig(strcat('./figures/EC50_',name,'.fig'));
    exportgraphics(gcf,strcat('./figures/EC50_',name,'.png'));
    exportgraphics(gcf,strcat('./figures/EC50_',name,'.eps'),'ContentType','vector');


    %%% helper functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [CYP] = hill_cyp(hill_pars,logC)
        CYP = CYPmin + (hill_pars(1) - CYPmin) ./ (1 + 10 .^ (hill_pars(2) * (hill_pars(3) - logC)));
    end

    function [output] = model_drif(pars,C,T) 
        solution_drif = ode23s(@ode_drif,[0 tmax],...
                               [0 1 1 1],...
                               [],...
                               pars,...
                               C);
        output = deval(solution_drif,T,2);
 
    end

    function [output] = model_rif(C,T) 
        solution_rif = ode23s(@ode_rif,[0 tmax],...
                               [0 1 1 1],...
                               [],...
                               C);
        output = deval(solution_rif,T,2);
 
    end

    function [dxdt] = ode_drif(t,x,pars,C)        
        dxdt = zeros(4,1);

        dxdt(1) = pars(1)*C*(1 - x(1))*exp(-k_r*t) - k_pxrdeg*x(1);                     % pxr
        dxdt(2) = pars(2)*x(1) + k_mRNAcyp3a4deg*(1 - x(2));                               % mRNA CYP3A4 
        dxdt(3) = k_mRNAcyp2c9fold*x(1) + k_mRNAcyp2c9deg*(1 - x(3));                      % mRNA CYP2C9
        dxdt(4) = pars(3)*x(1) + k_mRNAcyp2b6deg*(1 - x(4));                               % mRNA CYP2B6

    end

    function [dxdt] = ode_rif(t,x,C)        
        dxdt = zeros(4,1);

        dxdt(1) = k_pxrmax*C*(1 - x(1))*exp(-k_r*t) - k_pxrdeg*x(1);                     % pxr
        dxdt(2) = k_mRNAcyp3a4fold*x(1) + k_mRNAcyp3a4deg*(1 - x(2));                      % mRNA CYP3A4 
        dxdt(3) = k_mRNAcyp2c9fold*x(1) + k_mRNAcyp2c9deg*(1 - x(3));                      % mRNA CYP2C9
        dxdt(4) = k_mRNAcyp2b6fold*x(1) + k_mRNAcyp2b6deg*(1 - x(4));                      % mRNA CYP2B6
    end


end