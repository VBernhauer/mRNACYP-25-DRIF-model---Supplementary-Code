function [] = plot_EC50_exposure()

    %%% command:
    %%% plot_EC50_exposure()

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
    T           = [0.5:0.01:1,1.1:0.1:12,13:1:120]; 
    T_drif      = [5:0.1:12,13:1:120]; 

    logC        = -3:0.01:2;
    C           = 10 .^ logC;

    if ~exist('./figures', 'dir')
        mkdir('./figures')
    end

    
    EC50_rif    = [];
    EC50_rif_CImin = [];
    EC50_rif_CImax = [];

    EC50_drif   = [];
    EC50_drif_CImin = [];
    EC50_drif_CImax = [];

    for idx_T = 1:length(T)

        data_rif_min = load(strcat('./EC50_minmax/EC50_rif_min_T_',num2str(T(idx_T)),'.mat'));
        rif_min = data_rif_min.('data_rif_min');

        data_rif_max = load(strcat('./EC50_minmax/EC50_rif_max_T_',num2str(T(idx_T)),'.mat'));
        rif_max = data_rif_max.('data_rif_max');

        data_drif_min = load(strcat('./EC50_minmax/EC50_drif_min_T_',num2str(T(idx_T)),'.mat'));
        drif_min = data_drif_min.('data_drif_min');

        data_drif_max = load(strcat('./EC50_minmax/EC50_drif_max_T_',num2str(T(idx_T)),'.mat'));
        drif_max = data_drif_max.('data_drif_max');


        CYPmin  = 1;
        CYPmax  = 3;
        N = 2;
        logEC50 = 0.5;
        hill_pars = [CYPmax, N, logEC50];
        opts = optimset('Display','off');

        %%% RIF %%%
        var50_rif = [];
        for ii = 1:length(C)
            output_rif = model_rif(C(ii),T(idx_T));
            var50_rif  = [var50_rif, output_rif];
        end

        [EC50_rif_pars, resnorm] = lsqcurvefit(@hill_cyp,hill_pars,logC,var50_rif,[],[],opts);
        EC50_rif    = [EC50_rif, EC50_rif_pars(end)];

        cyp3a4_rif  = csapi(logC,hill_cyp(EC50_rif_pars,logC),EC50_rif_pars(end));        

        interp_cyp3a4_rif_min = csapi(rif_max,logC,cyp3a4_rif);
        interp_cyp3a4_rif_max = csapi(rif_min,logC,cyp3a4_rif);

        EC50_rif_CImin  = [EC50_rif_CImin, interp_cyp3a4_rif_min];
        EC50_rif_CImax  = [EC50_rif_CImax, interp_cyp3a4_rif_max];

        if T(idx_T) >= T_drif(1)

            %%% 25-DRIF %%%
            var50_drif = [];
            for ii = 1:length(C)
                output_drif = model_drif(pars, C(ii), T(idx_T));
                var50_drif  = [var50_drif, output_drif];
            end   

            [EC50_drif_pars, resnorm]   = lsqcurvefit(@hill_cyp,hill_pars,logC,var50_drif,[],[],opts);
            EC50_drif                   = [EC50_drif, EC50_drif_pars(end)];

            cyp3a4_drif = csapi(logC,hill_cyp(EC50_drif_pars,logC),EC50_drif_pars(end));

            interp_cyp3a4_drif_min = csapi(drif_max,logC,cyp3a4_drif);
            interp_cyp3a4_drif_max = csapi(drif_min,logC,cyp3a4_drif);

            EC50_drif_CImin  = [EC50_drif_CImin, interp_cyp3a4_drif_min];
            EC50_drif_CImax  = [EC50_drif_CImax, interp_cyp3a4_drif_max];
        end

        disp(strcat({'Time '},num2str(T(idx_T)),{' done.'}));

    end


    fontsize    = 14;

    figure();
    ax = gca;
    set(ax,...
        'box','on',...
        'XLim',[-0.05*120 1.05*120],...
        'XTick',[0,12,24,48,72,120],...
        'YLim',[-1-0.25 3+0.25],...
        'YTick',[-1,0,1,2,3],...
        'FontSize',fontsize);

    xlabel('Exposure time (h)','FontSize',fontsize);
    ylabel('Log_{10} EC_{50}','FontSize',fontsize);
    set(gca,'TickLength',[0.015, 0.015])
    hold on;  

    patch([T,fliplr(T)],[EC50_rif_CImin,fliplr(EC50_rif_CImax)],[0 0 1],'FaceAlpha',0.15,'EdgeColor',[0 0 1],'EdgeAlpha',0.25);

    rif=plot(T, EC50_rif,...
          'Color',[0 0 1], ...
          'LineWidth',1.5,...
          'DisplayName','RIF');

    patch([T_drif,fliplr(T_drif)],[EC50_drif_CImin,fliplr(EC50_drif_CImax)],[0 1 1],'FaceAlpha',0.15,'EdgeColor',[0 1 1],'EdgeAlpha',0.25);

    drif=plot(T_drif, EC50_drif,...
          'Color',[0 1 1], ...
          'LineWidth',1.5,...
          'DisplayName','25-DRIF');

    leg = legend(ax, [rif,drif]);
    pos = get(leg,'Position'); 
    pos(1) = pos(1) + 0.01;
    pos(2) = pos(2) - 0.01;
    set(leg,'Position',pos);
    set(leg,'IconColumnWidth',10);
    set(leg,'FontSize',fontsize); 

    savefig(strcat('./figures/EC50_exposure.fig'));
    exportgraphics(gcf,strcat('./figures/EC50_exposure.png'));
    exportgraphics(gcf,strcat('../LaTeX/figures/EC50_exposure.eps'),'ContentType','vector');


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