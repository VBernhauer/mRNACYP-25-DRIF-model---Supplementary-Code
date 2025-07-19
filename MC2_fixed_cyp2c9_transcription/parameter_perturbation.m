function [] = parameter_perturbation()

    %%% Command:
    %%% parameter_perturbation()  
   
    clc;
    close all;
    set(0, 'DefaultFigureVisible', 'on');

    metadata = load('../datamat.mat');

    %%% estimated parameters
    pars = readtable('maxLikValues.txt');
    pars = pars.Var2(2:end); % Best-fit parameters

    %%% fixed parameters 
    DRif                = 10;
    k_r                 = 0.049;
    k_pxrdeg            = 0.011;
    k_mRNAcyp3a4deg     = 0.044;
    k_mRNAcyp2c9fold    = 0.040;
    k_mRNAcyp2c9deg     = 0.036;
    k_mRNAcyp2b6deg     = 0.034;

    T = metadata.t_des;
    tmax = metadata.tmax_des;
    timeXaxis = [0, 12, 24, 48, 72, 120];
    YLims = {[-0.1 1.1]; [-0.35 4.35]; [-0.35 4.35]; [-0.35 4.35]; 
             [-0.1 1.1]; [-0.35 4.35]; [-0.35 4.35]; [-0.35 4.35];
             [-0.1 1.1]; [-0.35 4.35]; [-0.35 4.35]; [-0.35 4.35];
             [-0.1 1.1]; [-0.35 4.35]; [-0.35 4.35]; [-0.35 4.35];};    

    ylabels = {' activation';...
               ' mRNA fold change';...
               ' mRNA fold change';...
               ' mRNA fold change';...
               ' activation';...
               ' mRNA fold change';...
               ' mRNA fold change';...
               ' mRNA fold change';...
               ' activation';...
               ' mRNA fold change';...
               ' mRNA fold change';...
               ' mRNA fold change'};
    Text = {'PXR'; 'CYP3A4'; 'CYP2C9'; 'CYP2B6';
            'PXR'; 'CYP3A4'; 'CYP2C9'; 'CYP2B6';
            'PXR'; 'CYP3A4'; 'CYP2C9'; 'CYP2B6'};

    if ~exist('./figures', 'dir')
        mkdir('./figures');
    end

    param_names = {'$k_\mathrm{pxr,max}$',...
                   '$k_\mathrm{pxr,max}$',...
                   '$k_\mathrm{pxr,max}$',...
                   '$k_\mathrm{pxr,max}$',...
                   '$k_{mRNA_\mathrm{cyp3a4}^\mathrm{fold}}$',...
                   '$k_{mRNA_\mathrm{cyp3a4}^\mathrm{fold}}$',...
                   '$k_{mRNA_\mathrm{cyp3a4}^\mathrm{fold}}$',...
                   '$k_{mRNA_\mathrm{cyp3a4}^\mathrm{fold}}$',...
                   '$k_{mRNA_\mathrm{cyp2b6}^\mathrm{fold}}$',...
                   '$k_{mRNA_\mathrm{cyp2b6}^\mathrm{fold}}$',...
                   '$k_{mRNA_\mathrm{cyp2b6}^\mathrm{fold}}$',...
                   '$k_{mRNA_\mathrm{cyp2b6}^\mathrm{fold}}$'};

    %%% Sensitivity Analysis %%%
    sensitivity_range = [-0.5, 0, 0.5]; % -50%, best-fit, +50% 
    colors = [0 0 1; 0 0 0; 0 1 1]; % lower bound, best-fit, upper bound
    
    %%% Create one figure for all parameters
    tiledplot = tiledlayout(3,4,'TileSpacing','Compact');
    figlabs = {'(A)','(B)','(C)','(D)',...
               '(E)','(F)','(G)','(H)',...
               '(I)','(J)','(K)','(L)'};
    set(gcf, 'Position', [300, 100, 1300, 850]);
    num_rows = length(pars);
    num_cols = 4; % One for each variable (PXR, CYP3A4, etc.)
    
    mm = 0;
    nn = 0;
    for param_idx = 1:num_rows
        for var_idx = 1:num_cols
            mm = mm+1;
            ax(mm) = nexttile(mm);
            nn = nn + 1;
            set(ax(mm),...
                'box','on',...
                'XLim',[-0.05*max(timeXaxis) 1.05*max(timeXaxis)],...
                'XTick',timeXaxis,...
                'XTickLabel',timeXaxis,...
                'XTickLabelRotation',0,...
                'YLim',YLims{nn},...
                'FontSize',12);
            set(gca,'TickLength',[0.025, 0.01]);
            hold on;
        
            % Results storage
            outputs = cell(size(sensitivity_range));

            % Simulate for -50%, best-fit, and +50%  
            for i = 1:length(sensitivity_range)
                modified_pars = pars;
                modified_pars(param_idx) = pars(param_idx) * (1 + sensitivity_range(i));
                output = model(modified_pars, T, tmax);
                outputs{i} = output{var_idx};
            end

            % Plot curves
            lower_bound = outputs{1}; % -50%
            best_fit    = outputs{2}; % Best-fit
            upper_bound = outputs{3}; % +50%

            % Plot lower bound
            plot(T, lower_bound, '-','Color', colors(1, :), 'LineWidth', 1, 'DisplayName', '-50%');
            hold on;
            % Plot best-fit
            plot(T, best_fit, 'Color', colors(2, :), 'LineWidth', 1, 'DisplayName', 'Best-Fit');
            hold on;
            % Plot upper bound
            plot(T, upper_bound, '-', 'Color', colors(3, :), 'LineWidth', 1, 'DisplayName', '+50%');
            hold on;

            if var_idx == 1
                plot(ax(mm),T,zeros(length(T)),'Color',[0 0 0],'LineWidth',1,'LineStyle','--')
                hold on;
            else
                plot(ax(mm),T,ones(length(T)),'Color',[0 0 0],'LineWidth',1,'LineStyle','--')
                hold on;
            end

            xlabel(ax(mm), 'Time (hours)');
            title(ax(mm), strcat(Text{nn}, ylabels{nn}), 'FontSize', 12, 'FontWeight', 'Normal');
            hold on;

            text(0.1,0.9,figlabs{mm},...
                'Units','Normalized',...
                'HorizontalAlignment','center',...
                'FontSize',12,...
                'FontWeight','Normal');    
            hold on;

            text(0.3,0.88,param_names{mm},...
                'Units','Normalized',...
                'HorizontalAlignment','center',...
                'FontSize',12,...
                'FontWeight','Normal', ...
                'Interpreter','latex');    
            hold on;
            legend('-50%', 'best-fit', '+50%', 'Location', 'northeast', 'FontSize', 10);
            % end
        end
    end
    tiledplot.TileSpacing = 'compact';
    tiledplot.Padding = 'compact';

    % Save the combined figure
    savefig('figures/combined_parameter_perturbation.fig');
    exportgraphics(gcf, 'figures/combined_parameter_perturbation.png');
    exportgraphics(gcf, '../../LaTeX/figures/combined_parameter_perturbation.eps', 'ContentType', 'vector');


    %%% Helper Functions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [output] = model(pars, T, tmax)
        solution = ode23s(@ode, [0 tmax], [0 1 1 1], [], pars);
        for kk = 1:4
            output{kk} = deval(solution, T, kk);
        end
    end
    
    % ODE system
    function [dxdt] = ode(t, x, pars)
        dxdt = zeros(4, 1);

        dxdt(1) = pars(1) * DRif * (1 - x(1)) * exp(-k_r * t) - k_pxrdeg * x(1); % PXR
        dxdt(2) = pars(2) * x(1) + k_mRNAcyp3a4deg * (1 - x(2));               % mRNA CYP3A4
        dxdt(3) = k_mRNAcyp2c9fold * x(1) + k_mRNAcyp2c9deg * (1 - x(3));      % mRNA CYP2C9
        dxdt(4) = pars(3) * x(1) + k_mRNAcyp2b6deg * (1 - x(4));               % mRNA CYP2B6
    end

end