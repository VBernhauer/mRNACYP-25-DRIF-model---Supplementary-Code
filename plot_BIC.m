function [] = plot_BIC()

    %%% command:
    %%% plot_BIC()
   
    clc;
    close all
    set(0,'DefaultFigureVisible','on');

    if ~exist('./figures', 'dir')
        mkdir('./figures')
    end
    
    num_of_datapoints = 45;

	num_of_parameters_MC1 = 1;    
    num_of_parameters_MC2 = 3; 
    num_of_parameters_MC3 = 4;
    
    loglik_MC1 	= load('./fixed_transcription/mlvalues.mat');
    loglik_MC2 	= load('./fixed_cyp2c9_transcription/mlvalues.mat');
    loglik_MC3 	= load('./varying_transcription/mlvalues.mat');


	loglik_MC1 	= loglik_MC1.mlvalue_array(:,:);    
    loglik_MC2	= loglik_MC2.mlvalue_array(:,:);
    loglik_MC3 	= loglik_MC3.mlvalue_array(:,:);

    BIC_MC1 	= num_of_parameters_MC1 * log(num_of_datapoints) - 2 * loglik_MC1;    
    BIC_MC2 	= num_of_parameters_MC2 * log(num_of_datapoints) - 2 * loglik_MC2;
    BIC_MC3 	= num_of_parameters_MC3 * log(num_of_datapoints) - 2 * loglik_MC3; 

    [p,h,stats] = ranksum(BIC_MC1, BIC_MC2)
    [p,h,stats] = ranksum(BIC_MC1, BIC_MC3)
    [p,h,stats] = ranksum(BIC_MC2, BIC_MC3)

    fontsize = 14;

    figure();
    ax = gca;
    set(gcf,'Position', [300, 100, 700, 300]);
    set(ax,...
        'box','on',...
        'XLim',[60 130],...
        'FontSize',fontsize);
    set(gca,'TickLength',[0.025, 0.01])
    hold on;
    histogram(BIC_MC1,...
          'FaceColor', 'none',...
          'EdgeColor', [0 1 1],...
          'FaceAlpha', 1,...
          'EdgeAlpha', 1,...
          'Normalization', 'probability',...
          'DisplayStyle', 'stairs',...
          'LineWidth',1.5);
    hold on;
    histogram(BIC_MC2,...
          'FaceColor', 'none',...
          'EdgeColor', [0 0 0],...
          'FaceAlpha', 1,...
          'EdgeAlpha', 1,...
          'Normalization', 'probability',...
          'DisplayStyle', 'stairs',...
          'LineWidth',1.5);
    hold on;
    histogram(BIC_MC3,...
          'FaceColor', 'none',...
          'EdgeColor', [0.75 0.75 0.75],...
          'FaceAlpha', 1,...
          'EdgeAlpha', 1,...
          'Normalization', 'probability',...
          'DisplayStyle', 'stairs',...
          'LineWidth',1.5);
    hold on;
    l = legend('MC1','MC2','MC3','Interpreter','latex');
    set(l,'Location','northwest');
    set(l,'FontSize',fontsize);   
    xlabel('Bayesian Information Criterion','FontSize',fontsize);
    ylabel('Normalized frequency','FontSize',fontsize);

	savefig(strcat('./figures/BIC.fig'));
    exportgraphics(gcf,'./figures/BIC.png');
    exportgraphics(gcf,'./figures/BIC.eps','ContentType','vector');

    end