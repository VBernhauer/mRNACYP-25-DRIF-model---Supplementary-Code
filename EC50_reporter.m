function [] = EC50_reporter()

    close all;
    clear all;

    if ~exist('./figures', 'dir')
        mkdir('./figures')
    end

    metadata = load('./datamat.mat');

    mean_reporter_12h = metadata.mean_reporter_12h;
    std_reporter_12h = metadata.std_reporter_12h;

    mean_reporter_24h = metadata.mean_reporter_24h;
    std_reporter_24h = metadata.std_reporter_24h;

    mean_reporter_48h = metadata.mean_reporter_48h;
    std_reporter_48h = metadata.std_reporter_48h;

    mean_reporter_72h = metadata.mean_reporter_72h;
    std_reporter_72h = metadata.std_reporter_72h;

    data = [mean_reporter_12h; mean_reporter_24h; mean_reporter_48h; mean_reporter_72h];
    stdv = [std_reporter_12h; std_reporter_24h; std_reporter_48h; std_reporter_72h];

    logC_reporter = log10([0.1, 0.5, 1, 2.5, 5, 10, 20, 30]);

    logC = -3:0.01:2;
    logCfit = -3:0.01:log10(30);

    CYPmin  = 1;
    CYPmax  = 3;
    N = 1.5;
    logEC50 = 0.01;
    hill_pars = [CYPmax, N, logEC50];
    opts = optimset('Display','off');

    EC50_pars_out   = [];
    interp          = [];
    for ii = 1:size(data,1)
        [EC50_pars, resnorm]    = lsqnonlin(@(pars)error(pars,data(ii,:),stdv(ii,:),logC_reporter),hill_pars,[],[],opts);
        EC50_pars_out           = [EC50_pars_out; EC50_pars];
        interp                  = [interp; csapi(logC_reporter,hill(EC50_pars,logC_reporter))];
    end

    fontsize    = 12;
    ec_fontsize = 10;
    varmin      = 0;
    varmax      = 10;
    colors      = {[0.7176 0.2745 1], [0.9804 0.8667 0.2353], [0.9294 0.6941 0.1255], [0.5 0.5 0.5]};
    name        = {'12 h, EC_{50} =', '24 h, EC_{50} =', '48 h, EC_{50} =', '72 h, EC_{50} ='};
    T           = [12,24,48,72];
    yax         = [0.6,0.5,0.4,0.3];

    figure();
    ax = gca;
    set(ax,...
        'box','on',...
        'XLim',[logC(1)-0.15 logC(end)+0.15],...
        'XTick',-3:0.5:2,...
        'YLim',[varmin-0.45 varmax+0.45],...
        'YTick',[0,2,4,6,8,10],...
        'FontSize',fontsize);

    xlabel('Log_{10} concentration of RIF (µM)','FontSize',fontsize);
    ylabel('CYP3A4 gene reporter activity','FontSize',fontsize);
    set(gca,'TickLength',[0.015, 0.015])
    hold on;  

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
    line_1=plot([logC(1) logC(end)], [1 1],...
          'Color',[0 0 0], ...
          'LineWidth',1.0,...
          'LineStyle','--');


    %%% data + errors %%%
    for ii = 1:size(data,1)
        er(ii) = errorbar(logC_reporter, data(ii,:), stdv(ii,:),...
              'v',...
              'Color',colors{ii},...
              'MarkerFaceColor',colors{ii}, ...
              'MarkerSize',8,...
              'MarkerEdgeColor',[0.25 0.25 0.25],...
              'LineWidth',1.0,...
              'DisplayName',[name{ii}, ' ', num2str(10.^EC50_pars_out(ii,3),'%4.2f'),' µM']);
        hold on;
        plot(logCfit, hill(EC50_pars_out(ii,:),logCfit),...
          'Color',colors{ii},...
          'LineWidth',1.01);

        % text(0.35,yax(ii),strcat('EC_{50,', num2str(T(ii)),' h}^{RIF} = ',{' '},num2str(10.^EC50_pars_out(ii,3),'%4.2f'),' \muM'),...
        %     'Units','Normalized',...
        %     'HorizontalAlignment','right',...
        %     'FontSize',ec_fontsize,...
        %     'FontWeight','Normal',...
        %     'Color',colors{ii}); 
    
    end
    leg1 = legend(er);
    pos = get(leg1,'Position');
    pos(1) = 0.3*pos(1);
    pos(2) = 0.95*pos(2);
    set(leg1,'Position',pos);
    set(leg1,'FontSize',ec_fontsize);

    savefig(strcat('./figures/EC50_reporters.fig'));
    exportgraphics(gcf,strcat('./figures/EC50_reporters.png'));
    exportgraphics(gcf,strcat('./figures/EC50_reporters.eps'),'ContentType','vector');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function [CYP] = hill(pars,logC)
        CYP = CYPmin + (pars(1) - CYPmin) ./ (1 + 10 .^ (pars(2) * (pars(3) - logC)));
    end

    function [output] = error(pars,data,stdv,logC)
        CYP     = hill(pars,logC);
        output  = (CYP - data).^2 ./ stdv.^2;
    end


end