function [] = posterior_values_ci()

    %%% command:
    %%% posterior_values_ci()
   
    clc;
    close all
    set(0,'DefaultFigureVisible','on');
    
    metadata = load('../datamat.mat');

    data    = {metadata.dataCYP3A4_des;...
               metadata.dataCYP2C9_des;...
               metadata.dataCYP2B6_des};
    stdev   = {metadata.std_dataCYP3A4_des;...
               metadata.std_dataCYP2C9_des;...
               metadata.std_dataCYP2B6_des};
    time    = metadata.time_des;
    tmax    = metadata.tmax_des;
    
    %%% fixed parameters %%%
    DRif                = 10;
    k_r                 = 0.049;
    k_pxrdeg            = 0.011;
    k_mRNAcyp3a4deg     = 0.044;
    k_mRNAcyp2c9fold    = 0.040;
    k_mRNAcyp2c9deg     = 0.036;
    k_mRNAcyp2b6deg     = 0.034;

 
    if ~exist('./figures', 'dir')
        mkdir('./figures')
    end
    
    chains = [];
    for kk = 1:5
        kkchains = load(strcat('./chains/chains_',num2str(kk),'.mat'));
        chains = [chains; kkchains.chains(:,:);];
    end
    
    mlvalue = -99999.99;
    mlpars = zeros(1,size(chains,2));
    mlvalue_array = [];
    for ll = 1:size(chains,1)
        if mod(ll,100)==0
            disp(ll);
        end
        llvalue = loglike(chains(ll,:));
        mlvalue_array = [mlvalue_array; llvalue];
        if llvalue > mlvalue
            mlvalue = llvalue;
            mlpars = chains(ll,:);
        end
    end
    %%% save chains %%%
    save(strcat('./mlvalues.mat'),'mlvalue_array');
            
    parsmean = mean(chains);
    parsmedian = median(chains);
    parsquantile = [];
    for kk = 1:size(chains,2)
        parsquantile = [parsquantile; quantile(chains(:,kk),[0.05 0.95])];
    end
    
    
    a{1} = 'k_pxrmet';
    a{2} = 'k_mRNA_cyp3a4fold';
    a{3} = 'k_mRNA_cyp2b6fold';
    labs = {a{1},...
            a{2},...
            a{3}};
        
    names = {'parameter','mean','median','95CI lower', '95CI upper'};
    
    fid = fopen(strcat('posteriorValues.txt'),'w');

    fprintf(fid, '%2s %2s %2s %2s %2s\n', names{:});
    for kk = 1:length(labs)
        fprintf(fid,'%0s %.8f %.8f %.8f %.8f\n',labs{kk},parsmean(kk),parsmedian(kk),parsquantile(kk,1),parsquantile(kk,2));
    end
    fclose(fid);
    
    fidml = fopen(strcat('maxLikValues.txt'),'w');

    fprintf(fidml, '%2s %.8f \n', 'MLvalue', mlvalue);
    for kk = 1:length(labs)
        fprintf(fidml,'%0s %.8f\n',labs{kk},mlpars(kk));
    end
    fclose(fidml);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% log-likelihood formula %%%
    function [value] = lognormpdf(data,model,stdev)
        value = -0.5*((data-model)./stdev).^2  - log(sqrt(2*pi).*stdev);
    end

    %%% log-likelihood %%%
    function [value] = loglike(pars)        
        [solution] = model(pars);
        value = 0;
        for ii=1:length(data)
            for jj=1:size(data{ii},1)
                idx=find(~isnan(data{ii}(jj,:)));
                value = value+sum(lognormpdf(data{ii}(jj,idx),solution{ii}(1,idx),stdev{ii}(1,idx)));
            end        
        end
    end

    %%% output of the model %%%
    function [output] = model(pars)
        solution_des = ode23s(@odeDES,[0 tmax],...
                       [0 1 1 1],...
                       [],...
                       pars);
        output{1} = deval(solution_des,time,2);
        output{2} = deval(solution_des,time,3);
        output{3} = deval(solution_des,time,4);
    end

    %%% ODE system - DESRIF %%%
    function [dxdt] = odeDES(t,x,pars)        
        dxdt = zeros(4,1);
        
        dxdt(1) = pars(1)*DRif*(1 - x(1))*exp(-k_r*t) - k_pxrdeg*x(1);                     % pxr
        dxdt(2) = pars(2)*x(1) + k_mRNAcyp3a4deg*(1 - x(2));                               % mRNA CYP3A4 
        dxdt(3) = k_mRNAcyp2c9fold*x(1) + k_mRNAcyp2c9deg*(1 - x(3));                      % mRNA CYP2C9
        dxdt(4) = pars(3)*x(1) + k_mRNAcyp2b6deg*(1 - x(4));                               % mRNA CYP2B6

    end

end