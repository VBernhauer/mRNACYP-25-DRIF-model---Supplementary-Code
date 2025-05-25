function [results] = mcmc(nstart,nend,pars)

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
    k_mRNAcyp3a4fold    = 0.083;
    k_mRNAcyp3a4deg     = 0.044;
    k_mRNAcyp2c9fold    = 0.040;
    k_mRNAcyp2c9deg     = 0.036;
    k_mRNAcyp2b6fold    = 0.139;
    k_mRNAcyp2b6deg     = 0.034;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% SLICE SAMPLER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for iround=nstart:nend
        nSamples = 50000;
        burnIn = 25000;
        thin = 5;
        chains = slicesample(pars(iround,:),nSamples,'logpdf',@logprob);
        chains = chains(burnIn+1:thin:end,:);
        results{iround} = chains;
        
        %%% save chains %%%
        save(strcat('./chains/chains_',num2str(iround),'.mat'),'chains');
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% helper functions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% log-likelihood formula %%%
    function [value] = lognormpdf(data,model,stdev)
        value = -0.5*((data-model)./stdev).^2  - log(sqrt(2*pi).*stdev);
    end

    %%% log-probability %%%
    function [value] = logprob(pars)  
        lp = logprior(pars);
        if ~isfinite(lp)
            value = -inf;
        else
            value = lp+loglike(pars);
        end   
    end

    %%% priors %%%
    function [flag] = logprior(pars)
        if all(pars>0)
            flag = 0;
        else
            flag = -inf;
        end    
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
        
        dxdt(1) = pars*DRif*(1 - x(1))*exp(-k_r*t) - k_pxrdeg*x(1);                                 % pxr
        dxdt(2) = k_mRNAcyp3a4fold*x(1) + k_mRNAcyp3a4deg*(1 - x(2));                               % mRNA CYP3A4 
        dxdt(3) = k_mRNAcyp2c9fold*x(1) + k_mRNAcyp2c9deg*(1 - x(3));                               % mRNA CYP2C9
        dxdt(4) = k_mRNAcyp2b6fold*x(1) + k_mRNAcyp2b6deg*(1 - x(4));                               % mRNA CYP2B6

    end
    
end