function [results] = mcmc(nstart,nend,pars)

    metadata = load('datamat.mat');
    
    data    = {metadata.data_deg;metadata.data_l;metadata.data_h};
    stdev   = {metadata.stdev_deg;metadata.stdev_l;metadata.stdev_h};
    time    = {metadata.time_deg;metadata.time_l;metadata.time_h};
    tmax    = {metadata.tmax_l;metadata.tmax_h};
    Rif     = {1;10};

    data = vertcat(data{:});
    stdev = vertcat(stdev{:});

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% SLICE SAMPLER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for iround=nstart:nend
        nSamples = 100000;
        burnIn = 50000;
        thin = 10;
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
        lp = prior(pars);
        if ~isfinite(lp)
            value = -inf;
        else
            value = lp+loglike(pars);
        end   
    end

    %%% priors %%%
    function [flag] = prior(pars)
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
        %%% degradation of protein mRNA %%%
        output{1} = exp(-pars(5) .* time{1});
        output{2} = exp(-pars(7) .* time{1});
        output{3} = exp(-pars(9) .* time{1});
        
        solution_l = ode23s(@ode,[0 tmax{1}],...
                       [0 1 1 1 1],...
                       [],...
                       pars,...
                       Rif{1});
        solution_h = ode23s(@ode,[0 tmax{2}],...
                       [0 1 1 1 1],...
                       [],...
                       pars,...
                       Rif{2});           
        %%% mRNA kinetics for measured proteins %%%
        % low
        kk = 4;
        for ii = 2:5
            output{kk} = deval(solution_l,time{2},ii);
            kk = kk + 1;
        end
        % high
        for ii = 2:5
            output{kk} = deval(solution_h,time{3},ii);
            kk = kk + 1;
        end   
    end

    %%% ODE system %%%
    function [dxdt] = ode(t,x,pars,Xint)        
        dxdt = zeros(5,1);

        dxdt(1) = pars(1)*(1-x(1))*Xint*exp(-pars(2)*t) - pars(3)*x(1);    % activated PXR
        dxdt(2) = pars(4)*x(1) + pars(5)*(1-x(2));                         % CYP3A4
        dxdt(3) = pars(6)*x(1) + pars(7)*(1-x(3));                         % CYP2C9
        dxdt(4) = pars(8)*x(1) + pars(9)*(1-x(4));                         % CYP2B6
        dxdt(5) = pars(10)*x(1) + pars(11)*(1-x(5));                       % MDR1
    end

end