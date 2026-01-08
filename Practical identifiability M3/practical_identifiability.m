function [] = practical_identifiability()

    clc;
    clear all;
    close all;
    rng('shuffle');

    metadata        = load('../datamat.mat');
      
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
    k_mRNAcyp2c9deg     = 0.036;
    k_mRNAcyp2b6deg     = 0.034;

    MLpars_table    = readmatrix('../MC3_varying_transcription/maxLikValues.txt');
    MLpars          = transpose(MLpars_table(2:end,2));
    MLvalue         = loglike(MLpars([1]),MLpars([2,3,4]),([1]),([2,3,4]));
    MLvalues        = [MLvalue,MLpars];
    save('./MLvalues.mat','MLvalues');


    k_pxrmax_span        = [0.000005:0.000001:0.0005,0.0006:0.0001:0.015];
    k_mRNAcyp3a4_span    = 0.01:0.002:0.8;
    k_mRNAcyp2c9_span    = 0.01:0.002:0.8;
    k_mRNAcyp2b6_span    = 0.01:0.002:0.8;

    %%% k_pxrmax %%%
    funval = [];    
    idx_est = [2,3,4];
    idx_fix = [1];
    options = optimoptions(@fmincon,'Display','off','Algorithm','interior-point');
    for i=1:length(k_pxrmax_span)
        pars = [k_pxrmax_span(i), MLpars(2), MLpars(3), MLpars(4)];
        [x,fval,exitflag,output] = fmincon(@(params)loglike(params,pars(idx_fix),idx_est,idx_fix), pars(idx_est), [], [], [], [], 0*ones(length(idx_est),1), inf*ones(length(idx_est),1), [], options);
        funval = [funval; k_pxrmax_span(i), fval];
        disp(['k_pxr ',num2str(i),' done. Value = ',num2str(fval)])
    end    
    save('./k_pxrmax.mat','funval');
    disp('k_pxr done.')

    %%% k_mRNAcyp3a4 %%%
    funval = [];    
    idx_est = [1,3,4];
    idx_fix = [2];
    options = optimoptions(@fmincon,'Display','off','Algorithm','interior-point');
    for i=1:length(k_mRNAcyp3a4_span)
        pars = [MLpars(1), k_mRNAcyp3a4_span(i), MLpars(3), MLpars(4)];
        [x,fval,exitflag,output] = fmincon(@(params)loglike(params,pars(idx_fix),idx_est,idx_fix), pars(idx_est), [], [], [], [], 0*ones(length(idx_est),1), inf*ones(length(idx_est),1), [], options);
        funval = [funval; k_mRNAcyp3a4_span(i), fval];
        disp(['k_mRNAcyp3a4 ',num2str(i),' done. Value = ',num2str(fval)])
    end 
    save('./k_mRNAcyp3a4.mat','funval');
    disp('k_mRNAcyp3a4 done.')

    %%% k_mRNAcyp2c9 %%%
    funval = [];    
    idx_est = [1,2,4];
    idx_fix = [3];
    options = optimoptions(@fmincon,'Display','off','Algorithm','interior-point');
    for i=1:length(k_mRNAcyp2c9_span)
        pars = [MLpars(1), MLpars(2), k_mRNAcyp2c9_span(i), MLpars(4)];
        [x,fval,exitflag,output] = fmincon(@(params)loglike(params,pars(idx_fix),idx_est,idx_fix), pars(idx_est), [], [], [], [], 0*ones(length(idx_est),1), inf*ones(length(idx_est),1), [], options);
        funval = [funval;k_mRNAcyp2c9_span(i), fval];
        disp(['k_mRNAcyp2c9 ',num2str(i),' done. Value = ',num2str(fval)])
    end 
    save('./k_mRNAcyp2c9.mat','funval');
    disp('k_mRNAcyp2c9 done.')

    %%% k_mRNAcyp2b6 %%%
    funval = [];    
    idx_est = [1,2,3];
    idx_fix = [4];
    options = optimoptions(@fmincon,'Display','off','Algorithm','interior-point');
    for i=1:length(k_mRNAcyp2b6_span)
        pars = [MLpars(1), MLpars(2), MLpars(3), k_mRNAcyp2b6_span(i)];
        [x,fval,exitflag,output] = fmincon(@(params)loglike(params,pars(idx_fix),idx_est,idx_fix), pars(idx_est), [], [], [], [], 0*ones(length(idx_est),1), inf*ones(length(idx_est),1), [], options);
        funval = [funval; k_mRNAcyp2b6_span(i), fval];
        disp(['k_mRNAcyp2b6 ',num2str(i),' done. Value = ',num2str(fval)])
    end 
    save('./k_mRNAcyp2b6.mat','funval');
    disp('k_mRNAcyp2b6 done.')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% helper functions
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%% log-likelihood formula %%%
    function [value] = lognormpdf(data,model,stdev)
        value = -(1)*(-0.5*((data-model)./stdev).^2  - log(sqrt(2*pi).*stdev));
    end

    %%% log-likelihood %%%
    function [value] = loglike(pars_est,pars_fix,idx_est,idx_fix)  
        [solution] = model(pars_est,pars_fix,idx_est,idx_fix);
        value = 0;
        for ii=1:length(data)
            for jj=1:size(data{ii},1)
                idx=find(~isnan(data{ii}(jj,:)));
                value = value+sum(lognormpdf(data{ii}(jj,idx),solution{ii}(1,idx),stdev{ii}(1,idx)));
            end        
        end
    end

    %%% output of the model %%%
    function [output] = model(pars_est,pars_fix,idx_est,idx_fix)
        solution_des = ode23s(@odeDES,[0 tmax],...
                       [0 1 1 1],...
                       [],...
                       pars_est,...
                       pars_fix,...
                       idx_est,...
                       idx_fix);
        output{1} = deval(solution_des,time,2);
        output{2} = deval(solution_des,time,3);
        output{3} = deval(solution_des,time,4);
    end

    %%% ODE system - DESRIF %%%
    function [dxdt] = odeDES(t,x,pars_est,pars_fix,idx_est,idx_fix)        
        dxdt = zeros(4,1);

        [k_pxrmax,k_mRNAcyp3a4,k_mRNAcyp2c9,k_mRNAcyp2b6]=findpars(pars_est,pars_fix,idx_est,idx_fix);
        
        dxdt(1) = k_pxrmax*DRif*(1 - x(1))*exp(-k_r*t) - k_pxrdeg*x(1);                     % pxr
        dxdt(2) = k_mRNAcyp3a4*x(1) + k_mRNAcyp3a4deg*(1 - x(2));                               % mRNA CYP3A4 
        dxdt(3) = k_mRNAcyp2c9*x(1) + k_mRNAcyp2c9deg*(1 - x(3));                               % mRNA CYP2C9
        dxdt(4) = k_mRNAcyp2b6*x(1) + k_mRNAcyp2b6deg*(1 - x(4));                               % mRNA CYP2B6

    end

    function [k_pxrmax,k_mRNAcyp3a4,k_mRNAcyp2c9,k_mRNAcyp2b6]=findpars(pars_est,pars_fix,idx_est,idx_fix)
    
        if length(find(idx_fix==1))==0
            I=find(idx_est==1);
            k_pxrmax = pars_est(I);
        else
            I=find(idx_fix==1);
            k_pxrmax = pars_fix(I);
        end
        
        if length(find(idx_fix==2))==0
            I=find(idx_est==2);
            k_mRNAcyp3a4 = pars_est(I);
        else
            I=find(idx_fix==2);
            k_mRNAcyp3a4 = pars_fix(I);
        end
        
        if length(find(idx_fix==3))==0
            I=find(idx_est==3);
            k_mRNAcyp2c9 = pars_est(I);
        else
            I=find(idx_fix==3);
            k_mRNAcyp2c9 = pars_fix(I);
        end
        
        if length(find(idx_fix==4))==0
            I=find(idx_est==4);
            k_mRNAcyp2b6 = pars_est(I);
        else
            I=find(idx_fix==4);
            k_mRNAcyp2b6 = pars_fix(I);
        end        
       
    end



end