function [] = Mann_Whitney_test()

    %%% command:
    %%% Mann_Whitney_test()
   
    clc;
    close all
    set(0,'DefaultFigureVisible','on');

    rif_chains = [];
    for jj = 1:5
        jjchains = load(strcat('./rifampicin/chains/chains_',num2str(jj),'.mat'));
        jjchains = jjchains.chains(:,:);
        for kk = size(jjchains,1)/2:size(jjchains,1)
            rif_chains = [rif_chains;log10(jjchains(kk,:))];
        end
    end

    drif_chains = [];
    for jj = 1:5
        jjchains = load(strcat('./MC2_fixed_cyp2c9_transcription/chains/chains_',num2str(jj),'.mat'));
        jjchains = jjchains.chains(:,:);
        for kk = size(jjchains,1)/2:size(jjchains,1)
            drif_chains = [drif_chains;log10(jjchains(kk,:))];
        end
    end

    %%% k_pxr,max
    [p,h,stats] = ranksum(rif_chains(:,1), drif_chains(:,1))

    %%% k_cyp3a4,fold
    [p,h,stats] = ranksum(rif_chains(:,4), drif_chains(:,2))

    %%% k_cyp2b6,fold
    [p,h,stats] = ranksum(rif_chains(:,8), drif_chains(:,3))

end