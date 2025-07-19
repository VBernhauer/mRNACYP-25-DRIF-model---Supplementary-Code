function [results] = EC50_minmax(nstart,nend,T)

    clc;
    close all
    set(0,'DefaultFigureVisible','on');
    
    logC = -3:0.01:2;
    name = 'cyp3a4';


    for idx_T=nstart:nend

        data_rif_min  = [];
        data_rif_max  = [];
        data_drif_min = [];
        data_drif_max = [];
    
        for idx_C = 1:length(logC)
    
            % disp((num2str(logC(idx_C))));
    
            data_rif = load(strcat('./EC50/',name,'_time_rif_',num2str(logC(idx_C)),'.mat'));
            var_name = strcat('outputvec_rif_',name);
            data_rif = data_rif.(var_name);
    
            data_rif_min  = [data_rif_min, min(data_rif(:,idx_T))];
            data_rif_max  = [data_rif_max, max(data_rif(:,idx_T))];
    
            data_drif = load(strcat('./EC50/',name,'_time_drif_',num2str(logC(idx_C)),'.mat'));
            var_name = strcat('outputvec_drif_',name);
            data_drif = data_drif.(var_name);
    
            data_drif_min  = [data_drif_min, min(data_drif(:,idx_T))];
            data_drif_max  = [data_drif_max, max(data_drif(:,idx_T))];
    
        end
    
        save(strcat('./EC50_minmax/EC50_rif_min_T_',num2str(T(idx_T)),'.mat'),'data_rif_min');
        save(strcat('./EC50_minmax/EC50_rif_max_T_',num2str(T(idx_T)),'.mat'),'data_rif_max');
        save(strcat('./EC50_minmax/EC50_drif_min_T_',num2str(T(idx_T)),'.mat'),'data_drif_min');
        save(strcat('./EC50_minmax/EC50_drif_max_T_',num2str(T(idx_T)),'.mat'),'data_drif_max');

        results{idx_T} = [];
       
    end

end