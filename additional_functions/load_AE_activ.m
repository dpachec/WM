function [ACT layerNames] = load_AE_activ(cfg, sessi, subj_ch_fr, paths, net2load);

    %lays2load = cfg.lays2load;
    lays2load = [0 1 5 10 14 19 23 27 32 37 42 46 48 50 52 54 56 58 60 62 67 68] +1;
    
    currentFolder = pwd; 
    
    gPath = [paths.activations.AE]; 
    normalizedMatrices = strcmp(net2load(end), 'N');
    net2load = [net2load(3) '_' net2load(4:end-1)];
    
    
    if sessi < subj_ch_fr
        if normalizedMatrices
            cd ([gPath '\freiburg\' net2load '\mean_normalized_rsm\' ])
            f2load = dir('*.csv*'); f2load= {f2load.name}';
            count = 1; 
            for layi= lays2load%1:length(fold2load)
                rdm = csvread(f2load{layi});
                ACT(count,:,:) = rdm; 
                lN = erase(f2load{layi}, '.csv');
                layerNames{count,:} = lN(5:end);
                count = count+1; 
            end    
        else
            cd ([gPath '\freiburg\' net2load '\default_rsm\' ])
            f2load = dir('*.csv*'); f2load= {f2load.name}';
            count = 1; 
            for layi= lays2load%1:length(fold2load)
                rdm = csvread(f2load{layi});
                ACT(count,:,:) = rdm; 
                lN = erase(f2load{layi}, '.csv');
                layerNames{count,:} = lN(5:end);
                count = count+1; 
            end    
        end
    else
        if normalizedMatrices
            cd ([gPath '\china\' net2load '\mean_normalized_rsm'])
            f2load = dir('*.csv*'); f2load= {f2load.name}';
            count = 1; 
            for layi= lays2load%1:length(fold2load)
                rdm = csvread(f2load{layi});
                ACT(count,:,:) = rdm; 
                lN = erase(f2load{layi}, '.csv');
                layerNames{count,:} = lN(5:end);
                count = count+1; 
            end    
        else
            cd ([gPath '\china\' net2load '\default_rsm\'])
            f2load = dir('*.csv*'); f2load= {f2load.name}';
            count = 1; 
            for layi= lays2load%1:length(fold2load)
                rdm = csvread(f2load{layi});
                ACT(count,:,:) = rdm; 
                lN = erase(f2load{layi}, '.csv');
                layerNames{count,:} = lN(5:end);
                count = count+1; 
            end    

        end
        
    end



    cd (currentFolder)



end