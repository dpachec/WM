function [ACT layerNames] = load_AE_ECO_activ(cfg, sessi, paths, net2load);

    %lays2load = cfg.lays2load;
    lays2load = [0 1 5 10 14 19 23 27 32 37 42 46 48 50 52 54 56 58 60 62 67 68] +1;
    
    currentFolder = pwd; 
    
    gPath = [paths.activations.AE 'eco\']; 
    normalizedMatrices = strcmp(net2load(end), 'N')
    if normalizedMatrices
        net2load = [net2load(1) '_' net2load(2:end-2)];
    else
        net2load = [net2load(1) '_' net2load(2:end-1)];
    end

        
    

    if normalizedMatrices
        cd ([gPath net2load '\mean_normalized_rsm\'])
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
        cd ([gPath net2load '\default_rsm\'])
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
    




    cd (currentFolder)



end