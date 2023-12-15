function [ACT layerNames] = load_AE_activ(cfg, sessi, subj_ch_fr, paths, net2load);

    %lays2load = cfg.lays2load;
    lays2load = [1 5 10 14 19 23 27 32 37 42 44 46 48 50 52 54 56 58 60 66];
    
    currentFolder = pwd; 
    
    gPath = [paths.activations.AE net2load]; 
    
    if sessi < subj_ch_fr
        cd ([gPath '\freiburg'])
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
        cd ([gPath '\china'])
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