function [ACT] = load_AE_activ(cfg, sessi, subj_ch_fr, paths, net2load);

    lays2load = cfg.lays2load;
    
    currentFolder = pwd; 
    
    gPath = [paths.activations.AE net2load]; 
    
    if sessi < subj_ch_fr
        cd ([gPath '\freiburg'])
        f2load = dir('*.csv*'); f2load= {f2load.name}';
        count = 1; 
        for layi= cfg.lays2load%1:length(fold2load)
            rdm = csvread(f2load{layi});
            ACT(count,:,:) = rdm; 
            count = count+1; 
        end    
    else
        cd ([gPath '\china'])
        f2load = dir('*.csv*'); f2load= {f2load.name}';
        count = 1; 
        for layi= cfg.lays2load%1:length(fold2load)
            rdm = csvread(f2load{layi});
            ACT(count,:,:) = rdm; 
            count = count+1; 
        end    
        
    end



    cd (currentFolder)



end