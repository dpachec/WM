
function [ACT] = load_CORrt_activ(cfg, sessi, subj_ch_fr, paths);

    lays2load = cfg.lays2load;
    
    currentFolder = pwd; 
    
    gPath = paths.activations.corNet_rt; 
    
    if sessi < subj_ch_fr
        cd ([gPath 'freiburg'])
        fold2load = dir('*cornet*'); fold2load= {fold2load.name}';
        fold2load  = fold2load([4 3 6 5 8 7 2 1]);
        for layi=1:length(fold2load)
            cd (fold2load{layi})
            b = readNPY('features.npy');
            rdm = corr(b', 'type','s');
            ACT(layi,:,:) = rdm; 
            cd ..
        end            
    else
       cd ([gPath 'china'])
        fold2load = dir('*cornet*'); fold2load= {fold2load.name}';
        fold2load  = fold2load([4 3 6 5 8 7 2 1]);
        for layi=1:length(fold2load)
            cd (fold2load{layi})
            b = readNPY('features.npy');
            rdm = corr(b', 'type','s');
            ACT(layi,:,:) = rdm; 
            cd ..
        end   
        
    end

    cd (currentFolder)
      
end