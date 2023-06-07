
function [ACT] = load_alexECO_activ(cfg, sessi, subj_ch_fr, paths);

    lays2load = cfg.lays2load;
    
    currentFolder = pwd; 
    
    gPath = paths.activations.AlexEco; 
    
    if sessi < subj_ch_fr
        cd ([gPath 'FR'])
        fold2load = dir('*Alexnet*'); fold2load= {fold2load.name}';
        fold2load  = fold2load([4 6 7 8 5 1 2 3]);
        for layi=1:length(fold2load)
            cd (fold2load{layi})
            b = readNPY('features.npy');
            rdm = corr(b', 'type','s');
            ACT(layi,:,:) = rdm; 
            cd ..
        end            
    else
       cd ([gPath 'CH'])
        fold2load = dir('*Alexnet*'); fold2load= {fold2load.name}';
        fold2load  = fold2load([4 6 7 8 5 1 2 3]);
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