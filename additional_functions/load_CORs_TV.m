
function [ACT] = load_CORs_activ(cfg, sessi, subj_ch_fr, paths);

    lays2load = cfg.lays2load;
    
    currentFolder = pwd; 
    
    gPath = paths.activations.corNet_sTV; 
    
    if sessi < subj_ch_fr
        cd ([gPath 'freiburg'])
        fold2load = dir('*cornet*'); fold2load= {fold2load.name}';
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













% 
% function [ACT] = load_CORs_activ(cfg, sessi, subj_ch_fr, paths);
% 
%     lays2load = cfg.lays2load;
%     
%     gPath = paths.activations.corNet_s; 
%     
%     if sessi < subj_ch_fr
%         cd ([gPath 'freiburg'])
%         fold2load = dir('*activations*'); fold2load= {fold2load.name}';
%         fold2load1 = cellfun(@(x) strsplit(x, '_'), fold2load, 'un', 0); 
%         fold2load2 = cellfun(@(x) x{2}, fold2load1, 'un', 0); 
%         x = double(string(fold2load2));
%         [id1 id2] = sort(x);
%         fold2load3  = fold2load(id2);
%         
%         for layi=1:length(fold2load)
%             load(fold2load3{layi})
%             b = reshape(a, 60, []);
%             rdm = corr(b', 'type','s');
%             ACT(layi,:,:) = rdm; 
%         end            
%     else
%         cd ([gPath 'china'])
%         fold2load = dir('*activations*'); fold2load= {fold2load.name}';
%         fold2load1 = cellfun(@(x) strsplit(x, '_'), fold2load, 'un', 0); 
%         fold2load2 = cellfun(@(x) x{2}, fold2load1, 'un', 0); 
%         x = double(string(fold2load2));
%         [id1 id2] = sort(x);
%         fold2load3  = fold2load(id2);
%         for layi=1:length(fold2load)
%             load(fold2load3{layi})
%             b = reshape(a, 60, []);
%             rdm = corr(b', 'type','s');
%             ACT(layi,:,:) = rdm; 
%         end   
%         
%     end
% 
%     cd (paths.github)
%       
% end