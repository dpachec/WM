function[ACT] = load_CORrt_TL(cfg, sessi, subj_ch_fr, paths);%load network if not loaded yet


currentFolder = pwd; 
cd ([paths.activations.corNet_rtTL])



if sessi < subj_ch_fr
    sublist = dir('*FR.mat');
    sublist = {sublist.name};  
    load(sublist{1});
    if length(cfg.lays2load) > 1
         a{6} = nan(60); a{11} = nan(60); a{16} = nan(60);
         act_prev = a(cfg.lays2load);
         ACT = cat(3, act_prev{:});
         ACT = permute(ACT, [3, 1, 2]);
        %a(6, :, :) = nan(60); a(11, :, :) = nan(60); a(16, :, :) = nan(60);
        %act_prev = a(cfg.lays2load, :, :);
        %ACT = act_prev; 

    else
        ACT = a{cfg.lays2load}; 
    end
    

else
    
    sublist = dir('*CH.mat');
    sublist = {sublist.name};  
    load(sublist{1});
    if length(cfg.lays2load) > 1
        a{6} = nan(60); a{11} = nan(60); a{16} = nan(60);
        act_prev = a(cfg.lays2load);
        ACT = cat(3, act_prev{:});
        ACT = permute(ACT, [3, 1, 2]);
        
        %a(6, :, :) = nan(60); a(11, :, :) = nan(60); a(16, :, :) = nan(60);
        %act_prev = a(cfg.lays2load, :, :);
        %ACT = act_prev; 

    else
        ACT = a{cfg.lays2load}; 
    end
    
end

cd (currentFolder)


