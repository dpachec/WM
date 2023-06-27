function[ACT] = load_CORrt_TL(cfg, sessi, subj_ch_fr, paths);%load network if not loaded yet


currentFolder = pwd; 
cd ([paths.activations.corNet_rtTL])



if sessi < subj_ch_fr
    sublist = dir('*FR.mat');
    sublist = {sublist.name};  
    load(sublist{1});
    if length(cfg.lays2load) > 1
        ACT = squeeze(a(cfg.lays2load, :, :)); 
    else
        ACT = a(cfg.lays2load, :, :); 
    end
    

else
    
    sublist = dir('*CH.mat');
    sublist = {sublist.name};  
    load(sublist{1});
    if length(cfg.lays2load) > 1
        ACT = squeeze(a(cfg.lays2load, :, :)); 
    else
        ACT = a(cfg.lays2load, :, :); 
    end
    
end

cd (currentFolder)


