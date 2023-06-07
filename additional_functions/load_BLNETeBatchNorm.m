function[ACT] = load_BLNETeBatchNorm(cfg, sessi, subj_ch_fr, paths);%load network if not loaded yet


currentFolder = pwd; 
cd ([paths.activations.BLNETeBatchNorm])


if sessi < subj_ch_fr
    sublist = dir('*FR.mat');
    sublist = {sublist.name};  sublist = sort(sublist');
    load(sublist{1});
    fs2= squeeze(a(cfg.lays2load, :, :)); 
    ACT = fs2;   
else
    sublist = dir('*CH.mat');
    sublist = {sublist.name};  sublist = sort(sublist');
    load(sublist{1});
    fs2= squeeze(a(cfg.lays2load, :, :)); 
    ACT = fs2;   

end

cd (currentFolder)



end

