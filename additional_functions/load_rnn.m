function[ACT] = load_rnn(cfg, sessi, subj_ch_fr, paths);%load network if not loaded yet


currentFolder = pwd; 
cd ([paths.activations '\imageNet\spearman\'])

if sessi < subj_ch_fr
    cd FR
    sublist = dir('*FR.mat');
    sublist = {sublist.name};  sublist = sort(sublist');

    for layi=1:length(cfg.lays2load)
        load(sublist{cfg.lays2load(layi)});
        fs2= squeeze(a(1, :, :)); 
        ACT(layi,:,:) = fs2;   

    end            
else
    cd CH
    sublist = dir('*CH.mat');
    sublist = {sublist.name};  sublist = sort(sublist');

    tic

    for layi=1:length(cfg.lays2load)
        load(sublist{cfg.lays2load(layi)});
        fs2= squeeze(a(1, :, :)); 
        ACT(layi,:,:) = fs2;     

    end

    cd ..
end
cd ..
cd (currentFolder)



end

