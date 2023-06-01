function[ACT] = load_BLNETi(cfg, sessi, subj_ch_fr, paths);%load network if not loaded yet


currentFolder = pwd; 
cd ([paths.activations.BLNETi])



if sessi < subj_ch_fr
    sublist = dir('*FR.mat');
    sublist = {sublist.name};  sublist = sort(sublist');
    load(sublist{1});
    ACT = a(cfg.lays2load, :, :);
    

else
    
    sublist = dir('*CH.mat');
    sublist = {sublist.name};  sublist = sort(sublist');
    load(sublist{1});
    ACT = a(cfg.lays2load, :, :);   
    
end

cd (currentFolder)





% % % % 
% % % % if sessi < subj_ch_fr
% % % %     cd FR
% % % %     sublist = dir('*FR.mat');
% % % %     sublist = {sublist.name};  sublist = sort(sublist');
% % % % 
% % % %     for layi=1:length(cfg.lays2load)
% % % %         load(sublist{cfg.lays2load(layi)});
% % % %         fs2= squeeze(a(1, :, :)); 
% % % %         ACT(layi,:,:) = fs2;   
% % % % 
% % % %     end            
% % % % else
% % % %     cd CH
% % % %     sublist = dir('*CH.mat');
% % % %     sublist = {sublist.name};  sublist = sort(sublist');
% % % % 
% % % %     tic
% % % % 
% % % %     for layi=1:length(cfg.lays2load)
% % % %         load(sublist{cfg.lays2load(layi)});
% % % %         fs2= squeeze(a(1, :, :)); 
% % % %         ACT(layi,:,:) = fs2;     
% % % % 
% % % %     end
% % % % 
% % % %     cd ..
% % % % end
% % % % cd ..
% % % % cd (currentFolder)
% % % % 


end

