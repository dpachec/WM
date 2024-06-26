function [EEG] = select_electrodes(EEG, region)

atlas = my_ft_read_atlas('D:\Documents\GitHub\WM\additional_functions\aparc aseg-in-rawavg.nii'); %this is a specific version of 
% ft_read_atlas from SPM12, which was used in this project. Included in the additional functions folder because the new versino of fieldtrip is different
atlas.coordsys = 'mni';

S = struct2cell(EEG.chanlocs)'; 
elec_mni_frv.unit = 'mm';
elec_mni_frv.label = S(:,1);
elec_mni_frv.chanpos = cell2mat(S(:,2:4));

%remove channels without MNI coordinates
idx = cellfun(@(x) ~isempty(x), S(:,2));
EEG.chanlocs = EEG.chanlocs(idx);
EEG.data = EEG.data(idx,:);
EEG.nbchan = size(EEG.data,1);

clear myCoord;
for chani = 1:length(elec_mni_frv.chanpos)
    cfg            = [];
    cfg.roi        = elec_mni_frv.chanpos(chani,:);
    cfg.atlas      = atlas;
    cfg.inputcoord = 'mni';
    cfg.output     = 'label';
    cfg.maxqueryrange = 5;
    labels = ft_volumelookup(cfg, atlas);
    [~, indx] = max(labels.count); 
    % check https://www.fieldtriptoolbox.org/faq/how_can_i_determine_the_anatomical_label_of_a_source/
    myCoord(chani,:) = labels.name(indx);
end
[EEG.chanlocs.location] = myCoord{:}; 
coord = string(myCoord);



% % % % % This version below provides slightly different results because it use spm_read_vols function, which is similar to more recent ft_read_atlas versions
% % % % S = struct2cell(EEG.chanlocs)'; 
% % % % chanpos = cell2mat(S(:,2:4));
% % % % addpath('D:\Documents\MATLAB\spm12\');
% % % % V_aparc                 = spm_vol('D:\Documents\GitHub\WM\additional_functions\aparc aseg-in-rawavg.nii');
% % % % %V_aparc                 = spm_vol('D:\Documents\GitHub\WM\additional_functions\aparc+aseg.mgz');
% % % % 
% % % % [Y_aparc, XYZ_aparc]    = spm_read_vols(V_aparc);
% % % % load('fs_regionvalue2regionlabel.mat');
% % % % fprintf('You can evaluate %d different regions using this freesurfer template...\n', size(fs_regionvalue2regionlabel, 1));
% % % % 
% % % % clear myCoord;
% % % % for chani = 1:length(chanpos)
% % % %     MNI_coors           = round(chanpos(chani, :)); 
% % % %     tmpidx              = find(XYZ_aparc(1, :) == MNI_coors(1) & XYZ_aparc(2, :) == MNI_coors(2) & XYZ_aparc(3, :) == MNI_coors(3));
% % % %     this_region_value   = Y_aparc(tmpidx);
% % % %     tmpidx              = find(cell2mat(fs_regionvalue2regionlabel(:, 1)) == this_region_value);
% % % %     fs_region_name      = fs_regionvalue2regionlabel{tmpidx, 2};
% % % %     myCoord{chani}    = fs_region_name; 
% % % % end
% % % % 
% % % % [EEG.chanlocs.location] = myCoord{:}; 
% % % % coord = string(myCoord);


switch region
    
    case 'vvs'
        regions2use = { 'ctx-lh-inferiortemporal' , 'ctx-rh-inferiortemporal' ...
                        'ctx-lh-middletemporal', 'ctx-rh-middletemporal', ...
                        'ctx-lh-superiortemporal', 'ctx-rh-superiortemporal', ...
                        'ctx-lh-bankssts', 'ctx-rh-bankssts', ...                    
                        'ctx-lh-fusiform' , 'ctx-rh-fusiform', ...
                        'ctx-lh-temporalpole', 'ctx-rh-temporalpole', ...
                        'ctx-lh-inferiorparietal', 'ctx-rh-inferiorparietal',...
                        'ctx-lh-lateraloccipital', 'ctx-rh-lateraloccipital',...
                        'ctx-lh-lingual', 'ctx-lh-lingual', ...
                        'ctx-lh-parahippocampal', 'ctx-rh-parahippocampal',...
                        'ctx-lh-cuneus', 'ctx-rh-cuneus',...
                        'ctx-lh-pericalcarine', 'ctx-rh-pericalcarine',...
                        'ctx-lh-entorhinal', 'ctx-rh-entorhinal'};          
        regions2use = regions2use'
        
    case 'pfc'
        regions2use = {
        'ctx-lh-caudalanteriorcingulate', 'ctx-lh-caudalmiddlefrontal', 'ctx-lh-lateralorbitofrontal' ,...
        'ctx-lh-medialorbitofrontal', 'ctx-lh-parsopercularis' ,'ctx-lh-parsorbitalis' ,'ctx-lh-parstriangularis', ...
        'ctx-lh-rostralanteriorcingulate', 'ctx-lh-rostralmiddlefrontal', 'ctx-lh-superiorfrontal', ...
        'ctx-lh-frontalpole', ...
        
        
        'ctx-rh-caudalanteriorcingulate', 'ctx-rh-caudalmiddlefrontal', 'ctx-rh-lateralorbitofrontal' ,...
        'ctx-rh-medialorbitofrontal', 'ctx-rh-parsopercularis' ,'ctx-rh-parsorbitalis' ,'ctx-rh-parstriangularis', ...
        'ctx-rh-rostralanteriorcingulate', 'ctx-rh-rostralmiddlefrontal', 'ctx-rh-superiorfrontal', ...
        'ctx-rh-frontalpole', ...
        };
        regions2use = regions2use'
        
    case 'all'
        regions2use = {
            'ctx-lh-bankssts', 'ctx-lh-caudalanteriorcingulate', 'ctx-lh-caudalmiddlefrontal', ...
            'ctx-lh-corpuscallosum', 'ctx-lh-cuneus' , 'ctx-lh-entorhinal' ,'ctx-lh-fusiform', 'ctx-lh-inferiorparietal' ,...
            'ctx-lh-inferiortemporal' ,'ctx-lh-isthmuscingulate' ,'ctx-lh-lateraloccipital' ,'ctx-lh-lateralorbitofrontal' ,...
            'ctx-lh-lingual', 'ctx-lh-medialorbitofrontal', 'ctx-lh-middletemporal', 'ctx-lh-parahippocampal', 'ctx-lh-paracentral',...
            'ctx-lh-parsopercularis' ,'ctx-lh-parsorbitalis' ,'ctx-lh-parstriangularis', 'ctx-lh-pericalcarine' ,'ctx-lh-postcentral',...
            'ctx-lh-posteriorcingulate', 'ctx-lh-precentral', 'ctx-lh-precuneus', 'ctx-lh-rostralanteriorcingulate' ,...
            'ctx-lh-rostralmiddlefrontal', 'ctx-lh-superiorfrontal', 'ctx-lh-superiorparietal', 'ctx-lh-superiortemporal', ...
            'ctx-lh-supramarginal', 'ctx-lh-frontalpole', 'ctx-lh-temporalpole', 'ctx-lh-transversetemporal', 'ctx-lh-insula', ...
            'ctx-rh-bankssts', 'ctx-rh-caudalanteriorcingulate', 'ctx-rh-caudalmiddlefrontal', 'ctx-rh-corpuscallosum', ...
            'ctx-rh-cuneus', 'ctx-rh-entorhinal', 'ctx-rh-fusiform', 'ctx-rh-inferiorparietal', 'ctx-rh-inferiortemporal', 'ctx-rh-isthmuscingulate', ...
            'ctx-rh-lateraloccipital', 'ctx-rh-lateralorbitofrontal', 'ctx-rh-lingual', 'ctx-rh-medialorbitofrontal', 'ctx-rh-middletemporal', ...
            'ctx-rh-parahippocampal', 'ctx-rh-paracentral', 'ctx-rh-parsopercularis', 'ctx-rh-parsorbitalis', 'ctx-rh-parstriangularis', ...
            'ctx-rh-pericalcarine', 'ctx-rh-postcentral', 'ctx-rh-posteriorcingulate', 'ctx-rh-precentral', 'ctx-rh-precuneus',  ...
            'ctx-rh-rostralanteriorcingulate', 'ctx-rh-rostralmiddlefrontal', 'ctx-rh-superiorfrontal', 'ctx-rh-superiorparietal', ...
            'ctx-rh-superiortemporal', 'ctx-rh-supramarginal', 'ctx-rh-frontalpole', 'ctx-rh-temporalpole', 'ctx-rh-transversetemporal', ...
            'ctx-rh-insula', ... 
            'Left-Thalamus','Left-Thalamus-Proper', 'Left-Caudate', 'Left-Putamen', 'Left-Pallidum', 'Left-Hippocampus', 'Left-Amygdala', 'Left-Insula', ...
            'Left-Operculum', 'Left-Accumbens-area', 'Left-Substancia-Nigra', ...
            'Right-Thalamus','Right-Thalamus-Proper', 'Right-Caudate', 'Right-Putamen', 'Right-Pallidum', 'Right-Hippocampus', 'Right-Amygdala', 'Right-Insula', ...
            'Right-Operculum', 'Right-Accumbens-area', 'Right-Substancia-Nigra'
        };
        regions2use = regions2use'
        
    case 'hipp'
         regions2use = {
            'Left-Hippocampus', 'Right-Hippocampus'
        };
        regions2use = regions2use'
        
        


end



[id2] = ismember(coord,regions2use);
if isempty(id2)
   disp('No electrodes in the regions indicated'); 
end

EEG.chanlocs = EEG.chanlocs(id2);
EEG.data = EEG.data(id2,:);
EEG.nbchan = size(EEG.data,1);

end


