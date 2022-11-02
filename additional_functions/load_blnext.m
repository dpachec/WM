%%
function[ACT ids2rem] = load_blnext(net2load, lays2load, subji, subj_ch_fr, path, brainROI, oneListIDs);%load network if not loaded yet

    currentFolder = pwd; 
    cd (path)

%     % % % no need to load the activations (can be extracted from oneListIds_c)
%     if strcmp (brainROI, 'vvs')
%         cd image_sequences_vvs
%         sublist = dir('*.csv'); sublist = {sublist.name}'; 
%     elseif strcmp (brainROI, 'pfc')
%         cd image_sequences_pfc
%         sublist = dir('*.csv'); sublist = {sublist.name}'; 
%     end
%     all_seqs{subji,:} = csvread(sublist{subji}); 
%     cd .. 

    id3 = cellfun(@(x) strsplit(x), oneListIDs, 'un', 0);
    id4 = cellfun(@(x) double(string(x(:, 13:15))), id3, 'un', 0);
    seqSub = cat(1, id4{:});
    

    
    nLays = 1; % for now
    if strcmp(net2load, 'BLnext2')
        all_act{1} = load('Features_BLnext-2samples.mat');
        %all_act{2} = load('Features_BLnext-ReLU_Layer_6-2samples.mat');
        seqEnd = 6; 
    elseif strcmp(net2load, 'BLnext4')

        all_act{1} = load('Features_BLnext-4samples.mat');
        %all_act{2} = load('Features_BLnext-ReLU_Layer_6-4samples.mat');
        seqEnd = 12; 
    elseif strcmp(net2load, 'BLnext8')
        all_act{1} = load('Features_BLnext-8samples.mat');
        %all_act{2} = load('Features_BLnext-ReLU_Layer_6-8samples.mat');
        seqEnd = 24; 
    elseif strcmp(net2load, 'BLnext12')
        all_act{1} = load('Features_BLnext-12samples.mat');
        %all_act{1} = load('Features_BLnext-ReLU_Layer_6-12samples.mat');
        %all_act{2} = load('Features_BLnext-ReLU_Layer_6-12samples.mat');
        seqEnd = 36; 
    end

    
    for layi = 1:nLays
        clear acts
        f = fieldnames(all_act{layi});
        all_act_H = [f struct2cell(all_act{layi})];
        for seqi = 1:length(seqSub)
            seqStr = ['Seq' join(string(seqSub(seqi,:)), '_')]; seqStr = join(seqStr, '_');
            id2u = strmatch(seqStr, all_act_H(:,1)); 
            %actSub{subji,:}(seqi, :) = squeeze(all_act{id2u, 2}.predictions(6,:,:));
            if ~isempty(id2u)
                %disp (['L > ' length(all_act_H{id2u, 2}.events)]);
                for evi = 1:length(all_act_H{id2u, 2}.events)
                    acts(seqi,evi, :) = all_act_H{id2u, 2}.predictions(evi,:,:);
                end
            end
        end

        for evi = 1:size(acts, 2)
            out=squeeze(acts(:, evi, :));
            ids2rem = all(out==0,2); 
            out(ids2rem,:) = [];
            
            ACT_prev = corr(out', 'type', 's');
            ACT(evi, :, :) = ACT_prev;
        end
        
            
            
%         out=acts;
%         ids2rem = all(acts==0,2); 
%         out(ids2rem,:) = [];
%         
%         ACT_prev = corr(out', 'type', 's');
%         ACT(layi, :, :) = ACT_prev;
%         
        
        
        cd (currentFolder)

    
    end

    

% %     for subji = 1:32
% %         clear seqSub acts rdm
% %     
% %         seqSub = all_seqs{subji}; 
% %         for seqi = 1:length(seqSub)
% %             seqStr = ['Seq' join(string(seqSub(seqi,:)), '_')]; seqStr = join(seqStr, '_');
% %             id2u = strmatch(seqStr, all_act(:,1)); 
% %             %actSub{subji,:}(seqi, :) = squeeze(all_act{id2u, 2}.predictions(6,:,:));
% %             acts(seqi,:) = squeeze(all_act{id2u, 2}.predictions(seqEnd,:,:))';
% %             
% %     
% %         end
% %     
% %             rdm = corr(acts', 'type', 's');
% %     
% %             allRDMS{subji,:} = rdm; 
% %     
% %     
% %     end









