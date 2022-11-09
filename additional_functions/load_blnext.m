%%
function[ACT ids2rem] = load_blnext(net2load, lays2load, subji, subj_ch_fr, paths, brainROI, oneListIDs, period);%load network if not loaded yet

    cd (paths.multi_item_activations)

    
    nLays = 1; % for now
    if strcmp(net2load, 'BLnext2')
        all_act = load('Features_BLnext-2samples.mat');
        %all_act{2} = load('Features_BLnext-ReLU_Layer_6-2samples.mat');
        seqEnd = 6; 
    elseif strcmp(net2load, 'BLnext4')

        all_act = load('Features_BLnext-4samples.mat');
        %all_act{2} = load('Features_BLnext-ReLU_Layer_6-4samples.mat');
        seqEnd = 12; 
    elseif strcmp(net2load, 'BLnext8')
        all_act = load('Features_BLnext-8samples.mat');
        %all_act{2} = load('Features_BLnext-ReLU_Layer_6-8samples.mat');
        seqEnd = 24; 
    elseif strcmp(net2load, 'BLnext12')
        all_act = load('Features_BLnext-12samples.mat');
        %all_act{1} = load('Features_BLnext-ReLU_Layer_6-12samples.mat');
        %all_act{2} = load('Features_BLnext-ReLU_Layer_6-12samples.mat');
        seqEnd = 36; 
    end



     if strcmp(period(1:2), 'E4') | strcmp(period(1), 'M')
        id3 = cellfun(@(x) strsplit(x), oneListIDs, 'un', 0);
        id4 = cellfun(@(x) double(string(x(:, 13:15))), id3, 'un', 0);
        seqSub = cat(1, id4{:});

        f = fieldnames(all_act);
        all_act_H = [f struct2cell(all_act)];
        for seqi = 1:length(seqSub) % for each sequence presented to a subject
            seqStr = ['Seq' join(string(seqSub(seqi,:)), '_')]; seqStr = join(seqStr, '_'); %this is to match Lynns format 
            id2u = strmatch(seqStr, all_act_H(:,1)); 
            %actSub{subji,:}(seqi, :) = squeeze(all_act{id2u, 2}.predictions(6,:,:));
            if ~isempty(id2u)
                %disp (['L > ' length(all_act_H{id2u, 2}.events)]);
                for evi = 1:length(all_act_H{id2u, 2}.events) % This loops over time images were presented to the network
                    acts(seqi,evi, :) = all_act_H{id2u, 2}.predictions(evi,:,:);
                end
            end
        end
    
        % % % % acts now contains for each sequence the data presented in time
        % to the network (seqLen x netTime x units)
        
        for evi = 1:size(acts, 2)
            out=squeeze(acts(:, evi, :));
            ids2rem = all(out==0,2); % remove sequences not present in Lynns activations 
            out(ids2rem,:) = [];
            
            ACT_prev = corr(out', 'type', 's');
            ACT(evi, :, :) = ACT_prev;
        end
        
    
    
     elseif strcmp(period(1), 'E')
        id3 = cellfun(@(x) strsplit(x), oneListIDs, 'un', 0);
        id4 = cat(1, id3{:});
        itIDs = id4(:, 3);

        f = fieldnames(all_act);
        all_act_H = [f struct2cell(all_act)];
        for itemi = 1:length(itIDs) 
            
            itH = itIDs(itemi); 
            items2Cons = all_act_H(:,1); 
            items2Cons = cellfun(@(x) x(5:7), items2Cons, 'un', 0 );
            id2u = contains(items2Cons, itH); 
            j2check = all_act_H(id2u, 1);
            
            if ~isempty(id2u)
                allMatchingImages = all_act_H(id2u, 2); 
                currentImage = allMatchingImages{1}; 
                for evi = 1:length(currentImage.events) % This loops over time 
                    acts(itemi,evi, :) = currentImage.predictions(evi,:,:);
                end
            end
        end

        for evi = 1:size(acts, 2)
            out=squeeze(acts(:, evi, :));
            ids2rem = all(out==0,2); % remove items not present 
            out(ids2rem,:) = [];
            
            ACT_prev = corr(out', 'type', 's');
            ACT(evi, :, :) = ACT_prev;
        end

     end


   cd (paths.github)

 end

















 



    






