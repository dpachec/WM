%%
function[ACT] = load_blnext(cfg, subji, paths, oneListIDs);%load network if not loaded yet

    net2load = cfg.net2load; 
    lays2load = cfg.lays2load; 
    brainROI  = cfg.brainROI; 
    period = cfg.period; 
    

    cd (paths.activations.resNet)

    
    nLays = 1; % for now
    if strcmp(net2load, 'Res18-2-3-0')
        all_act = load('Features_resnet18_bl_8steps-layer3[0].record-2samples.mat');
    elseif strcmp(net2load, 'Res18-2-3-1')
        all_act = load('Features_resnet18_bl_8steps-layer3[1].record-2samples.mat');
    elseif strcmp(net2load, 'Res18-2-4-0')
        all_act = load('Features_resnet18_bl_8steps-layer4[0].record-2samples.mat');
    elseif strcmp(net2load, 'Res18-2-4-1')
        all_act = load('Features_resnet18_bl_8steps-layer4[1].record-2samples.mat');
    end


    if strcmp(net2load, 'Res18-4-3-0')
        all_act = load('Features_resnet18_bl_8steps-layer3[0].record-4samples.mat');
    elseif strcmp(net2load, 'Res18-4-3-1')
        all_act = load('Features_resnet18_bl_8steps-layer3[1].record-4samples.mat');
    elseif strcmp(net2load, 'Res18-4-4-0')
        all_act = load('Features_resnet18_bl_8steps-layer4[0].record-4samples.mat');
    elseif strcmp(net2load, 'Res18-4-4-1')
        all_act = load('Features_resnet18_bl_8steps-layer4[1].record-4samples.mat');
    end

    if strcmp(net2load, 'Res18-6-3-0')
        all_act = load('Features_resnet18_bl_8steps-layer3[0].record-6samples.mat');
    elseif strcmp(net2load, 'Res18-6-3-1')
        all_act = load('Features_resnet18_bl_8steps-layer3[1].record-6samples.mat');
    elseif strcmp(net2load, 'Res18-6-4-0')
        all_act = load('Features_resnet18_bl_8steps-layer4[0].record-6samples.mat');
    elseif strcmp(net2load, 'Res18-6-4-1')
        all_act = load('Features_resnet18_bl_8steps-layer4[1].record-6samples.mat');
    end
       

    if strcmp(net2load, 'Res18-8-3-0')
        all_act = load('Features_resnet18_bl_8steps-layer3[0].record-8samples.mat');
    elseif strcmp(net2load, 'Res18-8-3-1')
        all_act = load('Features_resnet18_bl_8steps-layer3[1].record-8samples.mat');
    elseif strcmp(net2load, 'Res18-8-4-0')
        all_act = load('Features_resnet18_bl_8steps-layer4[0].record-8samples.mat');
    elseif strcmp(net2load, 'Res18-8-4-1')
        all_act = load('Features_resnet18_bl_8steps-layer4[1].record-8samples.mat');
    end
    

    if strcmp(net2load(1:7), 'Res18-2')
        seqEnd = 6; 
    elseif strcmp(net2load(1:7), 'Res18-4')
        seqEnd = 12; 
    elseif strcmp(net2load(1:7), 'Res18-6')
        seqEnd = 18; 
    elseif strcmp(net2load(1:7), 'Res18-8')
        seqEnd = 24;
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
    
        % % % % acts now contains for each sequence the activations after each image presentation (seqLen x netTime x units)
        
        for evi = 1:size(acts, 2)
            out=squeeze(acts(:, evi, :));
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
            % take first second or third position in the item sequence
            if strcmp(cfg.BLIT, 'F')
                items2Cons = cellfun(@(x) x(5:7), items2Cons, 'un', 0 ); % 5:7 = 1 ; 9:11 = 2; 13:15 = 3
            elseif strcmp(cfg.BLIT, 'S')
                items2Cons = cellfun(@(x) x(9:11), items2Cons, 'un', 0 ); % 5:7 = 1 ; 9:11 = 2; 13:15 = 3
            elseif strcmp(cfg.BLIT, 'T')
                items2Cons = cellfun(@(x) x(13:15), items2Cons, 'un', 0 ); % 5:7 = 1 ; 9:11 = 2; 13:15 = 3
            end
            id2u = contains(items2Cons, itH); 
            j2check = all_act_H(id2u, 1);
            
            if ~isempty(id2u)
                allMatchingImages = all_act_H(id2u, 2); 
                %currentImage = allMatchingImages{1}; 
                for matchiMi = 1:length(allMatchingImages)
                    currentImage = allMatchingImages{matchiMi}; 
                    for evi = 1:length(currentImage.events) % This loops over time 
                        allMacts(matchiMi, itemi, evi, :) = currentImage.predictions(evi,:,:);
                    end
                end
            end
        end

        for allmi = 1:size(allMacts, 1)
            acts = squeeze(allMacts(allmi, :, :, :)); 
            for evi = 1:size(acts, 2)
                out=squeeze(acts(:, evi, :));
                ids2rem = all(out==0,2); % remove items not present 
                out(ids2rem,:) = nan;
            
                ACT_prev = corr(out', 'type', 's');
                ACT(allmi, evi, :, :) = ACT_prev;
                
            end
        end


     ACT = squeeze(mean(ACT, 'omitnan')); 

     end


   cd (paths.github)

 end






 



    






