
function [act_CH act_FR] = load_alex_activ(cfg, sessi, subj_ch_fr, paths);

lays2load = cfg.lays2load;


    lays = {'conv1' 'conv2' 'conv3' 'conv4' 'conv5' 'fc6' 'fc7' 'fc8'};
    %lays = {'relu1' 'relu2' 'relu3' 'relu4' 'relu5' 'relu6' 'relu7' 'fc8'};
    currentFolder = pwd; 
    net = alexnet; %equivalent to net = alexnet('Weights','imagenet')
    
    if sessi < subj_ch_fr
        imgPath = [paths.stim 'freiburg'];
        imageDS = imageDatastore(imgPath); 
        imageDS.ReadFcn = @customReadDatstoreImage;
        
        for layi=1:length(lays2load)
            fs1 = activations(net,imageDS,lays{lays2load(layi)});
            fs1 = permute(fs1, [4 1 2 3]);
            fs2 = reshape(fs1, size(fs1, 1), []);
            act_CH(layi, :, :) = corr(fs2', 'type', 's'); 
            ACT(layi,:,:) = corr(fs2', 'type', 's');    
        end            
    else
        imgPath = [paths.stim 'china'];
        imageDS = imageDatastore(imgPath); 
        imageDS.ReadFcn = @customReadDatstoreImage;
        for layi=1:length(lays2load)
            fs1 = activations(net,imageDS,lays{lays2load(layi)});
            fs1 = permute(fs1, [4 1 2 3]);
            fs2 = reshape(fs1, size(fs1, 1), []);
            act_CH(layi, :, :) = corr(fs2', 'type', 's'); 
            ACT(layi,:,:) = corr(fs2', 'type', 's');    
        end            
        
    end

    cd (currentFolder)
      
end