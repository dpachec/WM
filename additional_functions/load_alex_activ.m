
function [act_CH act_FR] = load_alex_activ(lays2load, subji, subj_ch_fr, path);%load network if not loaded yet

    lays = {'conv1' 'conv2' 'conv3' 'conv4' 'conv5' 'fc6' 'fc7' 'fc8'};
    currentFolder = pwd; 
    net = alexnet;
    
    if subji < subj_ch_fr
        imgPath = [path 'freiburg'];
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
        imgPath = [path 'china'];
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