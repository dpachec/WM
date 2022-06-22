
function [act_CH act_FR] = load_alex_activ(ver)

if strcmp(ver, 'imageNet')
    currentFolder = pwd; 
    net = alexnet;
    imgPath = 'D:\_WM\stimuli_experiment\china';
    imageDS = imageDatastore(imgPath); 
    imageDS.ReadFcn = @customReadDatstoreImage;

    clear act
    act{1} = activations(net,imageDS,'conv1');
    act{2} = activations(net,imageDS,'conv2');
    act{3} = activations(net,imageDS,'conv3');
    act{4} = activations(net,imageDS,'conv4');
    act{5} = activations(net,imageDS,'conv5');
    act{6} = activations(net,imageDS,'fc6');
    act{7} = activations(net,imageDS,'fc7');
    act{8} = activations(net,imageDS,'fc8');


    clear act_CH
    for layi = 1:length(act)
       clear fs1 fs2
       fs1 = act{layi};
       fs1 = permute(fs1, [4 1 2 3]);
       fs2 = reshape(fs1, size(fs1, 1), []);
       act_CH(layi, :, :) = corr(fs2', 'type', 's'); 
    end


    clearvars -except act_CH net currentFolder
    imgPath = 'D:\_WM\stimuli_experiment\freiburg';
    imageDS = imageDatastore(imgPath); 
    imageDS.ReadFcn = @customReadDatstoreImage;

    act{1} = activations(net,imageDS,'conv1');
    act{2} = activations(net,imageDS,'conv2');
    act{3} = activations(net,imageDS,'conv3');
    act{4} = activations(net,imageDS,'conv4');
    act{5} = activations(net,imageDS,'conv5');
    act{6} = activations(net,imageDS,'fc6');
    act{7} = activations(net,imageDS,'fc7');
    act{8} = activations(net,imageDS,'fc8');


    clear act_FR
    for layi = 1:length(act)
       clear fs1 fs2
       fs1 = act{layi};
       fs1 = permute(fs1, [4 1 2 3]);
       fs2 = reshape(fs1, size(fs1, 1), []);
       act_FR(layi, :, :) = corr(fs2', 'type', 's'); 
    end
    cd (currentFolder)
      


elseif  strcmp(ver, 'ecoset')
    
    currentFolder = pwd; 
    cd D:\_WM\analysis\out_contrasts\data\Alexnet\CH\training_seed_05\input_image_sets\activations
    sublist = dir('*_ecoset.mat');
    sublist = {sublist.name}';
    disp (['measurements -> ' num2str(length(sublist))]);
 
    count = 1; 
    for layi=1:8
        for imi = 1:60
            act = load(sublist{count});
            act = struct2cell(act); act = act{:};
            allV{imi} = act(:); 
            count = count+1;
        end
        allV1 = cell2mat(allV);
        act_CH(layi, :, :) = corr(allV1, 'type', 's');
    end
    
        
    cd D:\_WM\analysis\out_contrasts\data\Alexnet\FR\training_seed_05\input_image_sets\activations
    sublist = dir('*_ecoset.mat');
    sublist = {sublist.name}';
    disp (['measurements -> ' num2str(length(sublist))]);
 
    clear allV
    count = 1; 
    for layi=1:8
        for imi = 1:60
            act = load(sublist{count});
            act = struct2cell(act); act = act{:};
            allV{imi} = act(:); 
            count = count+1;
        end
        allV1 = cell2mat(allV);
        act_FR(layi, :, :) = corr(allV1, 'type', 's');
    end
    
    cd (currentFolder)
      
    
    
else
    
        
end 

end