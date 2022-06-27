
function data = customReadDatstoreImage(filename, network)

    onState = warning('off', 'backtrace'); 
    c = onCleanup(@() warning(onState)); 
    data = imread(filename); % added lines: 
    data = data(:,:,min(1:3, end)); 
    data = imresize(data,[227 227]);

    
    
end