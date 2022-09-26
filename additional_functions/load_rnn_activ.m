
function [act_CH act_FR] = load_rnn_activ()

disp (['Loading RNN activations ... ']);


% % % % % load RNN activations activations first
% get current path
currentFolder = pwd; 
cd ('/Users/danielpacheco/Documents/iEEG_data_analysis/WM/rnn_activations')

cd CH
sublist = dir('*actRNN_CH.mat');
sublist = {sublist.name}; sublist = sublist';

tic

for layi=1:length(sublist)
    load(sublist{layi});
    fs2 = reshape(FrameStack, size(FrameStack, 1), []);
    act_CH(layi,:,:) = corr(fs2', 'type', 's');    
            
end

cd ..

cd FR
sublist = dir('*actRNN_FR.mat');
sublist = {sublist.name}; sublist = sublist';

for layi=1:length(sublist)
    load(sublist{layi});
    fs2 = reshape(FrameStack, size(FrameStack, 1), []);
    act_FR(layi,:,:) = corr(fs2', 'type', 's');    
            
end            

cd ..
cd (currentFolder)

end