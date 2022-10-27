
function [act_CH act_FR] = load_rnn_activ()

disp (['Loading RNN activations ... ']);


% % % % % load RNN activations activations first
% get current path
currentFolder = pwd; 
cd ('/Users/danielpacheco/Documents/iEEG_data_analysis/WM/rnn_activations/RDMs/spearman/')
%cd ('/Users/danielpacheco/Documents/iEEG_data_analysis/WM/rnn_activations/RDMs/pearson/')

cd CH
sublist = dir('*_CH.mat');
sublist = {sublist.name}; sublist = sort(sublist');

tic

for layi=1:length(sublist)
    load(sublist{layi});
    fs2= squeeze(a(1, :, :)); 
    %fs2= a; 
    act_CH(layi,:,:) = fs2; 
            
end

cd ..

cd FR
sublist = dir('*_FR.mat');
sublist = {sublist.name}; sublist = sublist';

for layi=1:length(sublist)
    load(sublist{layi});
    fs2= squeeze(a(1, :, :)); 
    %fs2= a; 
    act_FR(layi,:,:) = fs2; 
            
end            

cd ..
cd (currentFolder)

end