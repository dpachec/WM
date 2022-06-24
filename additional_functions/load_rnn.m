function[ACT] = load_rnn(lays2load, subji, subj_ch_fr);%load network if not loaded yet


currentFolder = pwd; 
cd ('D:\_WM\analysis\out_contrasts\data\RNN')

if subji < subj_ch_fr
    cd FR
    sublist = dir('*actRNN_FR.mat');
    sublist = {sublist.name}; sublist = sublist';

    for layi=1:length(lays2load)
        load(sublist{lays2load(layi)});
        fs2 = reshape(FrameStack, size(FrameStack, 1), []);
        ACT(layi,:,:) = corr(fs2', 'type', 's');    

    end            
else
    cd CH
    sublist = dir('*actRNN_CH.mat');
    sublist = {sublist.name}; sublist = sublist';

    tic

    for layi=1:length(lays2load)
        load(sublist{lays2load(layi)});
        fs2 = reshape(FrameStack, size(FrameStack, 1), []);
        ACT(layi,:,:) = corr(fs2', 'type', 's');    

    end

    cd ..
end
cd ..
cd (currentFolder)



end