
    f2t = strsplit(f2sav, '_'); 
    if strcmp(f2t{2}, 'vvs') 
        subj_ch_fr = 17; 
        nSubj = 28;
        
    elseif strcmp(f2t{2}, 'pfc') 
        subj_ch_fr = 7;
        nSubj = 16;
        
     elseif strcmp(f2t{2}, 'hipp') 
        subj_ch_fr = 8;
        
    end

    if strcmp(f2t{3}, 'E') 
        it={'all'}; 
        pp1 = double(string(strsplit(f2t{7}, '-')));
        resTime = pp1(1):pp1(2);
    elseif strcmp(f2t{3}, 'M') 
        it={'7'}; 
        pp1 = double(string(strsplit(f2t{7}, '-')));
        resTime = pp1(1):pp1(2);
        
    elseif strcmp(f2t{3}, 'MA') 
        it={'8'}; 
        pp1 = double(string(strsplit(f2t{7}, '-')));
        resTime = pp1(1):pp1(2);
    elseif strcmp(f2t{3}, 'Probe') 
        it={'9'}; 
        pp1 = double(string(strsplit(f2t{7}, '-')));
        resTime = pp1(1):pp1(2); 
    end
    
    if strcmp(f2t{4}, 'Av') 
        avMeth= 'pow'; 
    else
        avMeth= 'none'; 
    end
    
    
   freqs2test          = [3:54]';
    

    if strcmp(f2t{1}, 'Alex')
        nLays = 8; 
    elseif strcmp(f2t{1}, 'RNN')
        nLays = 56; 
    else
        nLays = []; 
    end
        
    
    nTimes          = length(resTime)-4; %4 for the windows
    nFreqs          = length(freqs2test); 
    
        
    
    

