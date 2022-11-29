
function [oneListPow] = extract_power_WM (cfg_contrasts , cfg)

    oneListTraces_c = cfg_contrasts.oneListTraces; 
    timeRes = cfg.timeRes;  
    period = cfg.period; 
    
   % disp ('>> extracting power ... ');
    sr = 1000;
    data_ft = mat2ft(oneListTraces_c, sr);

    
    if timeRes == 0.01
%        lim_1   = 1; %
%        lim_2   = 901; 
%         disp ('10ms');  
         lim_1   = 101; %
         lim_2   = 650; 
%         disp(['lim_1 = ' num2str(lim_1) ' lim_2 = ' num2str(lim_2)]);
    end
    
    if timeRes == 0.05 %total length 240 
         disp ('50ms');  
         lim_1   = 91 - 5; %substract only here but not to the baseline
         lim_2   = 190 - 5; %substract only here but not to the baseline
         disp(['lim_1 = ' num2str(lim_1) ' lim_2 = ' num2str(lim_2)]);
    end

    if timeRes == 0.1
        
         lim_1   = 1; %
         lim_2   = 70; 
        
%         disp ('100ms');  
%         %lim_1   = 14; %
%         %lim_2   = 70; %for the locked to cue 
%         lim_1   = 1; % for the locked to probe
%         lim_2   = 70; %
%         
%         
%         disp(['lim_1 = ' num2str(lim_1) ' lim_2 = ' num2str(lim_2)]);
    end

    %fieldtrip analysis
    cfg              = [];
    cfg.method       = 'wavelet'; %%'mtmconvol'; % 'wavelet'; %
    cfg.width        = linspace(3, 6, 29);
    cfg.output       = 'pow';
    cfg.foi          = [1:1:29];   % analysis 2 to X Hz in steps of 2 Hz 
    %cfg.taper        = ['hanning']; %'dpss' %> slepian sequence of tapers
    cfg.pad          = 'nextpow2';
    if strcmp (timeRes, 'all')
        cfg.toi          = 'all'; % takes as reference the number of time windows defined above
    else
        %vector 1 x numtoi, the times on which the analysis windows should be centered (in seconds)
        cfg.toi          = 0:timeRes:(length(oneListTraces_c)/sr);
    end
    cfg.keeptrials   = 'yes'; % keep individual trials or average
    cfg.showcallinfo = 'no';% no log console
    cfg.feedback     = 'none'; 
    tf_data_L          = ft_freqanalysis(cfg, data_ft);
    dataL = tf_data_L.powspctrm;



    cfg              = [];
    cfg.method       = 'wavelet'; %'mtmconvol'; % 'wavelet'; %
    cfg.width        = linspace(6, 12, 25);
    %cfg.width        = linspace(6, 12, 121);
    cfg.output       = 'pow';
    %cfg.foi          = [30:1:150]   % analysis 2 to X Hz in steps of 2 Hz 
    cfg.foi          = [30:5:150];   % analysis 2 to X Hz in steps of 2 Hz 
    cfg.pad          = 'nextpow2';
    if strcmp (timeRes, 'all')
        cfg.toi          = 'all'; 
    else
        cfg.toi          = 0:timeRes:(length(oneListTraces_c)/sr); % 50ms; 
    end
    cfg.showcallinfo = 'no';% no log console
    cfg.feedback     = 'none'; 
    cfg.keeptrials   = 'yes'; % keep individual trials, if not, it makes an average 
    tf_data_H          = ft_freqanalysis(cfg, data_ft);
    dataH = tf_data_H.powspctrm;




    dataLH = cat (3, dataL, dataH);


    %cut edge artefacts from -6-6 to finalCut
    if ~strcmp(timeRes, 'all')
        %disp('before cut'); size(dataLH)
        dataLH = dataLH(:,:,:,lim_1:lim_2);
        %disp('after cut'); size(dataLH)
    end

    oneListPow = dataLH;
        
     if strcmp(period(1), 'E')
        oneListPow = oneListPow(:,:,:,16:40);
     elseif strcmp(period(1), 'M')
        oneListPow = oneListPow(:,:,:,16:65);
    end




end

