
function [oneListPow] = extract_power_WM (cfg_contrasts , cfg)

    oneListTraces_c = cfg_contrasts.oneListTraces; 
    timeRes = cfg.timeRes;  
    DNN_analysis = cfg.DNNs; 
    if DNN_analysis
        period = cfg.period; 
    end

    
   % disp ('>> extracting power ... ');
    sr = 1000;
    data_ft = mat2ft(oneListTraces_c, sr);

    
    if timeRes == 0.01
         disp ('10ms');  
         lim_1   = 101; 
         lim_2   = 650; 
    end
    
    if timeRes == 0.1        
         lim_1   = 1; %
         lim_2   = 70; 
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
    cfg.output       = 'pow';
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
        
    if DNN_analysis
        if strcmp(period(1), 'E') & timeRes == 0.1
           oneListPow = oneListPow(:,:,:,16:40);
        elseif strcmp(period(1), 'M') & timeRes == 0.1
            oneListPow = oneListPow(:,:,:,16:65);
        end
    end


end

