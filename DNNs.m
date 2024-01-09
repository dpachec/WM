%% Calculate from epoched raw traces
%% 
clear
%Network_ROI_layers_freqs_avRepet_avTimeFeatVect_freqResolv(0-1)_fitMode(0:noTrials; 1:Trials)_timeRes_win-width_mf

%f2sav = 'Alex_pfc_E123_[1-8]_3-54_1_0_1_0_.1_5_1'; 
%f2sav = 'CAT_pfc_M123_[1]_3-54_1_0_1_0_.1_5_1'; 
%f2sav = 'CORrtRELU_pfc_M123_[1-8]_3-54_1_0_1_0_.1_5_1'; 
%f2sav = 'BLNETe_pfc_M123_[8-8-56]_3-54_1_0_1_0_.1_5_1';
%f2sav = 'BLNETeBatchNorm_pfc_M123_[8-8-56]_3-54_1_0_1_0_.1_5_1';
%f2sav = 'AlexEco_pfc_E123_[1-8]_3-54_1_0_0_0_.1_5_1';
%
%f2sav = 'BLNETi_pfc_E123_[8-8-56]_3-8_1_0_0_0_.1_5_1';
%f2sav = 'CORrt_pfc_M123_[2-2-8]_3-54_1_0_1_0_.1_5_1'; 
%f2sav = 'Res18-2_pfc_MALL_[1]_3-54_0_0_1_0_.1_5_1'; 
%f2sav = 'CORrt_vvs_M123_[8]_3-54_1_0_1_0_.1_5_1';

%f2sav = 'BLNETe_pfc_M123NC_[8-8-56]_3-54_1_0_1_0_.1_5_1';
%f2sav = 'Res34-2-3-0_vvs_MALL_[1]_3-54_0_0_1_0_.1_5_1'; 
%f2sav =  'AE-t10_hipp_M123_[1-18]_3-54_1_0_1_0_.1_5_1'; 

f2sav = 'CAT_vvs_E123_[1]_3-54_1_0_1_0_.1_5_1';
%f2sav = 'ITM_vvs_E123_[1]_3-54_0_0_1_0_.1_5_1';


cfg = getParams(f2sav);
cfg.DNN_analysis = 1; 
paths = load_paths_WM(cfg.brainROI, cfg.net2load);
filelistSess = getFilesWM(paths.traces);

t1 = datetime; 

for sessi= 1:length(filelistSess) %this one starts at 1 and not at 3
    disp(['File > ' num2str(sessi)]);
    load([paths.traces filelistSess{sessi}]);   

    cfg_contrasts               = getIdsWM(cfg.period, cfg_contrasts);

    if length(cfg_contrasts.oneListIds) > 1
        cfg_contrasts.oneListPow    = extract_power_WM (cfg_contrasts, cfg); % 
    
        cfg_contrasts               = normalize_WM(cfg_contrasts, 1, 'sess', []);
    
        cfg_contrasts               = average_repetitions(cfg, cfg_contrasts);
        
        neuralRDMs                  = createNeuralRDMs(cfg, cfg_contrasts);
        networkRDMs                 = createNetworkRDMs(cfg, cfg_contrasts, sessi, paths);

        % just to test with shuffled labels
        sC = size(networkRDMs, 2);
        ids = randperm(sC);
        networkRDMs = networkRDMs(:, ids, ids); 

    
        nnFit{sessi,1}              = fitModel_WM(neuralRDMs, networkRDMs, cfg.fitMode); 
        nnFit{sessi,2}              = cfg_contrasts.oneListIds; 
    end
end

mkdir ([paths.results.DNNs]);
save([paths.results.DNNs f2sav '.mat'], 'nnFit');

t2 = datetime; 
etime(datevec(t2), datevec(t1))


%% IN LOOP 
clear, clc
%Network_ROI_EoM_layers_freqs_avRepet_avTimeFeatVect_freqResolv(0-1)__fitMode(0:noTrials; 1:Trials)__timeRes__win-width__mf
    
listF2sav = {

'Alex_vvs_E123CI_[1-8]_3-54_0_0_1_0_.1_5_1';
'Alex_vvs_E123CC_[1-8]_3-54_0_0_1_0_.1_5_1';
'Alex_vvs_E123II_[1-8]_3-54_0_0_1_0_.1_5_1';
'Alex_vvs_E123IC_[1-8]_3-54_0_0_1_0_.1_5_1';

'BLNETi_vvs_E123CI_[8-8-56]_3-54_0_0_1_0_.1_5_1';
'BLNETi_vvs_E123CC_[8-8-56]_3-54_0_0_1_0_.1_5_1';
'BLNETi_vvs_E123II_[8-8-56]_3-54_0_0_1_0_.1_5_1';
'BLNETi_vvs_E123IC_[8-8-56]_3-54_0_0_1_0_.1_5_1';

'CORrt_vvs_E123CI_[2-2-8]_3-54_0_0_1_0_.1_5_1';
'CORrt_vvs_E123CC_[2-2-8]_3-54_0_0_1_0_.1_5_1';
'CORrt_vvs_E123II_[2-2-8]_3-54_0_0_1_0_.1_5_1';
'CORrt_vvs_E123IC_[2-2-8]_3-54_0_0_1_0_.1_5_1';

};   


for listi = 1:length(listF2sav)
    disp(['File > ' num2str(listi) '      ' listF2sav{listi}]);
    clearvars -except listF2sav listi 
        
    f2sav       = listF2sav{listi}; 
    cfg = getParams(f2sav);
    cfg.DNN_analysis = 1; 
    paths = load_paths_WM(cfg.brainROI, cfg.net2load);
    filelistSess = getFiles(paths.traces);
    
    for sessi= 1:length(filelistSess) %this one starts at 1 and not at 3
        %disp(['File > ' num2str(sessi)]);
        load([paths.traces filelistSess{sessi}]);   
    
        cfg_contrasts               = getIdsWM(cfg.period, cfg_contrasts);
    
        if length(cfg_contrasts.oneListIds) > 1
            cfg_contrasts.oneListPow    = extract_power_WM (cfg_contrasts, cfg); % 
        
            cfg_contrasts               = normalize_WM(cfg_contrasts, 1, 'sess', []);
        
            cfg_contrasts               = average_repetitions(cfg, cfg_contrasts);
            
            neuralRDMs                  = createNeuralRDMs(cfg, cfg_contrasts);
            networkRDMs                 = createNetworkRDMs(cfg, cfg_contrasts, sessi, paths);
        
            nnFit{sessi,1}              = fitModel_WM(neuralRDMs, networkRDMs, cfg.fitMode); 
            %nnFit{sessi,1}              = fitModelPartialCorrelation(neuralRDMs, networkRDMs, cfg.fitMode); 
            nnFit{sessi,2}              = cfg_contrasts.oneListIds; 
        end

    end

    
    save([paths.results.DNNs f2sav '.mat'], 'nnFit');

end






%% PERMUTATIONS IN LOOP
%Network_ROI_ER_layers_freqs_avRepet_avTFV_fRes(0-1)_fitMode(0:noTrials; 1:Trials)__timeRes__win__mf

clear
nPerm = 100;

listF2sav = {


'ITM_vvs_E123CAT1_[1]_3-54_0_0_1_0_.1_5_1';
'ITM_vvs_E123CAT2_[1]_3-54_0_0_1_0_.1_5_1';
'ITM_vvs_E123CAT3_[1]_3-54_0_0_1_0_.1_5_1';
'ITM_vvs_E123CAT4_[1]_3-54_0_0_1_0_.1_5_1';
'ITM_vvs_E123CAT5_[1]_3-54_0_0_1_0_.1_5_1';
'ITM_vvs_E123CAT6_[1]_3-54_0_0_1_0_.1_5_1';

'ITM_pfc_E123CAT1_[1]_3-54_0_0_1_0_.1_5_1';
'ITM_pfc_E123CAT2_[1]_3-54_0_0_1_0_.1_5_1';
'ITM_pfc_E123CAT3_[1]_3-54_0_0_1_0_.1_5_1';
'ITM_pfc_E123CAT4_[1]_3-54_0_0_1_0_.1_5_1';
'ITM_pfc_E123CAT5_[1]_3-54_0_0_1_0_.1_5_1';
'ITM_pfc_E123CAT6_[1]_3-54_0_0_1_0_.1_5_1';
% 
'ITM_vvs_M123CAT1_[1]_3-54_0_0_1_0_.1_5_1';
'ITM_vvs_M123CAT2_[1]_3-54_0_0_1_0_.1_5_1';
'ITM_vvs_M123CAT3_[1]_3-54_0_0_1_0_.1_5_1';
'ITM_vvs_M123CAT4_[1]_3-54_0_0_1_0_.1_5_1';
'ITM_vvs_M123CAT5_[1]_3-54_0_0_1_0_.1_5_1';
'ITM_vvs_M123CAT6_[1]_3-54_0_0_1_0_.1_5_1';
% 
'ITM_pfc_M123CAT1_[1]_3-54_0_0_1_0_.1_5_1';
'ITM_pfc_M123CAT2_[1]_3-54_0_0_1_0_.1_5_1';
'ITM_pfc_M123CAT3_[1]_3-54_0_0_1_0_.1_5_1';
'ITM_pfc_M123CAT4_[1]_3-54_0_0_1_0_.1_5_1';
'ITM_pfc_M123CAT5_[1]_3-54_0_0_1_0_.1_5_1';
'ITM_pfc_M123CAT6_[1]_3-54_0_0_1_0_.1_5_1';

'CAT_vvs_M11_[1]_3-54_1_0_1_0_.1_5_1';
'CAT_pfc_M11_[1]_3-54_1_0_1_0_.1_5_1';
'CAT_vvs_M12_[1]_3-54_1_0_1_0_.1_5_1';
'CAT_pfc_M12_[1]_3-54_1_0_1_0_.1_5_1';
'CAT_vvs_M13_[1]_3-54_1_0_1_0_.1_5_1';
'CAT_pfc_M13_[1]_3-54_1_0_1_0_.1_5_1';

'CAT_vvs_M11_[1]_3-54_0_0_1_0_.1_5_1';
'CAT_pfc_M11_[1]_3-54_0_0_1_0_.1_5_1';
'CAT_vvs_M12_[1]_3-54_0_0_1_0_.1_5_1';
'CAT_pfc_M12_[1]_3-54_0_0_1_0_.1_5_1';
'CAT_vvs_M13_[1]_3-54_0_0_1_0_.1_5_1';
'CAT_pfc_M13_[1]_3-54_0_0_1_0_.1_5_1';

'ITM_vvs_M11_[1]_3-54_0_0_1_0_.1_5_1';
'ITM_pfc_M11_[1]_3-54_0_0_1_0_.1_5_1';
'ITM_vvs_M12_[1]_3-54_0_0_1_0_.1_5_1';
'ITM_pfc_M12_[1]_3-54_0_0_1_0_.1_5_1';
'ITM_vvs_M13_[1]_3-54_0_0_1_0_.1_5_1';
'ITM_pfc_M13_[1]_3-54_0_0_1_0_.1_5_1';

'Alex_pfc_M11_[1-8]_3-54_1_0_1_0_.1_5_1';
'Alex_pfc_M12_[1-8]_3-54_1_0_1_0_.1_5_1';
'Alex_pfc_M13_[1-8]_3-54_1_0_1_0_.1_5_1';
'Alex_vvs_M11_[1-8]_3-54_1_0_1_0_.1_5_1';
'Alex_vvs_M12_[1-8]_3-54_1_0_1_0_.1_5_1';
'Alex_vvs_M13_[1-8]_3-54_1_0_1_0_.1_5_1';

'BLNETi_pfc_M11_[8-8-56]_3-54_1_0_1_0_.1_5_1';
'BLNETi_pfc_M12_[8-8-56]_3-54_1_0_1_0_.1_5_1';
'BLNETi_pfc_M13_[8-8-56]_3-54_1_0_1_0_.1_5_1';
'BLNETi_vvs_M11_[8-8-56]_3-54_1_0_1_0_.1_5_1';
'BLNETi_vvs_M12_[8-8-56]_3-54_1_0_1_0_.1_5_1';
'BLNETi_vvs_M13_[8-8-56]_3-54_1_0_1_0_.1_5_1';

'CORrt_pfc_M11_[2-2-8]_3-54_1_0_1_0_.1_5_1';
'CORrt_pfc_M12_[2-2-8]_3-54_1_0_1_0_.1_5_1';
'CORrt_pfc_M13_[2-2-8]_3-54_1_0_1_0_.1_5_1';
'CORrt_vvs_M11_[2-2-8]_3-54_1_0_1_0_.1_5_1';
'CORrt_vvs_M12_[2-2-8]_3-54_1_0_1_0_.1_5_1';
'CORrt_vvs_M13_[2-2-8]_3-54_1_0_1_0_.1_5_1';

'Alex_pfc_M11_[1-8]_3-54_0_0_1_0_.1_5_1';
'Alex_pfc_M12_[1-8]_3-54_0_0_1_0_.1_5_1';
'Alex_pfc_M13_[1-8]_3-54_0_0_1_0_.1_5_1';
'Alex_vvs_M11_[1-8]_3-54_0_0_1_0_.1_5_1';
'Alex_vvs_M12_[1-8]_3-54_0_0_1_0_.1_5_1';
'Alex_vvs_M13_[1-8]_3-54_0_0_1_0_.1_5_1';

'BLNETi_pfc_M11_[8-8-56]_3-54_0_0_1_0_.1_5_1';
'BLNETi_pfc_M12_[8-8-56]_3-54_0_0_1_0_.1_5_1';
'BLNETi_pfc_M13_[8-8-56]_3-54_0_0_1_0_.1_5_1';
'BLNETi_vvs_M11_[8-8-56]_3-54_0_0_1_0_.1_5_1';
'BLNETi_vvs_M12_[8-8-56]_3-54_0_0_1_0_.1_5_1';
'BLNETi_vvs_M13_[8-8-56]_3-54_0_0_1_0_.1_5_1';

'CORrt_pfc_M11_[2-2-8]_3-54_0_0_1_0_.1_5_1';
'CORrt_pfc_M12_[2-2-8]_3-54_0_0_1_0_.1_5_1';
'CORrt_pfc_M13_[2-2-8]_3-54_0_0_1_0_.1_5_1';
'CORrt_vvs_M11_[2-2-8]_3-54_0_0_1_0_.1_5_1';
'CORrt_vvs_M12_[2-2-8]_3-54_0_0_1_0_.1_5_1';
'CORrt_vvs_M13_[2-2-8]_3-54_0_0_1_0_.1_5_1';


'Alex_vvs_E123CAT1_[1-8]_3-54_0_0_1_0_.1_5_1';
'Alex_vvs_E123CAT2_[1-8]_3-54_0_0_1_0_.1_5_1';
'Alex_vvs_E123CAT3_[1-8]_3-54_0_0_1_0_.1_5_1';
'Alex_vvs_E123CAT4_[1-8]_3-54_0_0_1_0_.1_5_1';
'Alex_vvs_E123CAT5_[1-8]_3-54_0_0_1_0_.1_5_1';
'Alex_vvs_E123CAT6_[1-8]_3-54_0_0_1_0_.1_5_1';

'Alex_pfc_E123CAT1_[1-8]_3-54_0_0_1_0_.1_5_1';
'Alex_pfc_E123CAT2_[1-8]_3-54_0_0_1_0_.1_5_1';
'Alex_pfc_E123CAT3_[1-8]_3-54_0_0_1_0_.1_5_1';
'Alex_pfc_E123CAT4_[1-8]_3-54_0_0_1_0_.1_5_1';
'Alex_pfc_E123CAT5_[1-8]_3-54_0_0_1_0_.1_5_1';
'Alex_pfc_E123CAT6_[1-8]_3-54_0_0_1_0_.1_5_1';
% 
'Alex_vvs_M123CAT1_[1-8]_3-54_0_0_1_0_.1_5_1';
'Alex_vvs_M123CAT2_[1-8]_3-54_0_0_1_0_.1_5_1';
'Alex_vvs_M123CAT3_[1-8]_3-54_0_0_1_0_.1_5_1';
'Alex_vvs_M123CAT4_[1-8]_3-54_0_0_1_0_.1_5_1';
'Alex_vvs_M123CAT5_[1-8]_3-54_0_0_1_0_.1_5_1';
'Alex_vvs_M123CAT6_[1-8]_3-54_0_0_1_0_.1_5_1';
% 
'Alex_pfc_M123CAT1_[1-8]_3-54_0_0_1_0_.1_5_1';
'Alex_pfc_M123CAT2_[1-8]_3-54_0_0_1_0_.1_5_1';
'Alex_pfc_M123CAT3_[1-8]_3-54_0_0_1_0_.1_5_1';
'Alex_pfc_M123CAT4_[1-8]_3-54_0_0_1_0_.1_5_1';
'Alex_pfc_M123CAT5_[1-8]_3-54_0_0_1_0_.1_5_1';
'Alex_pfc_M123CAT6_[1-8]_3-54_0_0_1_0_.1_5_1';


'BLNETi_vvs_E123CAT1_[8-8-56]_3-54_0_0_1_0_.1_5_1';
'BLNETi_vvs_E123CAT2_[8-8-56]_3-54_0_0_1_0_.1_5_1';
'BLNETi_vvs_E123CAT3_[8-8-56]_3-54_0_0_1_0_.1_5_1';
'BLNETi_vvs_E123CAT4_[8-8-56]_3-54_0_0_1_0_.1_5_1';
'BLNETi_vvs_E123CAT5_[8-8-56]_3-54_0_0_1_0_.1_5_1';
'BLNETi_vvs_E123CAT6_[8-8-56]_3-54_0_0_1_0_.1_5_1';

'BLNETi_pfc_E123CAT1_[8-8-56]_3-54_0_0_1_0_.1_5_1';
'BLNETi_pfc_E123CAT2_[8-8-56]_3-54_0_0_1_0_.1_5_1';
'BLNETi_pfc_E123CAT3_[8-8-56]_3-54_0_0_1_0_.1_5_1';
'BLNETi_pfc_E123CAT4_[8-8-56]_3-54_0_0_1_0_.1_5_1';
'BLNETi_pfc_E123CAT5_[8-8-56]_3-54_0_0_1_0_.1_5_1';
'BLNETi_pfc_E123CAT6_[8-8-56]_3-54_0_0_1_0_.1_5_1';
% 
'BLNETi_vvs_M123CAT1_[8-8-56]_3-54_0_0_1_0_.1_5_1';
'BLNETi_vvs_M123CAT2_[8-8-56]_3-54_0_0_1_0_.1_5_1';
'BLNETi_vvs_M123CAT3_[8-8-56]_3-54_0_0_1_0_.1_5_1';
'BLNETi_vvs_M123CAT4_[8-8-56]_3-54_0_0_1_0_.1_5_1';
'BLNETi_vvs_M123CAT5_[8-8-56]_3-54_0_0_1_0_.1_5_1';
'BLNETi_vvs_M123CAT6_[8-8-56]_3-54_0_0_1_0_.1_5_1';
% 
'BLNETi_pfc_M123CAT1_[8-8-56]_3-54_0_0_1_0_.1_5_1';
'BLNETi_pfc_M123CAT2_[8-8-56]_3-54_0_0_1_0_.1_5_1';
'BLNETi_pfc_M123CAT3_[8-8-56]_3-54_0_0_1_0_.1_5_1';
'BLNETi_pfc_M123CAT4_[8-8-56]_3-54_0_0_1_0_.1_5_1';
'BLNETi_pfc_M123CAT5_[8-8-56]_3-54_0_0_1_0_.1_5_1';
'BLNETi_pfc_M123CAT6_[8-8-56]_3-54_0_0_1_0_.1_5_1';


'CORrt_vvs_E123CAT1_[2-2-8]_3-54_0_0_1_0_.1_5_1';
'CORrt_vvs_E123CAT2_[2-2-8]_3-54_0_0_1_0_.1_5_1';
'CORrt_vvs_E123CAT3_[2-2-8]_3-54_0_0_1_0_.1_5_1';
'CORrt_vvs_E123CAT4_[2-2-8]_3-54_0_0_1_0_.1_5_1';
'CORrt_vvs_E123CAT5_[2-2-8]_3-54_0_0_1_0_.1_5_1';
'CORrt_vvs_E123CAT6_[2-2-8]_3-54_0_0_1_0_.1_5_1';

'CORrt_pfc_E123CAT1_[2-2-8]_3-54_0_0_1_0_.1_5_1';
'CORrt_pfc_E123CAT2_[2-2-8]_3-54_0_0_1_0_.1_5_1';
'CORrt_pfc_E123CAT3_[2-2-8]_3-54_0_0_1_0_.1_5_1';
'CORrt_pfc_E123CAT4_[2-2-8]_3-54_0_0_1_0_.1_5_1';
'CORrt_pfc_E123CAT5_[2-2-8]_3-54_0_0_1_0_.1_5_1';
'CORrt_pfc_E123CAT6_[2-2-8]_3-54_0_0_1_0_.1_5_1';
% 
'CORrt_vvs_M123CAT1_[2-2-8]_3-54_0_0_1_0_.1_5_1';
'CORrt_vvs_M123CAT2_[2-2-8]_3-54_0_0_1_0_.1_5_1';
'CORrt_vvs_M123CAT3_[2-2-8]_3-54_0_0_1_0_.1_5_1';
'CORrt_vvs_M123CAT4_[2-2-8]_3-54_0_0_1_0_.1_5_1';
'CORrt_vvs_M123CAT5_[2-2-8]_3-54_0_0_1_0_.1_5_1';
'CORrt_vvs_M123CAT6_[2-2-8]_3-54_0_0_1_0_.1_5_1';
% 
'CORrt_pfc_M123CAT1_[2-2-8]_3-54_0_0_1_0_.1_5_1';
'CORrt_pfc_M123CAT2_[2-2-8]_3-54_0_0_1_0_.1_5_1';
'CORrt_pfc_M123CAT3_[2-2-8]_3-54_0_0_1_0_.1_5_1';
'CORrt_pfc_M123CAT4_[2-2-8]_3-54_0_0_1_0_.1_5_1';
'CORrt_pfc_M123CAT5_[2-2-8]_3-54_0_0_1_0_.1_5_1';
'CORrt_pfc_M123CAT6_[2-2-8]_3-54_0_0_1_0_.1_5_1';

                    
             };
    

for listi = 6 % 1:length(listF2sav)
    
    clearvars -except listF2sav listi nPerm 
        
    f2sav       = listF2sav{listi}
        
    cfg = getParams(f2sav);
    cfg.DNN_analysis = 1; 
    paths = load_paths_WM(cfg.brainROI, cfg.net2load);
    filelistSess = getFiles(paths.traces);
    
    t1 = datetime; 
    
    for sessi= 1:length(filelistSess) 
        disp(['File > ' num2str(sessi)]);
        load([paths.traces filelistSess{sessi}]);   

        nChans = size(cfg_contrasts.chanNames, 1); 
    
        if nChans > 1

            cfg_contrasts = getIdsWM(cfg.period, cfg_contrasts);

            if length(cfg_contrasts.oneListIds) > 1
                cfg_contrasts.oneListPow    = extract_power_WM (cfg_contrasts, cfg); % 
                
                cfg_contrasts = normalize_WM(cfg_contrasts, 1, 'sess', []);
                cfg_contrasts               = average_repetitions(cfg, cfg_contrasts);
                
                neuralRDMs                  = createNeuralRDMs(cfg, cfg_contrasts);
                networkRDMs                 = createNetworkRDMs(cfg, cfg_contrasts, sessi, paths);
        
        
                % % % restrict time for permutation data
                if strcmp(cfg.period(1), 'M')
                    if ndims(neuralRDMs) == 4
                        neuralRDMs = neuralRDMs(:,:,:,6:40); %frequency-resolved
                    else
                        neuralRDMs = neuralRDMs(:,:,6:40); %band analysis
                    end
                else
                    if ndims(neuralRDMs) == 4
                        neuralRDMs = neuralRDMs(:,:,:,6:15); %frequency-resolved
                    else
                        neuralRDMs = neuralRDMs(:,:,6:15); %band analysis
                    end
                end
                
                for permi = 1:nPerm
                    sC = size(networkRDMs, 2);
                    ids = randperm(sC);
                    networkRDMs = networkRDMs(:, ids, ids); 
                    nnFitPerm(permi, sessi,:,:, :)              = fitModel_WM(neuralRDMs, networkRDMs, cfg.fitMode); 
                end
            end
        end

    end
    
    mkdir ([paths.results.DNNs]);
    save([paths.results.DNNs f2sav '_' num2str(nPerm) 'p.mat'], 'nnFitPerm', '-v7.3');
    

end
   
t2 = datetime; 
etime(datevec(t2), datevec(t1))


%%  plot all layers FREQUENCY RESOLVED
%Network_ROI_ER_layers_freqs_avRepet_avTFV_fRes(0-1)_fitMode(0:noTrials; 1:Trials)_timeRes_win_mf
clear , clc

%f2sav = 'CORrt_pfc_E123_[2-2-8]_3-54_1_0_1_0_.1_5_1';
%f2sav = 'BLNETe_pfc_M123_[8-8-56]_3-54_1_0_1_0_.1_5_1';
%f2sav = 'CORrt_vvs_E123_[2-2-8]_3-54_1_0_1_0_.1_5_1';
%f2sav = 'Res18-8_vvs_MALL_[3]_3-54_0_0_1_0_.1_5_1'; 
%f2sav = 'BLNETi_pfc_M123_[8-8-56]_3-54_1_0_1_0_.1_5_1';
%f2sav = 'AlexEco_pfc_M123_[1-8]_3-54_1_0_1_0_.1_5_1'; 
%f2sav = 'CAT_pfc_E123_[1]_3-54_1_0_1_0_.1_5_1';
%f2sav = 'Res18-8_vvs_MALL_[3]_3-54_0_0_1_0_.1_5_1'; 
%f2sav = 'Res34-2_vvs_MALL_[1]_3-54_0_0_1_0_.1_5_1'; 
%f2sav =  'AE-t10_hipp_M123_[1-18]_3-54_1_0_1_0_.1_5_1'; 
%f2sav = 'BLNETi_pfc_E123_[8-8-56]_3-54_0_0_1_0_.1_5_1';
%f2sav = 'BNETi_pfc_M123_[1-7]_3-54_1_0_1_0_.1_5_1'; 
%f2sav = 'BDNETi_pfc_M123_[1-13]_3-54_1_0_1_0_.1_5_1'; 
%f2sav = 'BFNETi_pfc_M123_[1-7]_3-54_1_0_1_0_.1_5_1'; 
%f2sav = 'BKNETi_pfc_M123_[1-7]_3-54_1_0_1_0_.1_5_1'; 
%f2sav = 'CORz_pfc_M123_[1-4]_3-54_1_0_1_0_.1_5_1'; 
%f2sav = 'CORrt_pfc_M123II_[2-2-8]_3-54_0_0_1_0_.1_5_1';
%f2sav = 'ITM_vvs_M123_[1]_3-54_0_0_1_0_.1_5_1'; 
%f2sav = 'Alex_vvs_M12_[1-8]_3-54_1_0_1_0_.1_5_1';
%f2sav = 'ITM_vvs_E123CAT1_[1]_3-54_0_0_1_0_.1_5_1';
%f2sav =  'BLNETi_pfc_M123_[8-8-56]_3-54_1_0_1_0_.1_5_1_PC';

f2sav = 'BLNETi_vvs_E123CI_[8-8-56]_3-54_0_0_1_0_.1_5_1';



cfg = getParams(f2sav);
if strcmp(cfg.brainROI, 'vvs')
    sub2exc = [18 22];
elseif strcmp(cfg.brainROI, 'pfc')
    sub2exc = [1];
elseif strcmp(cfg.brainROI, 'hipp')
    sub2exc = [2];
end

paths = load_paths_WM(cfg.brainROI, cfg.net2load);
load([paths.results.DNNs f2sav '.mat']);


tiledlayout(8,9, 'TileSpacing', 'tight', 'Padding', 'compact');
if strcmp(cfg.period(1), 'M')
    set(gcf, 'Position', [100 100 1200 1000])
else
    set(gcf, 'Position', [100 100 700 1200])
    %set(gcf, 'Position', [100 100 1800 1000])
end

for layi = 1:size(nnFit{1}, 1)
    ax1 = nexttile;
    clear nnH
    for subji = 1:length(nnFit)
       if strcmp(cfg.period(1), 'M')
         nnH(subji, : ,:) = atanh(nnFit{subji, 1}(layi,:,1:40));
       else
         nnH(subji, : ,:) = atanh(nnFit{subji, 1}(layi,:,1:15));
         %nnH(subji, : ,:) = atanh(nnFit{subji, 1}(layi,:,1:134));
       end
    end
    
    nnH(sub2exc, :, :) = []; 
    nnH = squeeze(nnH);
    %[h p ci ts] = ttest(nnH, 0, "Tail","right");
    [h p ci ts] = ttest(nnH);
    h = squeeze(h); t = squeeze(ts.tstat); 
    %h(:, 1:5) = 0; % only sum p-values in clusters after the baseline
    
    d2p = squeeze(mean(nnH, 'omitnan'));
    freqs = 1:52; 
    clustinfo = bwconncomp(h);
    allClustInfo{layi} = clustinfo; 

    % store allTOBS
    if ~isempty(clustinfo.PixelIdxList)
        for pixi = 1:length(clustinfo.PixelIdxList)
             %if length(clustinfo.PixelIdxList{pixi}) > 1
                allTObs(layi, pixi, :) = sum(t(clustinfo.PixelIdxList{pixi}));
             %end        
        end
    else
        allTObs(layi, :, :) = 0;
    end


    
    if strcmp(cfg.period(1), 'M')
        times = 1:size(t, 2); 
    else
        times = 1:15; 
        %times = 1:134;
    end
    myCmap = colormap(brewermap([],'YlOrRd'));
    colormap(myCmap)
    contourf(times, freqs, t, 100, 'linecolor', 'none'); hold on; %colorbar
    contour(times, freqs, h, 1, 'Color', [0, 0, 0], 'LineWidth', 2);
    
    if strcmp(cfg.period(1), 'M')
        set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
        set(gca, 'xlim', [1 40], 'clim', [-5 5], 'FontSize', 10);
        plot([5 5],get(gca,'ylim'), 'k:','lineWidth', 2);
    else
        set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
        set(gca, 'FontSize', 8, 'clim', [-5 5]);
        plot([5 5],get(gca,'ylim'), 'k:','lineWidth', 2);
        
    end
    

end

exportgraphics(gcf, [paths.results.DNNs 'myP.png'], 'Resolution', 300); 



%%  Save BL-NET eco activations only in clusters DURING MAINTENANCE

f2sav = 'BLNETe_pfc_M123_[8-8-56]_3-54_1_0_1_0_.1_5_1';
%f2sav = 'AlexEco_pfc_M123_[1-8]_3-54_1_0_1_0_.1_5_1';

cfg = getParams(f2sav);
if strcmp(cfg.brainROI, 'vvs')
    sub2exc = [18 22];
elseif strcmp(cfg.brainROI, 'pfc')
    sub2exc = [1];
end

paths = load_paths_WM(cfg.brainROI, cfg.net2load);
load([paths.results.DNNs f2sav '.mat']);

clear nnHClust_pfc7 nnHClust_vvs4 nnHClust_vvs5 nnHClust_vvs6
for layi = 1:size(nnFit{1}, 1)
    clear nnH
    for subji = 1:length(nnFit)
           if strcmp(cfg.period(1), 'M')
             nnH(subji, : ,:) = atanh(nnFit{subji, 1}(layi,:,1:40));
           else
             nnH(subji, : ,:) = atanh(nnFit{subji, 1}(layi,:,1:15));
           end
    end
    
    nnH(sub2exc, :, :) = []; 


    for subji = 1:size(nnH, 1)
        nnHSubj = squeeze(nnH(subji, :, :)); 
        if strcmp(cfg.period(1), 'M')
            %if layi == 4 & strcmp(cfg.brainROI, 'vvs')
            %    nnHClust_vvs4(subji, :) = mean(nnHSubj(allClustInfo{4}.PixelIdxList{14}));
            %end 
            %if layi == 5 & strcmp(cfg.brainROI, 'vvs')
            %    nnHClust_vvs5(subji, :) = mean(nnHSubj(allClustInfo{5}.PixelIdxList{25})); 
            %end 
            %if layi == 6 & strcmp(cfg.brainROI, 'vvs')
            %    nnHClust_vvs6(subji, :) = mean(nnHSubj(allClustInfo{6}.PixelIdxList{17})); 
            %end
            if layi == 7 & strcmp(cfg.brainROI, 'pfc')
                nnHClust_pfc7(subji, :) = mean(nnHSubj(allClustInfo{7}.PixelIdxList{2})); 
            end
            %if layi == 4 & strcmp(cfg.brainROI, 'vvs') & strcmp(cfg.net2load, 'AlexEco')
            %    nnHClust_vvs4(subji, :) = mean(nnHSubj(allClustInfo{4}.PixelIdxList{14}));
            %end 
            %if layi == 5 & strcmp(cfg.brainROI, 'vvs') & strcmp(cfg.net2load, 'AlexEco')
                %nnHClust_vvs5(subji, :) = mean(nnHSubj(allClustInfo{5}.PixelIdxList{25})); 
            %end 
            %if layi == 6 & strcmp(cfg.brainROI, 'vvs') & strcmp(cfg.net2load, 'AlexEco')
            %    nnHClust_vvs6(subji, :) = mean(nnHSubj(allClustInfo{6}.PixelIdxList{17})); 
            %end
            if layi == 8 & strcmp(cfg.brainROI, 'pfc') & strcmp(cfg.net2load, 'AlexEco')
                nnHClust_pfc7(subji, :) = mean(nnHSubj(allClustInfo{7}.PixelIdxList{2})); 
            end
        elseif strcmp(cfg.period(1), 'E')
            nnHClust_vvsE1(subji, :) = mean(nnHSubj(allClustInfo{1}.PixelIdxList{1})); 
            nnHClust_vvsE2(subji, :) = mean(nnHSubj(allClustInfo{2}.PixelIdxList{1})); 
            nnHClust_vvsE3(subji, :) = mean(nnHSubj(allClustInfo{3}.PixelIdxList{1})); 
            nnHClust_vvsE4(subji, :) = mean(nnHSubj(allClustInfo{4}.PixelIdxList{1})); 
            nnHClust_vvsE5(subji, :) = mean(nnHSubj(allClustInfo{5}.PixelIdxList{1})); 
            nnHClust_vvsE6(subji, :) = mean(nnHSubj(allClustInfo{6}.PixelIdxList{1})); 
            nnHClust_vvsE7(subji, :) = mean(nnHSubj(allClustInfo{7}.PixelIdxList{1})); 
        end
    end

    
   
end

%% 
%[h p ci ts] = ttest(nnHClust_vvs6);
[h p ci ts] = ttest(nnHClust_pfc7);
%[h p ci ts] = ttest(nnHClust_vvsE7);
h = squeeze(h); t = squeeze(ts.tstat); 

disp (['t = ' num2str(ts.tstat) '  ' ' p = ' num2str(p)]);

%% plot 7 bar
clear data
data.data = [nnHClust_vvsE1 nnHClust_vvsE2 nnHClust_vvsE3 nnHClust_vvsE4 ...
             nnHClust_vvsE5 nnHClust_vvsE6 nnHClust_vvsE7]; 

figure(2); set(gcf,'Position', [0 0 500 620]); 
mean_S = mean(data.data, 1);
hb = scatter([1 2 3 4 5 6 7], data.data, 'k'); hold on;
%set(hb, 'lineWidth', 0.01, 'Marker', '.', 'MarkerSize',45);hold on;

h = bar (mean_S);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 3);
set(gca,'XTick',[1 2 3 4 5 6 7],'XTickLabel',{'', ''}, 'FontSize', 30, 'linew',2, 'xlim', [0 8] );
set(gca, 'ylim', [-.1 .3])
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 3);

[h p ci t] = ttest (data.data(:,1), data.data(:,2));
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth', 3);


exportgraphics(gcf, 'allM.png', 'Resolution', 300);

%% plot 3 bar
clear data
data.data = [nnHClust_vvs4 nnHClust_vvs5 nnHClust_vvs6]; 

figure(2); set(gcf,'Position', [0 0 500 620]); 
mean_S = mean(data.data, 1);
hb = scatter([1 2 3], data.data, 'k'); hold on;
%set(hb, 'lineWidth', 0.01, 'Marker', '.', 'MarkerSize',45);hold on;

h = bar (mean_S);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 3);
set(gca,'XTick',[1 2 3],'XTickLabel',{'', ''}, 'FontSize', 30, 'linew',2, 'xlim', [0 4] );
set(gca, 'ylim', [-.05 .1])
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 3);

[h p ci t] = ttest (data.data(:,1), data.data(:,2));
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth', 3);


exportgraphics(gcf, 'allM.png', 'Resolution', 300);

%% plot 1 bar
clear data
data.data = [nnHClust_pfc7]; 

figure(2); set(gcf,'Position', [0 0 500 620]); 
mean_S = mean(data.data, 1);
hb = scatter([1], data.data, 100, 'k'); hold on;
%set(hb, 'lineWidth', 0.01, 'Marker', '.', 'MarkerSize',45);hold on;

h = bar (mean_S);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 3);
set(gca,'XTick',[1],'XTickLabel',{'', ''}, 'FontSize', 30, 'linew',2, 'xlim', [0 2] );
set(gca, 'ylim', [-.05 .1])
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 3);

[h p ci t] = ttest (data.data(:,1));
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth', 3);


exportgraphics(gcf, 'allM.png', 'Resolution', 300);

%%  plot all layers FREQUENCY RESOLVED FANCY PLOT
%Network_ROI_ER_layers_freqs_avRepet_avTFV_fRes(0-1)_fitMode(0:noTrials; 1:Trials)_timeRes_win_mf
clear , clc
%f2sav = 'CORrt_vvs_E123_[2-2-8]_3-54_1_0_1_0_.1_5_1';
%f2sav = 'Res18-8_vvs_MALL_[3]_3-54_0_0_1_0_.1_5_1'; 
%f2sav = 'AlexEco_pfc_M123_[1-8]_3-54_1_0_1_0_.1_5_1'; 

%f2sav = 'BLNETe_pfc_M123_[8-8-56]_3-54_1_0_1_0_.1_5_1';
f2sav = 'BLNETi_pfc_M123_[8-8-56]_3-54_1_0_1_0_.1_5_1';
%f2sav = 'CAT_vvs_M123_[1]_3-54_1_0_1_0_.1_5_1';
%f2sav = 'Alex_pfc_E123_[1-8]_3-54_1_0_1_0_.1_5_1';
%f2sav = 'ITM_vvs_E123_[1]_3-54_0_0_1_0_.1_5_1'; 

cfg = getParams(f2sav);
if strcmp(cfg.brainROI, 'vvs')
    sub2exc = [18 22];
elseif strcmp(cfg.brainROI, 'pfc')
    sub2exc = [1];
end

paths = load_paths_WM(cfg.brainROI, cfg.net2load);
load([paths.results.DNNs f2sav '.mat']);


tiledlayout(8,8, 'TileSpacing', 'compact', 'Padding', 'compact');
if strcmp(cfg.period(1), 'M')
    set(gcf, 'Position', [100 100 1800 1300])
else
    set(gcf, 'Position', [100 100 670 1300])
end

for layi = 1:size(nnFit{1}, 1)
    ax1 = nexttile;
    clear nnH
    for subji = 1:length(nnFit)
       if strcmp(cfg.period(1), 'M')
         nnH(subji, : ,:) = atanh(nnFit{subji, 1}(layi,:,1:40));
       else
         nnH(subji, : ,:) = atanh(nnFit{subji, 1}(layi,:,1:15));
       end
    end
    
    nnH(sub2exc, :, :) = []; 
    nnH = squeeze(nnH);
    %[h p ci ts] = ttest(nnH, 0, "Tail","right");
    [h p ci ts] = ttest(nnH);
    h = squeeze(h); t = squeeze(ts.tstat); 
    h(:, 1:5) = 0; % only sum p-values in clusters after the baseline

    d2p = squeeze(mean(nnH, 'omitnan'));
    freqs = 1:520; 
    clustinfo = bwconncomp(h);
    
    % % % % % % 
    if strcmp(cfg.period(1), 'M')
        h = zeros(52, 40); 
        %h(clustinfo.PixelIdxList{1}) = 1;
        
        if strcmp(cfg.net2load, 'BLNETi') & strcmp(cfg.brainROI, 'vvs') & layi == 4
            h(clustinfo.PixelIdxList{14}) = 1;
        end
        if strcmp(cfg.net2load, 'BLNETi') & strcmp(cfg.brainROI, 'vvs') & layi == 5
            h(clustinfo.PixelIdxList{25}) = 1;
        end
        if strcmp(cfg.net2load, 'BLNETi') & strcmp(cfg.brainROI, 'vvs') & layi == 6
            h(clustinfo.PixelIdxList{17}) = 1;
        end
        if strcmp(cfg.net2load, 'BLNETi') & strcmp(cfg.brainROI, 'pfc') & layi == 7
            h(clustinfo.PixelIdxList{2}) = 1;
        end
        
        if strcmp(cfg.net2load, 'CORrt') & strcmp(cfg.brainROI, 'vvs') & layi == 4
            h(clustinfo.PixelIdxList{15}) = 1;
        end
        if strcmp(cfg.net2load, 'CORrt') & strcmp(cfg.brainROI, 'pfc') & layi == 1
            h(clustinfo.PixelIdxList{1}) = 1;
        end
        if strcmp(cfg.net2load, 'CORrt') & strcmp(cfg.brainROI, 'pfc') & layi == 4
            h(clustinfo.PixelIdxList{3}) = 1;
        end

    else
        h = zeros(52, 15); 
        if strcmp (cfg.net2load, 'Alex') & strcmp(cfg.brainROI, 'vvs') 
            h(clustinfo.PixelIdxList{1}) = 1;
            if layi == 4
                h(clustinfo.PixelIdxList{2}) = 1;
            end
            if layi == 6
                h(clustinfo.PixelIdxList{2}) = 1;
            end
        end
        if strcmp (cfg.net2load, 'CORrt') & strcmp(cfg.brainROI, 'vvs') & layi == 1
            h(clustinfo.PixelIdxList{1}) = 1;
            h(clustinfo.PixelIdxList{2}) = 1;
            h(clustinfo.PixelIdxList{3}) = 1;
            h(clustinfo.PixelIdxList{4}) = 1;
        end
        if strcmp (cfg.net2load, 'CORrt') & strcmp(cfg.brainROI, 'vvs') & layi == 2
            h(clustinfo.PixelIdxList{1}) = 1;
        end
        if strcmp (cfg.net2load, 'CORrt') & strcmp(cfg.brainROI, 'vvs') & layi == 3
            h(clustinfo.PixelIdxList{1}) = 1;
        end
        if strcmp (cfg.net2load, 'CORrt') & strcmp(cfg.brainROI, 'vvs') & layi == 4
            h(clustinfo.PixelIdxList{1}) = 1;
            h(clustinfo.PixelIdxList{2}) = 1;
        end
        if strcmp (cfg.net2load, 'CAT') & strcmp(cfg.brainROI, 'vvs')
            h(clustinfo.PixelIdxList{1}) = 1;
        end
        if strcmp (cfg.net2load, 'CAT') & strcmp(cfg.brainROI, 'pfc')
            h(clustinfo.PixelIdxList{2}) = 1;
            h(clustinfo.PixelIdxList{3}) = 1;
        end
        if strcmp(cfg.net2load, 'BLNETe') & strcmp(cfg.brainROI, 'vvs') 
            h(clustinfo.PixelIdxList{1}) = 1;
            if layi ==4
                h(clustinfo.PixelIdxList{2}) = 1;
            end
            if layi == 7 
                h(clustinfo.PixelIdxList{2}) = 1;
            end
        end

    
    
    end

    % store allTOBS
    if ~isempty(clustinfo.PixelIdxList)
        for pixi = 1:length(clustinfo.PixelIdxList)
             %if length(clustinfo.PixelIdxList{pixi}) > 1
                allTObs(layi, pixi, :) = sum(t(clustinfo.PixelIdxList{pixi}));
             %end        
        end
    else
        allTObs(layi, :, :) = 0;
    end


    
    if strcmp(cfg.period(1), 'M')
        times = 1:size(t, 2)*10; 
    else
        times = 1:150; 
    end
    myCmap = colormap(brewermap([],'*spectral'));
    colormap(myCmap)
    contourf(times, freqs,myresizem(t, 10), 100, 'linecolor', 'none'); hold on; %colorbar
    contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 2);
    
    if strcmp(cfg.period(1), 'M')
        set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
        set(gca, 'xlim', [1 400], 'clim', [-5 5], 'FontSize', 10);
        plot([45 45],get(gca,'ylim'), 'k:','lineWidth', 4);
    else
        set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
        set(gca, 'FontSize', 8, 'clim', [-5 5]);
        plot([45 45],get(gca,'ylim'), 'k:','lineWidth', 4);
        
    end
    

end

exportgraphics(gcf, [paths.results.DNNs 'myP.png'], 'Resolution', 300); 
close all; 

%% plot final figure only last layer / time point ENCODING
times = 1:150;
freqs = 1:520; 
%h = zeros(52, 15); 

%h(clustinfo.PixelIdxList{1}) = 1; % Alexnet VVS Layer 1-3, 5

%h(clustinfo.PixelIdxList{1}) = 1; % Alexnet VVS Layer 4 and 6
%h(clustinfo.PixelIdxList{2}) = 1; % Alexnet VVS Layer 4 and 6


%h(clustinfo.PixelIdxList{5}) = 1; % BLNETi PFC
%h(clustinfo.PixelIdxList{1}) = 1; % BLNETi VVS

%h(clustinfo.PixelIdxList{1}) = 1; % BLNETe cluster 1
%h(clustinfo.PixelIdxList{2}) = 1; % BLNETe cluster 2

%h(clustinfo.PixelIdxList{2}) = 1; %Cornet
%h(clustinfo.PixelIdxList{4}) = 1; %Cornet cluster 2


figure; set(gcf, 'Position', [100 100 200 400])
contourf(times, freqs, myresizem(t, 10), 100, 'linecolor', 'none'); hold on; %colorbar
contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 4);

set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
set(gca, 'xlim', [1 150], 'clim', [-5 5], 'FontSize', 10);
set(gca, 'clim', [-4 4], 'FontSize', 10);
plot([45 45],get(gca,'ylim'), 'k:','lineWidth', 5); 
%colorbar

myCmap = colormap(brewermap([],'*Spectral'));
colormap(myCmap)

%exportgraphics(gcf, [paths.results.DNNs 'myP.png'], 'Resolution', 300); 
exportgraphics(gcf, ['myP.png'], 'Resolution', 300); 

%% plot final figure only last layer MAINTENANCE
times = 1:400;
freqs = 1:520; 
%h = zeros(52, 39); 

h(clustinfo.PixelIdxList{2}) = 1; %BLNETi PFC

%h(clustinfo.PixelIdxList{3}) = 1; %CORNET PFC
%h(clustinfo.PixelIdxList{15}) = 1; %CORNET VVS

%h(clustinfo.PixelIdxList{5}) = 1; %category model

%h(clustinfo.PixelIdxList{2}) = 1; %pfc - Cornet
%h(clustinfo.PixelIdxList{23}) = 1; %pfc - Cornet

%h(clustinfo.PixelIdxList{8}) = 1; %pfc - ecoset
%h(clustinfo.PixelIdxList{10}) = 1; %vvs1 - ecoset
%h(clustinfo.PixelIdxList{23}) = 1; %vvs2 - ecoset


figure; set(gcf, 'Position', [1000 918 560 420])
myCmap = colormap(brewermap([],'*Spectral'));
colormap(myCmap)
%contourf(times, freqs, t, 100, 'linecolor', 'none'); hold on; %colorbar
%contour(times, freqs, h, 1, 'Color', [0, 0, 0], 'LineWidth', 2);
contourf(times, freqs, myresizem(t, 10), 100, 'linecolor', 'none'); hold on; %colorbar
contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 4);

set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
set(gca, 'xlim', [1 390], 'clim', [-4 4], 'FontSize', 10);
%set(gca, 'clim', [-4 4], 'FontSize', 10);
plot([45 45],get(gca,'ylim'), 'k:','lineWidth', 5); 
%colorbar

%exportgraphics(gcf, [paths.results.DNNs 'myP.png'], 'Resolution', 300); 
exportgraphics(gcf, ['myP.png'], 'Resolution', 300); 



%%  plot all layers FREQUENCY RESOLVED
%Network_ROI_ER_layers_freqs_avRepet_avTFV_fRes(0-1)_fitMode(0:noTrials; 1:Trials)_timeRes_win_mf
clear , clc

%f2sav = 'CORrt_pfc_E123_[2-2-8]_3-54_1_0_1_0_.1_5_1';
%f2sav = 'BLNETe_pfc_M123_[8-8-56]_3-54_1_0_1_0_.1_5_1';
%f2sav = 'CORrt_vvs_E123_[2-2-8]_3-54_1_0_1_0_.1_5_1';
%f2sav = 'Res18-8_vvs_MALL_[3]_3-54_0_0_1_0_.1_5_1'; 
%f2sav = 'BLNETi_pfc_M123_[8-8-56]_3-54_1_0_1_0_.1_5_1';
%f2sav = 'AlexEco_pfc_M123_[1-8]_3-54_1_0_1_0_.1_5_1'; 
%f2sav = 'CAT_pfc_E123_[1]_3-54_1_0_1_0_.1_5_1';
%f2sav = 'Res18-8_vvs_MALL_[3]_3-54_0_0_1_0_.1_5_1'; 
%f2sav = 'Res34-2_vvs_MALL_[1]_3-54_0_0_1_0_.1_5_1'; 
%f2sav =  'AE-t10_hipp_M123_[1-18]_3-54_1_0_1_0_.1_5_1'; 
%f2sav = 'Alex_pfc_M123I_[1-8]_3-54_1_0_1_0_.1_5_1';
%f2sav = 'BLNETi_pfc_E123_[8-8-56]_3-54_0_0_1_0_.1_5_1';
%f2sav = 'BNETi_pfc_M123_[1-7]_3-54_1_0_1_0_.1_5_1'; 
%f2sav = 'BDNETi_vvs_M123_[1-13]_3-54_1_0_1_0_.1_5_1'; 
%f2sav = 'BFNETi_pfc_M123_[1-7]_3-54_1_0_1_0_.1_5_1'; 
%f2sav = 'BKNETi_pfc_M123_[1-7]_3-54_1_0_1_0_.1_5_1'; 
%f2sav = 'CORz_vvs_M123_[1-4]_3-54_1_0_1_0_.1_5_1'; 
%f2sav = 'CORrt_pfc_E123_[2-2-8]_3-54_0_0_1_0_.1_5_1';

f2sav =  'ITM_pfc_M123_[1]_3-54_0_0_1_0_.1_5_1';



cfg = getParams(f2sav);
if strcmp(cfg.brainROI, 'vvs')
    sub2exc = [18 22];
elseif strcmp(cfg.brainROI, 'pfc')
    sub2exc = [1];
elseif strcmp(cfg.brainROI, 'hipp')
    sub2exc = [2];
end

paths = load_paths_WM(cfg.brainROI, cfg.net2load);
load([paths.results.DNNs f2sav '.mat']);


tiledlayout(8,9, 'TileSpacing', 'tight', 'Padding', 'compact');
if strcmp(cfg.period(1), 'M')
    set(gcf, 'Position', [100 100 1200 1000])
else
    set(gcf, 'Position', [100 100 700 1200])
end

for layi = 1:size(nnFit{1}, 1)
    ax1 = nexttile;
    clear nnH
    for subji = 1:length(nnFit)
       if strcmp(cfg.period(1), 'M')
         nnH(subji, : ,:) = atanh(nnFit{subji, 1}(layi,:,1:40));
       else
         nnH(subji, : ,:) = atanh(nnFit{subji, 1}(layi,:,1:15));
       end
    end
    
    nnH(sub2exc, :, :) = []; 
    nnH = squeeze(nnH);
    %[h p ci ts] = ttest(nnH, 0, "Tail","right");
    [h p ci ts] = ttest(nnH);
    h = squeeze(h); t = squeeze(ts.tstat); 
    h(:, 1:5) = 0; % only sum p-values in clusters after the baseline
    
    d2p = squeeze(mean(nnH, 'omitnan'));
    freqs = 1:52; 
    clustinfo = bwconncomp(h);
    allClustInfo{layi} = clustinfo; 

    % store allTOBS
    if ~isempty(clustinfo.PixelIdxList)
        for pixi = 1:length(clustinfo.PixelIdxList)
             %if length(clustinfo.PixelIdxList{pixi}) > 1
                allTObs(layi, pixi, :) = sum(t(clustinfo.PixelIdxList{pixi}));
             %end        
        end
    else
        allTObs(layi, :, :) = 0;
    end


    
    if strcmp(cfg.period(1), 'M')
        times = 1:size(t, 2); 
    else
        times = 1:15; 
    end
    myCmap = colormap(brewermap([],'YlOrRd'));
    colormap(myCmap)
    contourf(times, freqs, t, 100, 'linecolor', 'none'); hold on; %colorbar
    contour(times, freqs, h, 1, 'Color', [0, 0, 0], 'LineWidth', 2);
    
    if strcmp(cfg.period(1), 'M')
        set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
        set(gca, 'xlim', [1 40], 'clim', [-5 5], 'FontSize', 10);
        plot([5 5],get(gca,'ylim'), 'k:','lineWidth', 2);
    else
        set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
        set(gca, 'FontSize', 8, 'clim', [-5 5]);
        plot([5 5],get(gca,'ylim'), 'k:','lineWidth', 2);
        
    end
    

end

exportgraphics(gcf, [paths.results.DNNs 'myP.png'], 'Resolution', 300); 


%% COMPUTE CLUSTERS in each permutation FREQUENCY RESOLVED
clc
clearvars -except allTObs f2sav
%Network_ROI_ER_layers_freqs_avRepet_avTFV_fRes(0-1)_fitMode(0:noTrials; 1:Trials)__timeRes__win__mf

% use the same name as the previously plotted file
f2sav = [f2sav '_1000p.mat'];

cfg = getParams(f2sav);
paths = load_paths_WM(cfg.brainROI, cfg.net2load);

cd ([paths.results.DNNs])
load(f2sav);

if strcmp(cfg.brainROI, 'pfc')
    sub2exc = [1];
elseif strcmp(cfg.brainROI, 'vvs')
    sub2exc = [18 22];
end

f2t = strsplit(f2sav, '_');
nPerm = double(string((f2t{end}(1:end-5))));
nLays = size(nnFitPerm, 3);

for permi = 1:nPerm
    
    for layi = 1:nLays
        
        dataP = atanh(squeeze(nnFitPerm(permi, :,layi, :,:)));
        %dataP = squeeze(nnFitPerm(permi, :,layi, :,:));
        dataP(sub2exc, :, :) = []; 
        [h p ci ts] = ttest(dataP);
        h = squeeze(h); t = squeeze(ts.tstat);
        h(:, 1:5) = 0; % only sum t-values in cluster after the baseline
        
        
        clear allSTs  
        clustinfo = bwconncomp(h);
        for pxi = 1:length(clustinfo.PixelIdxList)
           allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
        end
        
        %sort all clusters 
        if exist('allSTs')
            [max2u id] = max(allSTs);
        else
            allSTs = 0; 
            id = 1;
        end
        
        max_clust_sum_perm(permi,layi,:) = allSTs(id); 
    end

end

cd (paths.github)


%% compute p value bands for all layers FREQ RES
clc
clear p

for layi = 1:size(allTObs, 1)
    clear mcsR mcsP
    for pixi = 1:size(allTObs, 2)
        mcsR = allTObs(layi, pixi); 

        if mcsR ~= 0
            mcsP = squeeze(max_clust_sum_perm(:,layi));
        
            %allAb = mcsP(abs(mcsP) > abs(mcsR));
            allAb = mcsP(mcsP > mcsR);
            p(layi, pixi,:) = 1 - ((nPerm-1) - (length (allAb)))  / nPerm;
        else
            p(layi, pixi,:) = nan; 
        end
    end
end

p (p==1.0010 | p == 1) = nan; 
p_ranked = p; p_ranked(isnan(p_ranked)) = []; 
p_ranked = sort(p_ranked(:));
p

%% p last layer only
p = p (end,:);
p_ranked = p; p_ranked(p_ranked == 0 | isnan(p_ranked)) = []; 
p_ranked = sort(p_ranked(:))





%% SAVE ONLY NEURAL RDMS 
clear, clc
%Network_ROI_EoM_layers_freqs_avRepet_avTimeFeatVect_freqResolv(0-1)__fitMode(0:noTrials; 1:Trials)__timeRes__win-width__mf
    
listF2sav = {

%'Alex_pfc_E123_[1-8]_3-54_1_0_1_0_.1_5_1'; 
%'Alex_vvs_E123_[1-8]_3-54_1_0_1_0_.1_5_1'; 
%'Alex_pfc_M123_[1-8]_3-54_1_0_1_0_.1_5_1'; 
%'Alex_vvs_M123_[1-8]_3-54_1_0_1_0_.1_5_1'; 

%'Alex_pfc_E123_[1-8]_3-54_0_0_1_0_.1_5_1'; 
%'Alex_vvs_E123_[1-8]_3-54_0_0_1_0_.1_5_1'; 
%'Alex_pfc_M123_[1-8]_3-54_0_0_1_0_.1_5_1'; 
'Alex_vvs_M123_[1-8]_3-54_0_0_1_0_.1_5_1'; 

};   


for listi = 1:length(listF2sav)
    disp(['File > ' num2str(listi) '      ' listF2sav{listi}]);
    clearvars -except listF2sav listi 
        
    f2sav       = listF2sav{listi}; 
    cfg = getParams(f2sav);
    cfg.DNN_analysis = 1; 
    paths = load_paths_WM(cfg.brainROI, cfg.net2load);
    filelistSess = getFiles(paths.traces);
    
    for sessi= 1:length(filelistSess) %this one starts at 1 and not at 3
        load([paths.traces filelistSess{sessi}]);   
       
        [cfg_contrasts] = getIdsWM(cfg.period, cfg_contrasts);
        cfg_contrasts.oneListPow    = extract_power_WM (cfg_contrasts, cfg); % 
        cfg_contrasts               = normalize_WM(cfg_contrasts, 1, 'sess', []);
        cfg_contrasts               = average_repetitions(cfg, cfg_contrasts);
        neuralRDMs                  = createNeuralRDMs(cfg, cfg_contrasts);
        allNeuralRDMS{sessi,1}      = neuralRDMs; 
        allNeuralRDMS{sessi,2}      = cfg_contrasts.oneListIds; 
        
    end
    
    save([paths.results.neuralRDMS.IT f2sav '.mat'], 'allNeuralRDMS');

end

%% Process neural RDMs and compute CCI for CATEGORIES

clear
paths = load_paths_WM('vvs', 'none');
filelistSess = getFilesWM(paths.results.neuralRDMS.CAT);



t1 = datetime; 

for sessi= 1:length(filelistSess) %this one starts at 1 and not at 3
    disp(['File > ' num2str(sessi)]);
    clearvars -except paths filelistSess sessi

    load([paths.results.neuralRDMS.CAT filelistSess{sessi}]);   
    nSubjs = size(allNeuralRDMS, 1); 
    nFreqs = size(allNeuralRDMS{1}, 3); 
    nTimes = size(allNeuralRDMS{1}, 4); 

    clear CCI
    for subji = 1:nSubjs
        neuralRDMs  = allNeuralRDMS{subji, 1}; 
        ids         = allNeuralRDMS{subji, 2};
        CM          = squeeze(load_CATMODEL_activ(ids));
        

        for freqi  = 1:nFreqs
            for timei = 1:nTimes
                neuralRDM = neuralRDMs(:, :, freqi, timei); 
                W = neuralRDM(CM==1); 
                B = neuralRDM(CM==0); 
                CCI(subji, freqi, timei,:) = mean(W) - mean(B); 
            end
        end
    end


    save([paths.results.neuralRDMS.CAT filelistSess{sessi}(1:end-4) '_CCI.mat'], 'CCI');

end

%% plot mean CCI 


clear
paths = load_paths_WM('vvs', 'none');
cd(paths.results.neuralRDMS.CATPlts);

filelistSess = dir('*mat'); filelistSess = [{filelistSess.name}'];
for sessi= 1:length(filelistSess) 
    
    disp(['File > ' num2str(sessi)]);

    clearvars -except paths filelistSess sessi
    load(filelistSess{sessi});   
    
    mCCI = squeeze(mean(CCI)); 
    
    [h p ci ts] = ttest(CCI); 
    h = squeeze(h); t = squeeze(ts.tstat); 
    
    if size(mCCI, 2) == 21
        times = 1:210;
        freqs = 1:520; 
        figure; set(gcf, 'Position', [100 100 200 400])
        contourf(times, freqs, myresizem(mCCI, 10), 100, 'linecolor', 'none'); hold on; colorbar
        %contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 4);
        
        set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
        set(gca, 'xlim', [1 150], 'FontSize', 10); 
        %set(gca, 'clim', [-.5 .5], 'FontSize', 10);
        plot([45 45],get(gca,'ylim'), 'k:','lineWidth', 5); 
        %colorbar
    else
        times = 1:460;
        freqs = 1:520; 
        figure; set(gcf, 'Position', [1000 918 560 420])
        myCmap = colormap(brewermap([],'*Spectral'));
        colormap(myCmap)
        contourf(times, freqs, myresizem(mCCI, 10), 100, 'linecolor', 'none'); hold on; colorbar
        %contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 4);
        
        set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
        set(gca, 'xlim', [1 390], 'FontSize', 10);
        %set(gca, 'clim', [-4 4], 'FontSize', 10);
        plot([45 45],get(gca,'ylim'), 'k:','lineWidth', 5); 
        %colorbar


    end
    myCmap = colormap(brewermap([],'*Spectral'));
    colormap(myCmap)
    
    exportgraphics(gcf, [num2str(sessi) '_myP.png'], 'Resolution', 300); 


end





%% COMPUTE PERMUTATIONS

clear

nPerm = 100; 
paths = load_paths_WM('vvs', 'none');
filelistSess = getFilesWM(paths.results.neuralRDMS.CAT);

for sessi= 4 %1:length(filelistSess) %this one starts at 1 and not at 3
    disp(['File > ' num2str(sessi)]);
    clearvars -except paths filelistSess sessi nPerm

    load([paths.results.neuralRDMS.CAT filelistSess{sessi}]);   
    nSubjs = size(allNeuralRDMS, 1); 
    nFreqs = size(allNeuralRDMS{1}, 3); 
    nTimes = size(allNeuralRDMS{1}, 4); 
    
    clear CCIPerm
    for subji = 1:nSubjs
        subji
        neuralRDMs  = allNeuralRDMS{subji, 1}; 
        ids         = allNeuralRDMS{subji, 2};

        for permi = 1:nPerm
            id4perm = randperm(length(ids)); 
            idsPERM = ids(id4perm);
            
            CM = squeeze(load_CATMODEL_activ(idsPERM));
            CM(CM==0) = 2; CM = tril(CM, -1); CM(CM==0) = 3; CM(CM==2) = 0; 
            for freqi  = 1:nFreqs
                for timei = 1:nTimes
                    neuralRDM = neuralRDMs(:, :, freqi, timei); 
                    W = neuralRDM(CM==1); 
                    B = neuralRDM(CM==0); 
                    CCIPerm(permi, subji, freqi, timei,:) = mean(W) - mean(B); 
                end
            end
        end
    end


    save([paths.results.neuralRDMS.CATPlts filelistSess{sessi}(1:end-4) '_CCI_' num2str(nPerm) 'p.mat'], 'CCIPerm');

end


%% plot with OUTLINE FROM PERM


clear
nPerm = 100;
paths = load_paths_WM('vvs', 'none');
cd(paths.results.neuralRDMS.CATPlts);

filelistSess = dir('*mat'); filelistSess = [{filelistSess.name}'];
for sessi= 6 %:length(filelistSess) %this one starts at 1 and not at 3
    
    disp(['File > ' num2str(sessi)]);

    clearvars -except paths filelistSess sessi nPerm
    load(filelistSess{sessi});   
    load(filelistSess{sessi+1});   
    
    mCCI = squeeze(mean(CCI)); 
    mCCIPERM = squeeze(mean(CCIPerm, 2)); 
    
    nFreq = size(mCCI,1); nTimes = size(mCCI, 2); 
    clear p
    for freqi = 1:nFreq
        for timei = 1:nTimes
            allCCIPerm = mCCIPERM(:, freqi, timei); 
            obsCCI = mCCI(freqi, timei); 
            allAB = allCCIPerm(allCCIPerm > obsCCI);
            p(freqi, timei) = 1 - ((nPerm-1) - (length (allAB)))  / nPerm;

        end
    end
    h = p<.05; 
    
    

    
    if size(mCCI, 2) == 21
        times = 1:210;
        freqs = 1:520; 
        figure; set(gcf, 'Position', [100 100 200 400])
        contourf(times, freqs, myresizem(mCCI, 10), 100, 'linecolor', 'none'); hold on; colorbar
        contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 4);
        
        set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
        set(gca, 'xlim', [1 150], 'FontSize', 10); 
        set(gca, 'clim', [-.015 .015], 'FontSize', 10);
        plot([45 45],get(gca,'ylim'), 'k:','lineWidth', 5); 
        %colorbar
    else
        times = 1:460;
        freqs = 1:520; 
        figure; set(gcf, 'Position', [1000 918 560 420])
        myCmap = colormap(brewermap([],'*Spectral'));
        colormap(myCmap)
        contourf(times, freqs, myresizem(mCCI, 10), 100, 'linecolor', 'none'); hold on; colorbar
        contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 4);
        
        set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
        set(gca, 'xlim', [1 390], 'FontSize', 10);
        set(gca, 'clim', [-.015 .015], 'FontSize', 10);
        plot([45 45],get(gca,'ylim'), 'k:','lineWidth', 5); 
        %colorbar


    end
    myCmap = colormap(brewermap([],'*Spectral'));
    colormap(myCmap)
    
    exportgraphics(gcf, [num2str(sessi) '_myP.png'], 'Resolution', 300); 


end





%% Process neural RDMs and compute ICI (ITEMS)

clear
paths = load_paths_WM('vvs', 'none');
filelistSess = getFilesWM(paths.results.neuralRDMS.IT);

t1 = datetime; 

for sessi= 3%1:length(filelistSess) %this one starts at 1 and not at 3
    disp(['File > ' num2str(sessi)]);
    load([paths.results.neuralRDMS.IT filelistSess{sessi}]);   
    nSubjs = size(allNeuralRDMS, 1); 
    nFreqs = size(allNeuralRDMS{1}, 3); 
    nTimes = size(allNeuralRDMS{1}, 4); 

    clear ICI
    for subji = 1:nSubjs
        neuralRDMs  = allNeuralRDMS{subji, 1}; 
        ids         = allNeuralRDMS{subji, 2};
        
        ITM = squeeze(load_ITMODEL_activ(ids));
        
        for freqi  = 1:nFreqs
            for timei = 1:nTimes
                neuralRDM = neuralRDMs(:, :, freqi, timei); 
                W = neuralRDM(ITM==1); 
                B = neuralRDM(ITM==0); 
                ICI(subji, freqi, timei,:) = mean(W, 'omitnan') - mean(B, 'omitnan'); 
            end
        end
    end


    save([paths.results.neuralRDMS.ITPlts filelistSess{sessi}(1:end-4) '_ICI.mat'], 'ICI');



end





%% plot mean ICI 


clear
paths = load_paths_WM('vvs', 'none');
cd(paths.results.neuralRDMS.ITPlts);

filelistSess = dir('*mat'); filelistSess = [{filelistSess.name}'];
for sessi= 4 %:length(filelistSess) 
    
    disp(['File > ' num2str(sessi)]);

    clearvars -except paths filelistSess sessi
    load(filelistSess{sessi});   
    
    mICI = squeeze(mean(ICI, 'omitnan')); 
    
    if size(mICI, 2) == 21
        times = 1:210;
        freqs = 1:520; 
        figure; set(gcf, 'Position', [100 100 200 400])
        contourf(times, freqs, myresizem(mICI, 10), 100, 'linecolor', 'none'); hold on; colorbar
        %contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 4);
        
        set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
        set(gca, 'xlim', [1 150], 'FontSize', 10); 
        %set(gca, 'clim', [-.5 .5], 'FontSize', 10);
        plot([45 45],get(gca,'ylim'), 'k:','lineWidth', 5); 
        %colorbar
    else
        times = 1:460;
        freqs = 1:520; 
        figure; set(gcf, 'Position', [1000 918 560 420])
        myCmap = colormap(brewermap([],'*Spectral'));
        colormap(myCmap)
        contourf(times, freqs, myresizem(mICI, 10), 100, 'linecolor', 'none'); hold on; colorbar
        %contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 4);
        
        set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
        set(gca, 'xlim', [1 390], 'FontSize', 10);
        set(gca, 'clim', [-.05 .05], 'FontSize', 10);
        plot([45 45],get(gca,'ylim'), 'k:','lineWidth', 5); 
        %colorbar


    end
    myCmap = colormap(brewermap([],'*Spectral'));
    colormap(myCmap)
    
    exportgraphics(gcf, [num2str(sessi) '_myP.png'], 'Resolution', 300); 


end




%% COMPUTE PERMUTATIONS

clear

nPerm = 100; 
paths = load_paths_WM('vvs', 'none');
filelistSess = getFilesWM(paths.results.neuralRDMS.IT);

for sessi= 1 %1:length(filelistSess) 
    disp(['File > ' num2str(sessi)]);
    clearvars -except paths filelistSess sessi nPerm

    load([paths.results.neuralRDMS.IT filelistSess{sessi}]);   
    nSubjs = size(allNeuralRDMS, 1); 
    nFreqs = size(allNeuralRDMS{1}, 3); 
    nTimes = size(allNeuralRDMS{1}, 4); 
    
    clear ICIPerm
    for subji = 1:nSubjs
        subji
        neuralRDMs  = allNeuralRDMS{subji, 1}; 
        ids         = allNeuralRDMS{subji, 2};

        for permi = 1:nPerm
            id4perm = randperm(length(ids)); 
            idsPERM = ids(id4perm);
            
            ITM = squeeze(load_ITMODEL_activ(idsPERM));
            for freqi  = 1:nFreqs
                for timei = 1:nTimes
                    neuralRDM = neuralRDMs(:, :, freqi, timei); 
                    W = neuralRDM(ITM==1); 
                    B = neuralRDM(ITM==0); 
                    ICIPerm(permi, subji, freqi, timei,:) = mean(W) - mean(B); 
                end
            end
        end
    end


    save([paths.results.neuralRDMS.ITPlts filelistSess{sessi}(1:end-4) '_ICI_' num2str(nPerm) 'p.mat'], 'ICIPerm');

end

%% plot with OUTLINE FROM PERM


clear
nPerm = 100;
paths = load_paths_WM('vvs', 'none');
cd(paths.results.neuralRDMS.ITPlts);

filelistSess = dir('*mat'); filelistSess = [{filelistSess.name}'];
for sessi= 3 %:length(filelistSess) %this one starts at 1 and not at 3
    
    disp(['File > ' num2str(sessi)]);

    clearvars -except paths filelistSess sessi nPerm
    load(filelistSess{sessi});   
    load(filelistSess{sessi+1});   
    
    mICI = squeeze(mean(ICI, 'omitnan')); 
    mICIPERM = squeeze(mean(ICIPerm, 2, 'omitnan')); 
    
    nFreq = size(mICI,1); nTimes = size(mICI, 2); 
    clear p
    for freqi = 1:nFreq
        for timei = 1:nTimes
            allICIPerm = mICIPERM(:, freqi, timei); 
            obsICI = mICI(freqi, timei); 
            allAB = allICIPerm(allICIPerm > obsICI);
            p(freqi, timei) = 1 - ((nPerm-1) - (length (allAB)))  / nPerm;

        end
    end
    h = p<.05; 
    
    

    
    if size(mICI, 2) == 21
        times = 1:210;
        freqs = 1:520; 
        figure; set(gcf, 'Position', [100 100 200 400])
        contourf(times, freqs, myresizem(mICI, 10), 100, 'linecolor', 'none'); hold on; colorbar
        contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 4);
        
        set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
        set(gca, 'xlim', [1 150], 'FontSize', 10); 
        set(gca, 'clim', [-.015 .015], 'FontSize', 10);
        plot([45 45],get(gca,'ylim'), 'k:','lineWidth', 5); 
        %colorbar
    else
        times = 1:460;
        freqs = 1:520; 
        figure; set(gcf, 'Position', [1000 918 560 420])
        myCmap = colormap(brewermap([],'*Spectral'));
        colormap(myCmap)
        contourf(times, freqs, myresizem(mICI, 10), 100, 'linecolor', 'none'); hold on; colorbar
        contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 4);
        
        set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
        set(gca, 'xlim', [1 390], 'FontSize', 10);
        set(gca, 'clim', [-.045 .045], 'FontSize', 10);
        plot([45 45],get(gca,'ylim'), 'k:','lineWidth', 5); 
        %colorbar


    end
    myCmap = colormap(brewermap([],'*Spectral'));
    colormap(myCmap)
    
    exportgraphics(gcf, [num2str(sessi) '_myP.png'], 'Resolution', 300); 


end










%% Partial correlation between neural RDMs and DNNs after partialling out the categorical model IN LOOP 
clear, clc
%Network_ROI_EoM_layers_freqs_avRepet_avTimeFeatVect_freqResolv(0-1)__fitMode(0:noTrials; 1:Trials)__timeRes__win-width__mf
    
listF2sav = {

'BLNETi_vvs_M123_[8-8-56]_3-54_1_0_1_0_.1_5_1_PC';
'BLNETi_vvs_E123_[8-8-56]_3-54_1_0_1_0_.1_5_1_PC';
'BLNETi_pfc_M123_[8-8-56]_3-54_1_0_1_0_.1_5_1_PC';
'BLNETi_pfc_E123_[8-8-56]_3-54_1_0_1_0_.1_5_1_PC';

'CORrt_vvs_M123_[2-2-8]_3-54_1_0_1_0_.1_5_1_PC';
'CORrt_vvs_E123_[2-2-8]_3-54_1_0_1_0_.1_5_1_PC';
'CORrt_pfc_M123_[2-2-8]_3-54_1_0_1_0_.1_5_1_PC';
'CORrt_pfc_E123_[2-2-8]_3-54_1_0_1_0_.1_5_1_PC';

};   


for listi = 1:length(listF2sav)
    disp(['File > ' num2str(listi) '      ' listF2sav{listi}]);
    clearvars -except listF2sav listi 
        
    f2sav       = listF2sav{listi}; 
    cfg = getParams(f2sav);
    cfg.DNN_analysis = 1; 
    paths = load_paths_WM(cfg.brainROI, cfg.net2load);
    filelistSess = getFiles(paths.traces);
    
    for sessi= 1:length(filelistSess) %this one starts at 1 and not at 3
        %disp(['File > ' num2str(sessi)]);
        load([paths.traces filelistSess{sessi}]);   
    
        cfg_contrasts               = getIdsWM(cfg.period, cfg_contrasts);
    
        if length(cfg_contrasts.oneListIds) > 1
            cfg_contrasts.oneListPow    = extract_power_WM (cfg_contrasts, cfg); % 
        
            cfg_contrasts               = normalize_WM(cfg_contrasts, 1, 'sess', []);
        
            cfg_contrasts               = average_repetitions(cfg, cfg_contrasts);
            
            neuralRDMs                  = createNeuralRDMs(cfg, cfg_contrasts);
            networkRDMs                 = createNetworkRDMs(cfg, cfg_contrasts, sessi, paths);
        
            nnFit{sessi,1}              = fitModelPartialCorrelation(cfg_contrasts, neuralRDMs, networkRDMs, cfg.fitMode); 
            nnFit{sessi,2}              = cfg_contrasts.oneListIds; 
        end

    end

    
    save([paths.results.DNNs f2sav '.mat'], 'nnFit');

end

%%  plot all layers FREQUENCY RESOLVED
%Network_ROI_ER_layers_freqs_avRepet_avTFV_fRes(0-1)_fitMode(0:noTrials; 1:Trials)_timeRes_win_mf
clear , clc

f2sav = 'CORrt_vvs_E123_[2-2-8]_3-54_1_0_1_0_.1_5_1_PC';

% 'Alex_vvs_E123_[1-8]_3-54_1_0_1_0_.1_5_1_PC';
% 'BLNETi_vvs_M123_[8-8-56]_3-54_1_0_1_0_.1_5_1_PC';
% 'BLNETi_vvs_E123_[8-8-56]_3-54_1_0_1_0_.1_5_1_PC';
% 'BLNETi_pfc_M123_[8-8-56]_3-54_1_0_1_0_.1_5_1_PC';
% 'BLNETi_pfc_E123_[8-8-56]_3-54_1_0_1_0_.1_5_1_PC';
% 
% 'CORrt_vvs_M123_[2-2-8]_3-54_1_0_1_0_.1_5_1_PC';
% 'CORrt_vvs_E123_[2-2-8]_3-54_1_0_1_0_.1_5_1_PC';
% 'CORrt_pfc_M123_[2-2-8]_3-54_1_0_1_0_.1_5_1_PC';
% 'CORrt_pfc_E123_[2-2-8]_3-54_1_0_1_0_.1_5_1_PC';



cfg = getParams(f2sav);
if strcmp(cfg.brainROI, 'vvs')
    sub2exc = [18 22];
elseif strcmp(cfg.brainROI, 'pfc')
    sub2exc = [1];
elseif strcmp(cfg.brainROI, 'hipp')
    sub2exc = [2];
end

paths = load_paths_WM(cfg.brainROI, cfg.net2load);
load([paths.results.DNNs f2sav '.mat']);


tiledlayout(8,9, 'TileSpacing', 'tight', 'Padding', 'compact');
if strcmp(cfg.period(1), 'M')
    set(gcf, 'Position', [100 100 1200 1000])
else
    set(gcf, 'Position', [100 100 700 1200])
end

for layi = 1:size(nnFit{1}, 1)
    ax1 = nexttile;
    clear nnH
    for subji = 1:length(nnFit)
       if strcmp(cfg.period(1), 'M')
         nnH(subji, : ,:) = atanh(nnFit{subji, 1}(layi,:,1:40));
       else
         nnH(subji, : ,:) = atanh(nnFit{subji, 1}(layi,:,1:15));
       end
    end
    
    nnH(sub2exc, :, :) = []; 
    nnH = squeeze(nnH);
    %[h p ci ts] = ttest(nnH, 0, "Tail","right");
    [h p ci ts] = ttest(nnH);
    h = squeeze(h); t = squeeze(ts.tstat); 
    h(:, 1:5) = 0; % only sum p-values in clusters after the baseline
    
    d2p = squeeze(mean(nnH, 'omitnan'));
    freqs = 1:52; 
    clustinfo = bwconncomp(h);
    allClustInfo{layi} = clustinfo; 

    % store allTOBS
    if ~isempty(clustinfo.PixelIdxList)
        for pixi = 1:length(clustinfo.PixelIdxList)
             %if length(clustinfo.PixelIdxList{pixi}) > 1
                allTObs(layi, pixi, :) = sum(t(clustinfo.PixelIdxList{pixi}));
             %end        
        end
    else
        allTObs(layi, :, :) = 0;
    end


    
    if strcmp(cfg.period(1), 'M')
        times = 1:size(t, 2); 
    else
        times = 1:15; 
    end
    myCmap = colormap(brewermap([],'YlOrRd'));
    colormap(myCmap)
    contourf(times, freqs, t, 100, 'linecolor', 'none'); hold on; %colorbar
    contour(times, freqs, h, 1, 'Color', [0, 0, 0], 'LineWidth', 2);
    
    if strcmp(cfg.period(1), 'M')
        set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
        set(gca, 'xlim', [1 40], 'clim', [-5 5], 'FontSize', 10);
        plot([5 5],get(gca,'ylim'), 'k:','lineWidth', 2);
    else
        set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
        set(gca, 'FontSize', 8, 'clim', [-5 5]);
        plot([5 5],get(gca,'ylim'), 'k:','lineWidth', 2);
        
    end
    

end

exportgraphics(gcf, [paths.results.DNNs 'myP.png'], 'Resolution', 300); 




%% PERMUTATIONS IN LOOP
%Network_ROI_ER_layers_freqs_avRepet_avTFV_fRes(0-1)_fitMode(0:noTrials; 1:Trials)__timeRes__win__mf

clear
nPerm = 1000;

listF2sav = {

'BLNETi_vvs_M123_[8-8-56]_3-54_1_0_1_0_.1_5_1_PC';
'BLNETi_vvs_E123_[8-8-56]_3-54_1_0_1_0_.1_5_1_PC';
'CORrt_vvs_M123_[2-2-8]_3-54_1_0_1_0_.1_5_1_PC';
'CORrt_vvs_E123_[2-2-8]_3-54_1_0_1_0_.1_5_1_PC';

'CORrt_pfc_M123_[2-2-8]_3-54_1_0_1_0_.1_5_1_PC';
'CORrt_pfc_E123_[2-2-8]_3-54_1_0_1_0_.1_5_1_PC';
'BLNETi_pfc_M123_[8-8-56]_3-54_1_0_1_0_.1_5_1_PC';
'BLNETi_pfc_E123_[8-8-56]_3-54_1_0_1_0_.1_5_1_PC';

                    
             };
    

for listi = 1:length(listF2sav)
    
    clearvars -except listF2sav listi nPerm 
        
    f2sav       = listF2sav{listi}
        
    cfg = getParams(f2sav);
    cfg.DNN_analysis = 1; 
    paths = load_paths_WM(cfg.brainROI, cfg.net2load);
    filelistSess = getFiles(paths.traces);
    
    t1 = datetime; 
    
    for sessi= 1:length(filelistSess) 
        disp(['File > ' num2str(sessi)]);
        load([paths.traces filelistSess{sessi}]);   

        nChans = size(cfg_contrasts.chanNames, 1); 
    
        if nChans > 1

            cfg_contrasts = getIdsWM(cfg.period, cfg_contrasts);

            if length(cfg_contrasts.oneListIds) > 1
                cfg_contrasts.oneListPow    = extract_power_WM (cfg_contrasts, cfg); % 
                
                cfg_contrasts = normalize_WM(cfg_contrasts, 1, 'sess', []);
                cfg_contrasts               = average_repetitions(cfg, cfg_contrasts);
                
                neuralRDMs                  = createNeuralRDMs(cfg, cfg_contrasts);
                networkRDMs                 = createNetworkRDMs(cfg, cfg_contrasts, sessi, paths);
        
        
                % % % restrict time for permutation data
                if strcmp(cfg.period(1), 'M')
                    if ndims(neuralRDMs) == 4
                        neuralRDMs = neuralRDMs(:,:,:,6:40); %frequency-resolved
                    else
                        neuralRDMs = neuralRDMs(:,:,6:40); %band analysis
                    end
                else
                    if ndims(neuralRDMs) == 4
                        neuralRDMs = neuralRDMs(:,:,:,6:15); %frequency-resolved
                    else
                        neuralRDMs = neuralRDMs(:,:,6:15); %band analysis
                    end
                end
                
                for permi = 1:nPerm
                    sC = size(networkRDMs, 2);
                    ids = randperm(sC);
                    networkRDMs = networkRDMs(:, ids, ids); 
                    nnFitPerm(permi, sessi,:,:, :)              = fitModelPartialCorrelation(cfg_contrasts, neuralRDMs, networkRDMs, cfg.fitMode); 
                end
            end
        end

    end
    
    mkdir ([paths.results.DNNs]);
    save([paths.results.DNNs f2sav '_' num2str(nPerm) 'p.mat'], 'nnFitPerm', '-v7.3');
    

end
   
t2 = datetime; 
etime(datevec(t2), datevec(t1))




%% COMPUTE ONE MEAN ITEM LEVEL DNN FIT PER SUBJECT AND PERFORM STATS


clear, clc
%Network_ROI_EoM_layers_freqs_avRepet_avTimeFeatVect_freqResolv(0-1)__fitMode(0:noTrials; 1:Trials)__timeRes__win-width__mf
    
listF2sav = {


% 'BLNETi_vvs_E123CAT1_[8-8-56]_3-54_0_0_1_0_.1_5_1';
% 'BLNETi_vvs_E123CAT2_[8-8-56]_3-54_0_0_1_0_.1_5_1';
% 'BLNETi_vvs_E123CAT3_[8-8-56]_3-54_0_0_1_0_.1_5_1';
% 'BLNETi_vvs_E123CAT4_[8-8-56]_3-54_0_0_1_0_.1_5_1';
% 'BLNETi_vvs_E123CAT5_[8-8-56]_3-54_0_0_1_0_.1_5_1';
% 'BLNETi_vvs_E123CAT6_[8-8-56]_3-54_0_0_1_0_.1_5_1';

% 'BLNETi_vvs_M123CAT1_[8-8-56]_3-54_0_0_1_0_.1_5_1';
% 'BLNETi_vvs_M123CAT2_[8-8-56]_3-54_0_0_1_0_.1_5_1';
% 'BLNETi_vvs_M123CAT3_[8-8-56]_3-54_0_0_1_0_.1_5_1';
% 'BLNETi_vvs_M123CAT4_[8-8-56]_3-54_0_0_1_0_.1_5_1';
% 'BLNETi_vvs_M123CAT5_[8-8-56]_3-54_0_0_1_0_.1_5_1';
% 'BLNETi_vvs_M123CAT6_[8-8-56]_3-54_0_0_1_0_.1_5_1';

% 'Alex_vvs_E123CAT1_[1-8]_3-54_0_0_1_0_.1_5_1';
% 'Alex_vvs_E123CAT2_[1-8]_3-54_0_0_1_0_.1_5_1';
% 'Alex_vvs_E123CAT3_[1-8]_3-54_0_0_1_0_.1_5_1';
% 'Alex_vvs_E123CAT4_[1-8]_3-54_0_0_1_0_.1_5_1';
% 'Alex_vvs_E123CAT5_[1-8]_3-54_0_0_1_0_.1_5_1';
% 'Alex_vvs_E123CAT6_[1-8]_3-54_0_0_1_0_.1_5_1';

% 'Alex_pfc_E123CAT1_[1-8]_3-54_0_0_1_0_.1_5_1';
% 'Alex_pfc_E123CAT2_[1-8]_3-54_0_0_1_0_.1_5_1';
% 'Alex_pfc_E123CAT3_[1-8]_3-54_0_0_1_0_.1_5_1';
% 'Alex_pfc_E123CAT4_[1-8]_3-54_0_0_1_0_.1_5_1';
% 'Alex_pfc_E123CAT5_[1-8]_3-54_0_0_1_0_.1_5_1';
% 'Alex_pfc_E123CAT6_[1-8]_3-54_0_0_1_0_.1_5_1';


'Alex_vvs_M123CAT1_[1-8]_3-54_0_0_1_0_.1_5_1';
'Alex_vvs_M123CAT2_[1-8]_3-54_0_0_1_0_.1_5_1';
'Alex_vvs_M123CAT3_[1-8]_3-54_0_0_1_0_.1_5_1';
'Alex_vvs_M123CAT4_[1-8]_3-54_0_0_1_0_.1_5_1';
'Alex_vvs_M123CAT5_[1-8]_3-54_0_0_1_0_.1_5_1';
'Alex_vvs_M123CAT6_[1-8]_3-54_0_0_1_0_.1_5_1';

% 'Alex_pfc_M123CAT1_[1-8]_3-54_0_0_1_0_.1_5_1';
% 'Alex_pfc_M123CAT2_[1-8]_3-54_0_0_1_0_.1_5_1';
% 'Alex_pfc_M123CAT3_[1-8]_3-54_0_0_1_0_.1_5_1';
% 'Alex_pfc_M123CAT4_[1-8]_3-54_0_0_1_0_.1_5_1';
% 'Alex_pfc_M123CAT5_[1-8]_3-54_0_0_1_0_.1_5_1';
% 'Alex_pfc_M123CAT6_[1-8]_3-54_0_0_1_0_.1_5_1';

% 'ITM_vvs_E123CAT1_[1]_3-54_0_0_1_0_.1_5_1';
% 'ITM_vvs_E123CAT2_[1]_3-54_0_0_1_0_.1_5_1';
% 'ITM_vvs_E123CAT3_[1]_3-54_0_0_1_0_.1_5_1';
% 'ITM_vvs_E123CAT4_[1]_3-54_0_0_1_0_.1_5_1';
% 'ITM_vvs_E123CAT5_[1]_3-54_0_0_1_0_.1_5_1';
% 'ITM_vvs_E123CAT6_[1]_3-54_0_0_1_0_.1_5_1';

% 'ITM_pfc_E123CAT1_[1]_3-54_0_0_1_0_.1_5_1';
% 'ITM_pfc_E123CAT2_[1]_3-54_0_0_1_0_.1_5_1';
% 'ITM_pfc_E123CAT3_[1]_3-54_0_0_1_0_.1_5_1';
% 'ITM_pfc_E123CAT4_[1]_3-54_0_0_1_0_.1_5_1';
% 'ITM_pfc_E123CAT5_[1]_3-54_0_0_1_0_.1_5_1';
% 'ITM_pfc_E123CAT6_[1]_3-54_0_0_1_0_.1_5_1';
% % 
% 'ITM_vvs_M123CAT1_[1]_3-54_0_0_1_0_.1_5_1';
% 'ITM_vvs_M123CAT2_[1]_3-54_0_0_1_0_.1_5_1';
% 'ITM_vvs_M123CAT3_[1]_3-54_0_0_1_0_.1_5_1';
% 'ITM_vvs_M123CAT4_[1]_3-54_0_0_1_0_.1_5_1';
% 'ITM_vvs_M123CAT5_[1]_3-54_0_0_1_0_.1_5_1';
% 'ITM_vvs_M123CAT6_[1]_3-54_0_0_1_0_.1_5_1';
% % 
% 'ITM_pfc_M123CAT1_[1]_3-54_0_0_1_0_.1_5_1';
% 'ITM_pfc_M123CAT2_[1]_3-54_0_0_1_0_.1_5_1';
% 'ITM_pfc_M123CAT3_[1]_3-54_0_0_1_0_.1_5_1';
% 'ITM_pfc_M123CAT4_[1]_3-54_0_0_1_0_.1_5_1';
% 'ITM_pfc_M123CAT5_[1]_3-54_0_0_1_0_.1_5_1';
% 'ITM_pfc_M123CAT6_[1]_3-54_0_0_1_0_.1_5_1';


};   


for listi = 1:length(listF2sav)
    disp(['File > ' num2str(listi) '      ' listF2sav{listi}]);
    clearvars -except listF2sav listi allnnFIT
        
    f2sav       = listF2sav{listi}; 
    cfg = getParams(f2sav);
    paths = load_paths_WM(cfg.brainROI, cfg.net2load);
    load([paths.results.DNNs f2sav '.mat']);
    x =  nnFit(:, 1);
    lay2u = 8; 
    x = cellfun(@(x) x(lay2u,:,:), x, 'un', 0); 
    allnnFIT{listi,:} = x; 

end


allNF = cat(2, allnnFIT{:}); 
allNF = cellfun(@(x) squeeze(x), allNF, 'un', 0);
allSNF = cell2mat(permute(allNF, [3 4 1 2])); 
size(allSNF);
mAallSNF = squeeze(mean(allSNF, 4));
mAallSNF = permute(mAallSNF, [3 1 2]); 

ax1 = nexttile;

if strcmp(cfg.period(1), 'M')
 nnH = mAallSNF(:,:,1:40);
else
 nnH = mAallSNF(:,:,1:15);
end

if strcmp(cfg.brainROI, 'vvs')
    sub2exc = [18 22];
elseif strcmp(cfg.brainROI, 'pfc')
    sub2exc = [1];
end


nnH(sub2exc, :, :) = []; 
nnH = squeeze(nnH);

% for subji = 1:size(nnH, 1)
%     nnHS = squeeze(nnH(subji, :, :));
%     mB = squeeze(mean(nnHS(:, 1:5),2));
%     nnHS = nnHS - mB; 
%     nnH(subji, :, :) = nnHS; 
% end

[h p ci ts] = ttest(nnH);
h = squeeze(h); t = squeeze(ts.tstat); 
%h(:, 1:5) = 0; % only sum p-values in clusters after the baseline

d2p = squeeze(mean(nnH, 'omitnan'));
freqs = 1:52; 
clustinfo = bwconncomp(h);

if strcmp(cfg.period(1), 'M')
    times = 1:size(t, 2); 
else
    times = 1:15; 
end
myCmap = colormap(brewermap([],'YlOrRd'));
colormap(myCmap)
contourf(times, freqs, t, 100, 'linecolor', 'none'); hold on; %colorbar
contour(times, freqs, h, 1, 'Color', [0, 0, 0], 'LineWidth', 2);
colorbar

if strcmp(cfg.period(1), 'M')
    set(gcf, 'Position', [100 100 1200 1000])
    set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
    set(gca, 'xlim', [1 40], 'FontSize', 10); %, 'clim', [-2 2]
    plot([5 5],get(gca,'ylim'), 'k:','lineWidth', 2);
else
    set(gcf, 'Position', [100 100 200 300])
    set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
    %set(gca, 'FontSize', 8, 'clim', [-5 5]);
    plot([5 5],get(gca,'ylim'), 'k:','lineWidth', 2);
    
end

exportgraphics(gcf, ['myP.png'], 'Resolution', 300); 

%% plot all nnH

for subji = 1:size(nnH, 1)
    
    figure()
    imagesc(squeeze(nnH(subji, :, :))); colorbar

end







%% Ttest at every time-frequency point comparing VVS and PFC fits 

%Network_ROI_ER_layers_freqs_avRepet_avTFV_fRes(0-1)_fitMode(0:noTrials; 1:Trials)_timeRes_win_mf
clear , clc
%f2sav1 = 'Alex_vvs_E123_[1-8]_3-54_1_0_1_0_.1_5_1';
%f2sav1 = 'BLNETi_vvs_E123_[8-8-56]_3-54_1_0_1_0_.1_5_1';
f2sav1 = 'CORrt_vvs_E123_[2-2-8]_3-54_1_0_1_0_.1_5_1';
cfg = getParams(f2sav1);
paths = load_paths_WM(cfg.brainROI, cfg.net2load);
load([paths.results.DNNs f2sav1 '.mat']);
nnFitVVS = nnFit; 

%f2sav2 = 'Alex_pfc_E123_[1-8]_3-54_1_0_1_0_.1_5_1';
%f2sav2 = 'BLNETi_pfc_E123_[8-8-56]_3-54_1_0_1_0_.1_5_1';
f2sav2 = 'CORrt_pfc_E123_[2-2-8]_3-54_1_0_1_0_.1_5_1';
cfg = getParams(f2sav2);
paths = load_paths_WM(cfg.brainROI, cfg.net2load);
load([paths.results.DNNs f2sav2 '.mat']);
nnFitPFC = nnFit; 



tiledlayout(8,9, 'TileSpacing', 'tight', 'Padding', 'compact');
if strcmp(cfg.period(1), 'M')
    set(gcf, 'Position', [100 100 1200 1000])
else
    set(gcf, 'Position', [100 100 700 1200])
end

for layi = 1 %:size(nnFit{1}, 1)
    ax1 = nexttile;
    clear nnHVVS nnHPFC
    for subji = 1:length(nnFitVVS)
       if strcmp(cfg.period(1), 'M')
         nnHVVS(subji, : ,:) = atanh(nnFitVVS{subji, 1}(layi,:,1:40));
       else
         nnHVVS(subji, : ,:) = atanh(nnFitVVS{subji, 1}(layi,:,1:15));
       end
    end
    for subji = 1:length(nnFitPFC)
       if strcmp(cfg.period(1), 'M')
         nnHPFC(subji, : ,:) = atanh(nnFitPFC{subji, 1}(layi,:,1:40));
       else
         nnHPFC(subji, : ,:) = atanh(nnFitPFC{subji, 1}(layi,:,1:15));
       end
    end
    
    nnHVVS([2 18], :, :) = []; 
    nnHPFC([1], :, :) = []; 
    nnHVVS = squeeze(nnHVVS);
    nnHPFC = squeeze(nnHPFC);

    [h p ci ts] = ttest2(nnHVVS, nnHPFC);
    h = squeeze(h); t = squeeze(ts.tstat); 
    h(:, 1:5) = 0; % only sum p-values in clusters after the baseline
    
    freqs = 1:52; 
    clustinfo = bwconncomp(h);
    allClustInfo{layi} = clustinfo; 

    % store allTOBS
    if ~isempty(clustinfo.PixelIdxList)
        for pixi = 1:length(clustinfo.PixelIdxList)
             allTObs(layi, pixi, :) = sum(t(clustinfo.PixelIdxList{pixi}));
             allTObs1(pixi, :) = sum(t(clustinfo.PixelIdxList{pixi}));
        end
    else
        allTObs(layi, :, :) = 0;
    end

    [max2u id] = max(abs(allTObs1));
    tObs = allTObs1(id); 


    if strcmp(cfg.period(1), 'M')
        times = 1:size(t, 2); 
    else
        times = 1:15; 
    end
    myCmap = colormap(brewermap([],'YlOrRd'));
    colormap(myCmap)
    contourf(times, freqs, t, 100, 'linecolor', 'none'); hold on; %colorbar
    contour(times, freqs, h, 1, 'Color', [0, 0, 0], 'LineWidth', 2);
    
    if strcmp(cfg.period(1), 'M')
        set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
        set(gca, 'xlim', [1 40], 'clim', [-5 5], 'FontSize', 10);
        plot([5 5],get(gca,'ylim'), 'k:','lineWidth', 2);
    else
        set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
        set(gca, 'FontSize', 8, 'clim', [-5 5]);
        plot([5 5],get(gca,'ylim'), 'k:','lineWidth', 2);
        
    end
    

end

exportgraphics(gcf, [paths.results.DNNs 'myP.png'], 'Resolution', 300); 


%% permutations
nPerm = 1000; 

junts = cat(1, nnHVVS, nnHPFC);
idsPerm = [zeros(1, size(nnHVVS, 1)), ones(1, size(nnHPFC, 1))]'; 

for permi = 1:nPerm

    idsPerm = idsPerm(randperm(length(idsPerm))); 
    nnHVVSPerm = junts(idsPerm == 0, :,:); 
    nnHPFCPerm = junts(idsPerm == 1, :,:); 

    [hPerm p ci ts] = ttest2(nnHVVSPerm, nnHPFCPerm);
    hPerm = squeeze(hPerm); tPerm = squeeze(ts.tstat); 

    clustinfoPerm = bwconncomp(hPerm);
    clear allSTs
    for pxi = 1:length(clustinfoPerm.PixelIdxList)
       allSTs(pxi) = sum(tPerm(clustinfoPerm.PixelIdxList{pxi}));% 
    end
    
    if exist('allSTs')
        [max2u id] = max(abs(allSTs));
    else
        disp('does not exist')
        allSTs = 0; 
        id = 1;
    end
    
    max_clust_sum_perm(permi,:) = allSTs(id); 
    
    

end

allAB = max_clust_sum_perm(max_clust_sum_perm> abs(tObs));
p = 1 - ((nPerm-1) - (length (allAB)))  / nPerm;

disp(['p = ' num2str(p)])


%% Ttest at every time-frequency point comparing VVS and PFC ICIs MAINTENANCE

clear
paths = load_paths_WM('vvs', 'none');
cd(paths.results.neuralRDMS.ITPlts);

filelistSess = dir('*mat'); filelistSess = [{filelistSess.name}'];
load('vvs_M123_[1-8]_3-54_0_0_1_0_.1_5_1_ICI');
ICIVVS = ICI; 
load('pfc_M123_[1-8]_3-54_0_0_1_0_.1_5_1_ICI');
ICIPFC = ICI; 


ICIVVS([2 18], :, :) = []; 
ICIPFC([1], :, :) = []; 

[h p ci ts] = ttest2(ICIVVS, ICIPFC);
h = squeeze(h); t = squeeze(ts.tstat); 

freqs = 1:52; 
times = 1:46
clustinfo = bwconncomp(h);

for pxi = 1:length(clustinfo.PixelIdxList)
   allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
end        
[max2u id] = max(abs(allSTs));

tObs =  allSTs(id); 


myCmap = colormap(brewermap([],'YlOrRd'));
colormap(myCmap)
contourf(times, freqs, t, 100, 'linecolor', 'none'); hold on; %colorbar
contour(times, freqs, h, 1, 'Color', [0, 0, 0], 'LineWidth', 2);
set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
set(gca, 'xlim', [1 40], 'clim', [-5 5], 'FontSize', 10);
plot([5 5],get(gca,'ylim'), 'k:','lineWidth', 2);

%% permutations
nPerm = 1000; 

junts = cat(1, ICIVVS, ICIPFC);
idsPerm = [zeros(1, size(ICIVVS, 1)), ones(1, size(ICIPFC, 1))]'; 

for permi = 1:nPerm

    idsPerm = idsPerm(randperm(length(idsPerm))); 
    ICIVVSPerm = junts(idsPerm == 0, :,:); 
    ICIPFCPerm = junts(idsPerm == 1, :,:); 

    [hPerm p ci ts] = ttest2(ICIVVSPerm, ICIPFCPerm);
    hPerm = squeeze(hPerm); tPerm = squeeze(ts.tstat); 

    clustinfoPerm = bwconncomp(hPerm);
    clear allSTs
    for pxi = 1:length(clustinfoPerm.PixelIdxList)
       allSTs(pxi) = sum(tPerm(clustinfoPerm.PixelIdxList{pxi}));% 
    end
    
    if exist('allSTs')
        [max2u id] = max(abs(allSTs));
    else
        disp('does not exist')
        allSTs = 0; 
        id = 1;
    end
    
    max_clust_sum_perm(permi,:) = allSTs(id); 
    
    

end

allAB = max_clust_sum_perm(max_clust_sum_perm> abs(tObs));
p = 1 - ((nPerm-1) - (length (allAB)))  / nPerm;

disp(['p = ' num2str(p)])


%% Ttest at every time-frequency point comparing VVS and PFC ICIs ENCODING

clear
paths = load_paths_WM('vvs', 'none');
cd(paths.results.neuralRDMS.ITPlts);

filelistSess = dir('*mat'); filelistSess = [{filelistSess.name}'];
load('vvs_E123_[1-8]_3-54_0_0_1_0_.1_5_1_ICI');
ICIVVS = ICI; 
load('pfc_E123_[1-8]_3-54_0_0_1_0_.1_5_1_ICI');
ICIPFC = ICI; 


ICIVVS([2 18], :, :) = []; 
ICIPFC([1], :, :) = []; 

[h p ci ts] = ttest2(ICIVVS, ICIPFC);
h = squeeze(h); t = squeeze(ts.tstat); 

freqs = 1:52; 
times = 1:21
clustinfo = bwconncomp(h);

myCmap = colormap(brewermap([],'YlOrRd'));
colormap(myCmap)
contourf(times, freqs, t, 100, 'linecolor', 'none'); hold on; %colorbar
contour(times, freqs, h, 1, 'Color', [0, 0, 0], 'LineWidth', 2);
set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
set(gca, 'xlim', [1 15], 'clim', [-5 5], 'FontSize', 10);
plot([5 5],get(gca,'ylim'), 'k:','lineWidth', 2);




































%%




































%%

