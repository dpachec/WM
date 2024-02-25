%% 
%% Export power from epoched rawTraces 
% % note that multi-item trials are discarded, and data is normalized across trials after removing these trials
clear 
paths = load_paths_WM('vvs', 'none');
filelistSess = getFilesWM(paths.traces);

for sessi= 1:length(filelistSess) %this one starts at 1 and not at 3
    disp(['Subj > ' num2str(sessi)]);
    clearvars -except sessi filelistSess paths cfg
    load([paths.traces filelistSess{sessi}]);   

    cfg.timeRes = .1; 
    cfg.DNN_analysis = 1; 
    cfg.period = 'M'; 
    oneListPow                  = extract_power_WM (cfg_contrasts, cfg); % 
    cfg_contrasts.oneListPow    = oneListPow; 
    
    cfg_contrastsE               = getIdsWM('E123', cfg_contrasts);
    cfg_contrastsE               = normalize_WM(cfg_contrastsE, 1, 'sess', []);

    cfg_contrastsM               = getIdsWM('M123', cfg_contrasts);
    cfg_contrastsM               = normalize_WM(cfg_contrastsM, 1, 'sess', []);

    cfg_contrasts = []; 
    cfg_contrasts.oneListPow = cat(1, cfg_contrastsE.oneListPow, cfg_contrastsM.oneListPow); 
    cfg_contrasts.chanNames = cfg_contrastsE.chanNames; 
    cfg_contrasts.oneListIds = cat(1, cfg_contrastsE.oneListIds, cfg_contrastsM.oneListIds); 



    save([filelistSess{sessi}(1:3) '_normPow.mat'], 'cfg_contrasts');
    
end

%% Calculate from epoched raw traces
clear
%Network_ROI_layers_freqs_avRepet_avTimeFeatVect_freqResolv(0-1)_fitMode(0:noTrials; 1:Trials)_timeRes_win-width_mf

%f2sav = 'Alex_pfc_E123_[1-8]_3-54_1_0_1_0_.1_5_1'; 
f2sav = 'BLNETi_pfc_M123_[8-8-56]_3-54_1_0_1_0_.1_5_1';


cfg = getParams(f2sav);
cfg.DNN_analysis = 1; 
paths = load_paths_WM(cfg.brainROI, cfg.net2load);
filelistSess = getFilesWM(paths.traces);

t1 = datetime; 

for sessi= 1:length(filelistSess) %this one starts at 1 and not at 3
    disp(['Subj > ' num2str(sessi)]);
    load([paths.traces filelistSess{sessi}]);   

    cfg_contrasts.oneListPow    = extract_power_WM (cfg_contrasts, cfg); % 
    cfg_contrasts               = normalize_WM(cfg_contrasts, 1, 'sess', []);
    cfg_contrasts               = getIdsWM(cfg.period, cfg_contrasts);

    if length(cfg_contrasts.oneListIds) > 1
        
        cfg_contrasts               = average_repetitions(cfg, cfg_contrasts);
        
        neuralRDMs                  = createNeuralRDMs(cfg, cfg_contrasts);
        networkRDMs                 = createNetworkRDMs(cfg, cfg_contrasts, sessi, paths);

   
        nnFit{sessi,1}              = fitModel_WM(neuralRDMs, networkRDMs, cfg.fitMode); 
        nnFit{sessi,2}              = cfg_contrasts.oneListIds; 
    end
end

mkdir ([paths.results.DNNs]);
save([paths.results.DNNs f2sav '.mat'], 'nnFit');

t2 = datetime; 
etime(datevec(t2), datevec(t1))



%% IN LOOP LOADING POWER DATA
clear, clc
%Network_ROI_EoM_layers_freqs_avRepet_avTimeFeatVect_freqResolv(0-1)__fitMode(0:noTrials; 1:Trials)__timeRes__win-width__mf
    
listF2sav = {

%'BLNETi_pfc_E123_[8-8-56]_3-54_0_0_1_1_.1_5_1';
'Alex_vvs_E123_[1-8]_3-54_0_0_1_0_.1_5_1'; 


};   

t1 = datetime; 
for listi = 1:length(listF2sav)
    disp(['File > ' num2str(listi) '      ' listF2sav{listi}]);
    clearvars -except listF2sav listi t1
        
    f2sav       = listF2sav{listi}; 
    cfg = getParams(f2sav);
    cfg.DNN_analysis = 1; 
    paths = load_paths_WM(cfg.brainROI, cfg.net2load);
    filelistSess = getFilesWM(paths.powerFromRT);
    
    for sessi= 1:length(filelistSess) 
        disp(['File > ' num2str(sessi)]);
        load([paths.powerFromRT filelistSess{sessi}]);   
        
        cfg_contrasts               = getIdsWM(cfg.period, cfg_contrasts);
    
        if length(cfg_contrasts.oneListIds) > 1 & size(cfg_contrasts.chanNames, 1) > 1
            cfg_contrasts               = average_repetitions(cfg, cfg_contrasts);
            neuralRDMs                  = createNeuralRDMs(cfg, cfg_contrasts);
            networkRDMs                 = createNetworkRDMs(cfg, cfg_contrasts, sessi, paths);
            if strcmp(cfg.meth, 'MASK')
                networkRDMs                 = restrictBetCatCorr(cfg_contrasts, networkRDMs); 
            end

            if ~strcmp(cfg.meth(1:2), 'PC')
                tic
                nnFit{sessi,1}              = fitModel_WM(neuralRDMs, networkRDMs, cfg.fitMode); 
                toc
            elseif strcmp(cfg.meth(1:2), 'PC')
                nnFit{sessi,1}              = fitModelPartialCorrelation(cfg_contrasts, neuralRDMs, networkRDMs, f2sav(end)); 
            end
            nnFit{sessi,2}              = cfg_contrasts.oneListIds; 
        end
    end
    
    save([paths.results.DNNs f2sav '.mat'], 'nnFit');

end
t2 = datetime; 
etime(datevec(t2), datevec(t1))




%% IN LOOP LOADING POWER DATA PERMUTATIONS
clear, clc
%Network_ROI_EoM_layers_freqs_avRepet_avTimeFeatVect_freqResolv(0-1)__fitMode(0:noTrials; 1:Trials)__timeRes__win-width__mf
    
nPerm = 1000; 
listF2sav = {


%'Alex_vvs_E123_[1-8]_3-54_0_0_1_0_.1_5_1_MASK'; 
'Alex_pfc_E123_[1-8]_3-54_0_0_1_0_.1_5_1_MASK'; 
'Alex_vvs_M123_[1-8]_3-54_0_0_1_0_.1_5_1_MASK'; 
'Alex_pfc_M123_[1-8]_3-54_0_0_1_0_.1_5_1_MASK'; 

'BLNETi_vvs_M123_[8-8-56]_3-54_0_0_1_0_.1_5_1_MASK'; 
'BLNETi_vvs_E123_[8-8-56]_3-54_0_0_1_0_.1_5_1_MASK'; 
'BLNETi_pfc_E123_[8-8-56]_3-54_0_0_1_0_.1_5_1_MASK'; 
'BLNETi_pfc_M123_[8-8-56]_3-54_0_0_1_0_.1_5_1_MASK'; 

'BLNETi_vvs_M123_[8-8-56]_3-54_1_0_1_0_.1_5_1_PCC';
'BLNETi_pfc_M123_[8-8-56]_3-54_1_0_1_0_.1_5_1_PCC';
'BLNETi_vvs_E123_[8-8-56]_3-54_1_0_1_0_.1_5_1_PCC';
'BLNETi_pfc_E123_[8-8-56]_3-54_1_0_1_0_.1_5_1_PCC';

'Alex_vvs_M123_[1-8]_3-54_1_0_1_0_.1_5_1_PCC';
'Alex_pfc_M123_[1-8]_3-54_1_0_1_0_.1_5_1_PCC';
'Alex_vvs_E123_[1-8]_3-54_1_0_1_0_.1_5_1_PCC';
'Alex_pfc_E123_[1-8]_3-54_1_0_1_0_.1_5_1_PCC';

'CORrt_vvs_M123_[2-2-8]_3-54_1_0_1_0_.1_5_1_PCC';
'CORrt_pfc_M123_[2-2-8]_3-54_1_0_1_0_.1_5_1_PCC';
'CORrt_vvs_E123_[2-2-8]_3-54_1_0_1_0_.1_5_1_PCC';
'CORrt_pfc_E123_[2-2-8]_3-54_1_0_1_0_.1_5_1_PCC';

'CORrt_vvs_M123_[2-2-8]_3-54_0_0_1_0_.1_5_1_MASK'; 
'CORrt_vvs_E123_[2-2-8]_3-54_0_0_1_0_.1_5_1_MASK'; 
'CORrt_pfc_E123_[2-2-8]_3-54_0_0_1_0_.1_5_1_MASK'; 
'CORrt_pfc_M123_[2-2-8]_3-54_0_0_1_0_.1_5_1_MASK'; 

'Alex_vvs_M11_[1-8]_3-54_1_0_1_0_.1_5_1'; 
'Alex_vvs_M12_[1-8]_3-54_1_0_1_0_.1_5_1'; 
'Alex_vvs_M13_[1-8]_3-54_1_0_1_0_.1_5_1'; 

'Alex_pfc_M11_[1-8]_3-54_1_0_1_0_.1_5_1'; 
'Alex_pfc_M12_[1-8]_3-54_1_0_1_0_.1_5_1'; 
'Alex_pfc_M13_[1-8]_3-54_1_0_1_0_.1_5_1'; 

'BLNETi_vvs_M11_[8-8-56]_3-54_1_0_1_0_.1_5_1'; 
'BLNETi_vvs_M12_[8-8-56]_3-54_1_0_1_0_.1_5_1'; 
'BLNETi_vvs_M13_[8-8-56]_3-54_1_0_1_0_.1_5_1'; 

'BLNETi_pfc_M11_[8-8-56]_3-54_1_0_1_0_.1_5_1'; 
'BLNETi_pfc_M12_[8-8-56]_3-54_1_0_1_0_.1_5_1'; 
'BLNETi_pfc_M13_[8-8-56]_3-54_1_0_1_0_.1_5_1'; 

};   

t1 = datetime; 

for listi = 1:length(listF2sav)
    disp(['File > ' num2str(listi) '      ' listF2sav{listi}]);
    clearvars -except listF2sav listi nPerm t1
        
    f2sav       = listF2sav{listi}; 
    cfg = getParams(f2sav);
    cfg.DNN_analysis = 1; 
    paths = load_paths_WM(cfg.brainROI, cfg.net2load);
    filelistSess = getFilesWM(paths.powerFromRT);
    
    for sessi= 1:length(filelistSess) %this one starts at 1 and not at 3
        disp(['File > ' num2str(sessi)]);
        load([paths.powerFromRT filelistSess{sessi}]);   
        
        cfg_contrasts               = getIdsWM(cfg.period, cfg_contrasts);
    
        if length(cfg_contrasts.oneListIds) > 1 & size(cfg_contrasts.chanNames, 1) > 1
            cfg_contrasts               = average_repetitions(cfg, cfg_contrasts);
            tic
            neuralRDMs                  = createNeuralRDMs(cfg, cfg_contrasts);
            toc
            networkRDMs                 = createNetworkRDMs(cfg, cfg_contrasts, sessi, paths);
            if strcmp(cfg.meth, 'MASK')
                networkRDMs                 = restrictBetCatCorr(cfg_contrasts, networkRDMs); 
            end
            neuralRDMs                  = restrictTime4Perm(cfg, neuralRDMs); 

            parfor permi = 1:nPerm
                sC = size(networkRDMs, 2);
                ids = randperm(sC);
                networkRDMs1 = networkRDMs(:, ids, ids); 
                if ~strcmp(cfg.meth(1:2), 'PC')
                    nnFitPerm(permi, sessi,:,:, :,:)     = fitModel_WM(neuralRDMs, networkRDMs1, cfg.fitMode); 
                elseif strcmp(cfg.meth(1:2), 'PC')
                    nnFitPerm(permi, sessi,:,:, :,:)     = fitModelPartialCorrelation(cfg_contrasts, neuralRDMs, networkRDMs, f2sav(end)); 
                end
            end
        end
    end
    
    save([paths.results.DNNs f2sav num2str(nPerm) 'p.mat'], 'nnFitPerm');

end

t2 = datetime; 
etime(datevec(t2), datevec(t1))


%% IN LOOP LOADING POWER DATA PERMUTATIONS for trial level analysis
clear, clc
%Network_ROI_EoM_layers_freqs_avRepet_avTimeFeatVect_freqResolv(0-1)__fitMode(0:noTrials; 1:Trials)__timeRes__win-width__mf
    
nPerm = 1000; 
listF2sav = {

'BLNETi_vvs_M123_[8-8-56]_3-54_0_0_1_1_.1_5_1';
'BLNETi_pfc_M123_[8-8-56]_3-54_0_0_1_1_.1_5_1';
'BLNETi_vvs_E123_[8-8-56]_3-54_0_0_1_1_.1_5_1';
'BLNETi_pfc_E123_[8-8-56]_3-54_0_0_1_1_.1_5_1';


};   

t1 = datetime; 

for listi = 1:length(listF2sav)
    disp(['File > ' num2str(listi) '      ' listF2sav{listi}]);
    clearvars -except listF2sav listi nPerm
        
    f2sav       = listF2sav{listi}; 
    cfg = getParams(f2sav);
    cfg.DNN_analysis = 1; 
    paths = load_paths_WM(cfg.brainROI, cfg.net2load);
    filelistSess = getFilesWM(paths.powerFromRT);
    
    for sessi= 1:length(filelistSess) %this one starts at 1 and not at 3
        disp(['Sessi > ' num2str(sessi)]);
        load([paths.powerFromRT filelistSess{sessi}]);   
        
        cfg_contrasts               = getIdsWM(cfg.period, cfg_contrasts);
    
        if length(cfg_contrasts.oneListIds) > 1 & size(cfg_contrasts.chanNames, 1) > 1
            cfg_contrasts               = average_repetitions(cfg, cfg_contrasts);
            neuralRDMs                  = createNeuralRDMs(cfg, cfg_contrasts);
            networkRDMs                 = createNetworkRDMs(cfg, cfg_contrasts, sessi, paths);
            neuralRDMs                  = restrictTime4Perm(cfg, neuralRDMs); 

            sC = size(networkRDMs, 2);
            ids = cellfun(@(x) strsplit(x, ' '), cfg_contrasts.oneListIds, 'un', 0); %str2num does not work
            ids = double(string(cat(1, ids{:})));%str2num does not work
          
            CC = find(ids(:,20)==1); 
            IC = find(ids(:,20)==0); 
            CI = find(ids(:,19)==1); 
            II = find(ids(:,19)==0); 
            
            parfor permi = 1:nPerm
                if length(CC) > 2
                    networkRDMsCC = networkRDMs(:, CC, CC); 
                    neuralRDMsCC = neuralRDMs(CC, CC, :, :); 
                    sC = size(neuralRDMsCC, 2);
                    ids = randperm(sC);
                    neuralRDMsCC = neuralRDMsCC(:, ids, ids); 
                    nnFitPerm{permi, sessi, 1} = fitModel_WM(neuralRDMsCC, networkRDMsCC, cfg.fitMode); 
                end
            end
            parfor permi = 1:nPerm
                if length(IC) > 2
                    networkRDMsIC = networkRDMs(:, IC, IC); 
                    neuralRDMsIC = neuralRDMs(IC, IC, :, :); 
                    sC = size(neuralRDMsIC, 2);
                    ids = randperm(sC);
                    neuralRDMsIC = neuralRDMsIC(:, ids, ids); 
                    nnFitPerm{permi, sessi, 2} = fitModel_WM(neuralRDMsIC, networkRDMsIC, cfg.fitMode); 
                end
            end
            parfor permi = 1:nPerm
                if length(CI) > 2
                    networkRDMsCI = networkRDMs(:, CI, CI); 
                    neuralRDMsCI = neuralRDMs(CI, CI, :, :);
                    sC = size(neuralRDMsCI, 2);
                    ids = randperm(sC);
                    neuralRDMsCI = neuralRDMsCI(:, ids, ids); 
                    nnFitPerm{permi, sessi, 3} = fitModel_WM(neuralRDMsCI, networkRDMsCI, cfg.fitMode); 
                end
            end
            parfor permi = 1:nPerm
                if length(II) > 2
                    networkRDMsII = networkRDMs(:, II, II); 
                    neuralRDMsII = neuralRDMs(II, II, :, :); 
                    sC = size(neuralRDMsII, 2);
                    ids = randperm(sC);
                    neuralRDMsII = neuralRDMsII(:, ids, ids); 
                    nnFitPerm{permi, sessi, 4} = fitModel_WM(neuralRDMsII, networkRDMsII, cfg.fitMode); 
                end                
            end
        end
    end
    
    save([paths.results.DNNs f2sav '.mat'], 'nnFitPerm');

end

t2 = datetime; 
etime(datevec(t2), datevec(t1))



%%  plot all layers FREQUENCY RESOLVED
%Network_ROI_ER_layers_freqs_avRepet_avTFV_fRes(0-1)_fitMode(0:noTrials; 1:Trials)_timeRes_win_mf
clear , clc

%f2sav = 'CAT_pfc_E123_[1]_3-54_1_0_1_0_.1_5_1';
%f2sav = 'Alex_vvs_E123_[1-8]_3-54_0_0_1_0_.1_5_1_MASK'
%f2sav = 'BLNETi_vvs_M123_[8-8-56]_3-54_0_0_1_0_.1_5_1_MASK'; 
%f2sav = 'Alex_vvs_E123_[1-8]_3-54_0_0_1_0_.1_5_1_MASK'; 
%f2sav = 'BLNETi_pfc_E123_[8-8-56]_3-54_1_0_1_0_.1_5_1_PCC'
f2sav = 'Alex_vvs_M13_[1-8]_3-54_1_0_1_0_.1_5_1'



cfg = getParams(f2sav);
if strcmp(cfg.brainROI, 'vvs')
    sub2exc = [18 22];
    %sub2exc = [18 22 23]; % for the incorrect trials
elseif strcmp(cfg.brainROI, 'pfc')
    sub2exc = [1];
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
       if ~isempty(nnFit{subji, 1})
           if strcmp(cfg.period(1), 'M')
             nnH(subji, : ,:) = atanh(nnFit{subji, 1}(layi,:,1:40));
           else
             nnH(subji, : ,:) = atanh(nnFit{subji, 1}(layi,:,1:15));
             %nnH(subji, : ,:) = atanh(nnFit{subji, 1}(layi,:,1:134));
           end
       end
    end
    
    %nnH(nnH==inf) = nan; 
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

%exportgraphics(gcf, [paths.results.DNNs 'myP.png'], 'Resolution', 300); 
exportgraphics(gcf, ['myP.png'], 'Resolution', 300); 




%% Extract activity in specific clusters from NNH


paths = load_paths_WM('vvs', 'none');
load ([paths.results.clusters 'ITM_VVS_enc']);

for subji = 1:size(nnH, 1)
    nnHSubj = squeeze(nnH(subji, :, :)); 
    nnHClust_vvsE1(subji, :) = mean(nnHSubj(clustinfo.PixelIdxList{6})); 
end

[h p ci ts] = ttest(nnHClust_vvsE1);

h = squeeze(h); t = squeeze(ts.tstat); 

disp (['t = ' num2str(ts.tstat) '  ' ' p = ' num2str(p)]);


%%  plot SEPARATELY FOR CORRECT AND INCORRECT starting from the trial level fits
%Network_ROI_ER_layers_freqs_avRepet_avTFV_fRes(0-1)_fitMode(0:noTrials; 1:Trials)_timeRes_win_mf
clear , clc

%f2sav = 'Alex_pfc_M123_[1-8]_3-54_0_0_1_1_.1_5_1';
f2sav = 'BLNETi_vvs_M123_[8-8-56]_3-54_0_0_1_1_.1_5_1';

cond2plot = 'CC'; 

cfg = getParams(f2sav);
if strcmp(cfg.brainROI, 'vvs')
    sub2exc = [18 22];
elseif strcmp(cfg.brainROI, 'pfc')
    sub2exc = [1 6 11];
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
       if ~isempty(nnFit{subji, 1})
           if strcmp(cfg.period(1), 'M')
             avTR = atanh(nnFit{subji, 1}(layi,:,:,1:40));
             %ids= str2num(cell2mat(nnFit{subji, 2}));
             ids = cellfun(@(x) strsplit(x, ' '), nnFit{subji, 2}, 'un', 0); %str2num does not work
             ids = double(string(cat(1, ids{:})));%str2num does not work
           elseif strcmp(cfg.period(1), 'E') % for clarity
             avTR = atanh(nnFit{subji, 1}(layi,:,:,1:15));
             ids= str2num(cell2mat(nnFit{subji, 2}));
           end
            
            if strcmp(cond2plot, 'CC')
                %ids_IC = ids(:,20)==0; 
                %nTR_IC = sum(ids_IC==1); 
                ids = find(ids(:,20)==1); 
                %if length(ids) > nTR_IC
                %    ids = ids(randperm(length(ids), nTR_IC), :);
                %end
            elseif strcmp(cond2plot, 'IC')
                ids = ids(:,20)==0; 
                nIncTR(subji, :) = sum(ids==1); 
            elseif strcmp(cond2plot, 'CI')
                ids = ids(:,19)==1; 
            elseif strcmp(cond2plot, 'II')
                ids = ids(:,19)==0; 
            elseif strcmp(cond2plot, 'allT')
                ids = 1:size(ids, 1); %just takes them all
            end

             nnH(subji, : ,:) = squeeze(mean(avTR(:, ids,:,:), 2)); 

       end
    end
    
    nnH(nnH==inf) = nan; 
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

%exportgraphics(gcf, [paths.results.DNNs 'myP.png'], 'Resolution', 300); 
exportgraphics(gcf, ['myP.png'], 'Resolution', 300); 



%%  COMPARE CORRECT VS INCORRECT
%Network_ROI_ER_layers_freqs_avRepet_avTFV_fRes(0-1)_fitMode(0:noTrials; 1:Trials)_timeRes_win_mf
clear , clc

%f2sav = 'Alex_vvs_E123_[1-8]_3-54_0_0_1_1_.1_5_1';
f2sav = 'BLNETi_vvs_M123_[8-8-56]_3-54_0_0_1_1_.1_5_1';

cfg = getParams(f2sav);
if strcmp(cfg.brainROI, 'vvs')
    sub2exc = [18 22];
    %sub2exc = [18 22 10 20 16];
elseif strcmp(cfg.brainROI, 'pfc')
    sub2exc = [1 6 11];
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
       if ~isempty(nnFit{subji, 1})
           if strcmp(cfg.period(1), 'M')
             avTR = atanh(nnFit{subji, 1}(layi,:,:,1:40));
             ids = cellfun(@(x) strsplit(x, ' '), nnFit{subji, 2}, 'un', 0); %str2num does not work
             ids = double(string(cat(1, ids{:})));%str2num does not work
           elseif strcmp(cfg.period(1), 'E') % for clarity
             avTR = atanh(nnFit{subji, 1}(layi,:,:,1:15));
             ids= str2num(cell2mat(nnFit{subji, 2}));
           end
            idsCC = ids(:,20)==1; 
            idsIC = ids(:,20)==0; 
            nIncTR(subji, :) = sum(idsIC==1); 
            idsCI = ids(:,19)==1; 
            idsII = ids(:,19)==0; 
            
            nnH1(subji, : ,:) = squeeze(mean(avTR(:, idsCC,:,:), 2)); 
            nnH2(subji, : ,:) = squeeze(mean(avTR(:, idsIC,:,:), 2)); 

       end
    end
    
    nnH1(nnH1==inf) = nan; 
    nnH1(sub2exc, :, :) = []; 
    nnH1 = squeeze(nnH1);

    nnH2(nnH2==inf) = nan; 
    nnH2(sub2exc, :, :) = []; 
    nnH2 = squeeze(nnH2);

    [h p ci ts] = ttest(nnH1, nnH2);
    h = squeeze(h); t = squeeze(ts.tstat); 
    %h(:, 1:5) = 0; % only sum p-values in clusters after the baseline
    
    d2p = squeeze(mean(nnH1 - nnH2, 'omitnan'));
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

%exportgraphics(gcf, [paths.results.DNNs 'myP.png'], 'Resolution', 300); 
exportgraphics(gcf, ['myP.png'], 'Resolution', 300); 





%% plot all layers BANDS > recover method for plotting the category model in individual bands

%Network_ROI_ER_layers_freqs_avRepet_avTFV_fRes(0-1)_fitMode(0:noTrials; 1:Trials)_timeRes_win_mf
clear , clc

f2sav = 'CAT_vvs_E123_[1-8]_3-8_1_0_0_0_.1_5_1';



cfg = getParams(f2sav);
if strcmp(cfg.brainROI, 'vvs')
    %sub2exc = [18 22];
    sub2exc = [18 22 10 20]; % for the incorrect trials
elseif strcmp(cfg.brainROI, 'pfc')
    sub2exc = [1];
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
       if ~isempty(nnFit{subji, 1})
           if strcmp(cfg.period(1), 'M')
             nnH(subji, : ,:) = atanh(nnFit{subji, 1}(layi,1:40));
           else
             nnH(subji, : ,:) = atanh(nnFit{subji, 1}(layi,1:15));
             %nnH(subji, : ,:) = atanh(nnFit{subji, 1}(layi,:,1:134));
           end
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
    plot(times, t); hold on; %colorbar
    
    if strcmp(cfg.period(1), 'M')
        
    else
        
    end
    

end

%exportgraphics(gcf, [paths.results.DNNs 'myP.png'], 'Resolution', 300); 
exportgraphics(gcf, ['myP.png'], 'Resolution', 300); 





%%  plot all layers FREQUENCY RESOLVED
%Network_ROI_ER_layers_freqs_avRepet_avTFV_fRes(0-1)_fitMode(0:noTrials; 1:Trials)_timeRes_win_mf
clear , clc

%f2sav = 'CAT_vvs_E123_[1-8]_3-8_1_0_0_0_.1_5_1';
f2sav =  'Alex_pfc_E123_[1-8]_3-54_0_0_1_0_.1_5_1_MASK'; 


cfg = getParams(f2sav);
if strcmp(cfg.brainROI, 'vvs')
    %sub2exc = [18 22];
    sub2exc = [18 22 10 20]; % for the incorrect trials
elseif strcmp(cfg.brainROI, 'pfc')
    sub2exc = [1];
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
       if ~isempty(nnFit{subji, 1})
           if strcmp(cfg.period(1), 'M')
             nnH(subji, : ,:) = atanh(nnFit{subji, 1}(layi,:,1:40));
           else
             nnH(subji, : ,:) = atanh(nnFit{subji, 1}(layi,:,1:15));
             %nnH(subji, : ,:) = atanh(nnFit{subji, 1}(layi,:,1:134));
           end
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

%exportgraphics(gcf, [paths.results.DNNs 'myP.png'], 'Resolution', 300); 
exportgraphics(gcf, ['myP.png'], 'Resolution', 300); 


%% Extract activity in specific clusters DURING MAINTENANCE (FOR PERFORMANCE; OR BL-NET FITS)
clear
%f2sav = 'BLNETi_vvs_M123_[8-8-56]_3-54_0_0_1_0_.1_5_1';
%f2sav = 'AlexEco_pfc_M123_[1-8]_3-54_1_0_1_0_.1_5_1';
f2sav = 'ITM_vvs_M123_[1]_3-54_0_0_1_0_.1_5_1';

cfg = getParams(f2sav);
if strcmp(cfg.brainROI, 'vvs')
    sub2exc = [18 22];
elseif strcmp(cfg.brainROI, 'pfc')
    sub2exc = [1];
end

paths = load_paths_WM(cfg.brainROI, cfg.net2load);
load([paths.results.DNNs f2sav '.mat']);

%load ([paths.results.clusters 'allClustInfo_PFC_BLNET']);
load ([paths.results.clusters 'all_clustinfo_VVS']);

clear nnHClust_pfc7 nnHClust_vvs4 nnHClust_vvs5 nnHClust_vvs6
for layi = 1 %1:size(nnFit{1}, 1)
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
            if strcmp(cfg.brainROI, 'vvs')
               nnHClust_vvs4(subji, :) = mean(nnHSubj(allClustInfo{4}.PixelIdxList{14}));
            end 
            if strcmp(cfg.brainROI, 'vvs')
               nnHClust_vvs5(subji, :) = mean(nnHSubj(allClustInfo{5}.PixelIdxList{25})); 
            end 
            if strcmp(cfg.brainROI, 'vvs')
               nnHClust_vvs6(subji, :) = mean(nnHSubj(allClustInfo{6}.PixelIdxList{17})); 
            end
            if strcmp(cfg.brainROI, 'pfc')
                nnHClust_pfc7(subji, :) = mean(nnHSubj(allClustInfo{7}.PixelIdxList{7})); 
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

%% Extract activity in specific clusters from NNH 

paths = load_paths_WM('vvs', 'none');
load ([paths.results.clusters 'clustinfo_PFC_px2']);

for subji = 1:size(nnH, 1)
    nnHSubj = squeeze(nnH(subji, :, :)); 
    nnHClust_vvsE1(subji, :) = mean(nnHSubj(clustinfo.PixelIdxList{2})); 
end

    
   



%% 
[h p ci ts] = ttest(nnHClust_vvsE1);

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

%% plot one bar
clear data
data.data = [nnHClust_vvsE1]; 

figure(2); set(gcf,'Position', [0 0 500 620]); 
mean_S = mean(data.data, 1);
hb = scatter([1], data.data, 100, 'k'); hold on;
%set(hb, 'lineWidth', 0.01, 'Marker', '.', 'MarkerSize',45);hold on;

h = bar (mean_S);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 3);
set(gca,'XTick',[1],'XTickLabel',{'', ''}, 'FontSize', 30, 'linew',2, 'xlim', [0 2] );
set(gca, 'ylim', [-.1 .15])
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





%% SAVE ONLY NEURAL RDMS 
clear, clc
%Network_ROI_EoM_layers_freqs_avRepet_avTimeFeatVect_freqResolv(0-1)__fitMode(0:noTrials; 1:Trials)__timeRes__win-width__mf
    
listF2sav = {

%'Alex_pfc_E123_[1-8]_3-54_1_0_1_0_.1_5_1'; 
%'Alex_vvs_E123_[1-8]_3-54_1_0_1_0_.1_5_1'; 
'Alex_pfc_M123_[1-8]_3-54_1_0_1_0_.1_5_1'; 
'Alex_vvs_M123_[1-8]_3-54_1_0_1_0_.1_5_1'; 

%'Alex_pfc_E123_[1-8]_3-54_0_0_1_0_.1_5_1'; 
%'Alex_vvs_E123_[1-8]_3-54_0_0_1_0_.1_5_1'; 
'Alex_pfc_M123_[1-8]_3-54_0_0_1_0_.1_5_1'; 
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
       
        [cfg_contrasts]             = getIdsWM(cfg.period, cfg_contrasts);
        cfg_contrasts.oneListPow    = extract_power_WM (cfg_contrasts, cfg); % 
        cfg_contrasts               = normalize_WM(cfg_contrasts, 1, 'sess', []);
        cfg_contrasts               = average_repetitions(cfg, cfg_contrasts);
        neuralRDMs                  = createNeuralRDMs(cfg, cfg_contrasts);
        allNeuralRDMS{sessi,1}      = neuralRDMs; 
        allNeuralRDMS{sessi,2}      = cfg_contrasts.oneListIds; 
        
    end
    
    save([paths.results.neuralRDMS.IT f2sav '.mat'], 'allNeuralRDMS');

end


%% %% SAVE ONLY NEURAL RDMS FROM POWER DATA
clear, clc
%Network_ROI_EoM_layers_freqs_avRepet_avTimeFeatVect_freqResolv(0-1)__fitMode(0:noTrials; 1:Trials)__timeRes__win-width__mf
    
listF2sav = {
'Alex_pfc_E123_[1-8]_3-54_1_0_1_0_.1_5_1'; 
'Alex_vvs_E123_[1-8]_3-54_1_0_1_0_.1_5_1'; 
'Alex_pfc_M123_[1-8]_3-54_1_0_1_0_.1_5_1'; 
'Alex_vvs_M123_[1-8]_3-54_1_0_1_0_.1_5_1'; 

'Alex_pfc_E123_[1-8]_3-54_0_0_1_0_.1_5_1'; 
'Alex_vvs_E123_[1-8]_3-54_0_0_1_0_.1_5_1'; 
'Alex_pfc_M123_[1-8]_3-54_0_0_1_0_.1_5_1'; 
'Alex_vvs_M123_[1-8]_3-54_0_0_1_0_.1_5_1'; 

};   

t1 = datetime; 
for listi = 1:length(listF2sav)
    disp(['File > ' num2str(listi) '      ' listF2sav{listi}]);
    clearvars -except listF2sav listi t1
        
    f2sav       = listF2sav{listi}; 
    cfg = getParams(f2sav);
    cfg.DNN_analysis = 1; 
    paths = load_paths_WM(cfg.brainROI, cfg.net2load);
    filelistSess = getFilesWM(paths.powerFromRT);
    
    for sessi= 1:length(filelistSess) 
        disp(['File > ' num2str(sessi)]);
        load([paths.powerFromRT filelistSess{sessi}]);   
        
        cfg_contrasts               = getIdsWM(cfg.period, cfg_contrasts);
    
        if length(cfg_contrasts.oneListIds) > 1
            cfg_contrasts               = average_repetitions(cfg, cfg_contrasts);
            neuralRDMs                  = createNeuralRDMs(cfg, cfg_contrasts);
            allNeuralRDMS{sessi,1}      = neuralRDMs; 
            allNeuralRDMS{sessi,2}      = cfg_contrasts.oneListIds; 
        end
    end
    
    save([paths.results.neuralRDMS f2sav(6:end) '.mat'], 'allNeuralRDMS');

end
t2 = datetime; 
etime(datevec(t2), datevec(t1))




%% Process and plot RDM in the PFC cluster during encoding
clear

f2load = 'pfc_E123_[1-8]_3-54_1_0_1_0_.1_5_1'; 
paths = load_paths_WM('vvs', 'none');
filelistSess = getFilesWM(paths.results.neuralRDMS);
load([paths.results.neuralRDMS f2load]);   

t1 = datetime; 
nSubjs = size(allNeuralRDMS, 1); 
nFreqs = size(allNeuralRDMS{1}, 3); 
nTimes = size(allNeuralRDMS{1}, 4); 

load ([paths.results.clusters 'CAT_PFC_ENC_px8']);
idsClust = clustinfo.PixelIdxList{8}; 

for subji = 1:nSubjs
    neuralRDMs  = allNeuralRDMS{subji, 1}; 
    ids         = allNeuralRDMS{subji, 2};
    nRows       = size(neuralRDMs, 1); 
    neuralRDM1 = reshape(neuralRDMs, nRows, nRows, nFreqs*nTimes); 
    rdm         = mean(neuralRDM1(:, :,idsClust ), 3); 
    meanRDM{subji,:} = rdm; 
    CM = load_CATMODEL_activ(ids); 
    rdm = vectorizeRDM(rdm); 
    CM = vectorizeRDM(CM); 
    allR(subji, :) = corr(CM', rdm, 'type', 's');    
    
    %z-score RDM
    rdm         = mean(neuralRDM1(:, :,idsClust ), 3); 
    CM = squeeze(load_M6_activ(ids));
    mW1 = mean(rdm(CM==1), 'omitnan'); mB1 = mean(rdm(CM ==-1), 'omitnan'); 
    mW2 = mean(rdm(CM==2), 'omitnan'); mB2 = mean(rdm(CM ==-2), 'omitnan'); 
    mW3 = mean(rdm(CM==3), 'omitnan'); mB3 = mean(rdm(CM ==-3), 'omitnan'); 
    mW4 = mean(rdm(CM==4), 'omitnan'); mB4 = mean(rdm(CM ==-4), 'omitnan'); 
    mW5 = mean(rdm(CM==5), 'omitnan'); mB5 = mean(rdm(CM ==-5), 'omitnan'); 
    mW6 = mean(rdm(CM==6), 'omitnan'); mB6 = mean(rdm(CM ==-6), 'omitnan');

    rdm(CM==1) = rdm(CM==1) - mW1; rdm(CM==-1) = rdm(CM==-1) - mB1; 
    rdm(CM==2) = rdm(CM==2) - mW2; rdm(CM==-2) = rdm(CM==-2) - mB2; 
    rdm(CM==3) = rdm(CM==3) - mW3; rdm(CM==-3) = rdm(CM==-3) - mB3; 
    rdm(CM==4) = rdm(CM==4) - mW4; rdm(CM==-4) = rdm(CM==-4) - mB4; 
    rdm(CM==5) = rdm(CM==5) - mW5; rdm(CM==-5) = rdm(CM==-5) - mB5; 
    rdm(CM==6) = rdm(CM==6) - mW6; rdm(CM==-6) = rdm(CM==-6) - mB6; 

    meanRDMZ{subji, :} = rdm;
    CM = load_CATMODEL_activ(ids); 
    rdmZ = vectorizeRDM(rdm); 
    CM = vectorizeRDM(CM); 
    allRZ(subji, :) = corr(CM', rdmZ, 'type', 's');    
end

[h p ci t] = ttest (allR);
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);


[h p ci t] = ttest (allRZ);
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);


%% plot for every subject 


for subji = 1:length(meanRDMZ)

    figure()
    tiledlayout(1,2)
    nexttile
    rdm = meanRDMZ{subji}; 
    rdm(logical(eye(size(rdm, 1)))) = nan; 
    imagesc(rdm); colorbar; axis square

    nexttile
    rdm = meanRDM{subji}; 
    rdm(logical(eye(size(rdm, 1)))) = nan; 
    imagesc(rdm); colorbar;axis square

end

%% plot mean
mRDM = cat(3, meanRDMZ{:}); 
mRDM = squeeze(mean(mRDM, 3)); 
mRDM(logical(eye(size(mRDM, 1)))) = nan; 
imagesc(mRDM); axis square; 

%% correlate with category model

for subji = 1:length(meanRDM)

    rdm = meanRDM{subji}; 
    CM = kron(eye(6), ones(10));
    rdm = vectorizeRDM(rdm); 
    CM = vectorizeRDM(CM); 
    allR(subji, :) = corr(CM, rdm, 'type', 's');    
end

[h p ci t] = ttest (allR);
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);



%% Process and plot RDM in the PFC cluster during maintenance
clear, clc

f2load = 'pfc_M123_[1-8]_3-54_1_0_1_0_.1_5_1'; 
paths = load_paths_WM('vvs', 'none');
filelistSess = getFilesWM(paths.results.neuralRDMS);
load([paths.results.neuralRDMS f2load]);   

t1 = datetime; 
nSubjs = size(allNeuralRDMS, 1); 
nFreqs = size(allNeuralRDMS{1}, 3); 
nTimes = size(allNeuralRDMS{1}, 4); 

load ([paths.results.clusters 'allClustInfo_PFC_BLNET']);
idsClust = allClustInfo{1,7}.PixelIdxList{7}; 

for subji = 1:nSubjs
    neuralRDMs  = allNeuralRDMS{subji, 1}; 
    ids         = allNeuralRDMS{subji, 2};
    nRows       = size(neuralRDMs, 1); 
    neuralRDM1 = reshape(neuralRDMs, nRows, nRows, nFreqs*nTimes); 
    rdm         = mean(neuralRDM1(:, :,idsClust ), 3); 
    meanRDM{subji,:} = rdm; 
    CM = load_CATMODEL_activ(ids); 
    rdm = vectorizeRDM(rdm); 
    CM = vectorizeRDM(CM); 
    allR(subji, :) = corr(CM', rdm, 'type', 's');    
    
    %z-score RDM
    rdm         = mean(neuralRDM1(:, :,idsClust ), 3); 
    CM = squeeze(load_M6_activ(ids));
    mW1 = mean(rdm(CM==1), 'omitnan'); mB1 = mean(rdm(CM ==-1), 'omitnan'); 
    mW2 = mean(rdm(CM==2), 'omitnan'); mB2 = mean(rdm(CM ==-2), 'omitnan'); 
    mW3 = mean(rdm(CM==3), 'omitnan'); mB3 = mean(rdm(CM ==-3), 'omitnan'); 
    mW4 = mean(rdm(CM==4), 'omitnan'); mB4 = mean(rdm(CM ==-4), 'omitnan'); 
    mW5 = mean(rdm(CM==5), 'omitnan'); mB5 = mean(rdm(CM ==-5), 'omitnan'); 
    mW6 = mean(rdm(CM==6), 'omitnan'); mB6 = mean(rdm(CM ==-6), 'omitnan');

    rdm(CM==1) = rdm(CM==1) - mW1; rdm(CM==-1) = rdm(CM==-1) - mB1; 
    rdm(CM==2) = rdm(CM==2) - mW2; rdm(CM==-2) = rdm(CM==-2) - mB2; 
    rdm(CM==3) = rdm(CM==3) - mW3; rdm(CM==-3) = rdm(CM==-3) - mB3; 
    rdm(CM==4) = rdm(CM==4) - mW4; rdm(CM==-4) = rdm(CM==-4) - mB4; 
    rdm(CM==5) = rdm(CM==5) - mW5; rdm(CM==-5) = rdm(CM==-5) - mB5; 
    rdm(CM==6) = rdm(CM==6) - mW6; rdm(CM==-6) = rdm(CM==-6) - mB6; 

    meanRDMZ{subji, :} = rdm;
    CM = load_CATMODEL_activ(ids); 
    rdmZ = vectorizeRDM(rdm); 
    CM = vectorizeRDM(CM); 
    allRZ(subji, :) = corr(CM', rdmZ, 'type', 's');    
end

[h p ci t] = ttest (allR);
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);


[h p ci t] = ttest (allRZ);
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

%% plot RDM in cluster for every subject 

for subji = 1:length(meanRDMZ)

    figure()
    tiledlayout(1,2)
    nexttile
    rdm = meanRDMZ{subji}; 
    rdm(logical(eye(size(rdm, 1)))) = nan; 
    imagesc(rdm); colorbar; axis square

    nexttile
    rdm = meanRDM{subji}; 
    rdm(logical(eye(size(rdm, 1)))) = nan; 
    imagesc(rdm); colorbar;axis square

end


%% plot mean
clear mRDM
meanRDMZ = meanRDMZ(2:end); 
mRDM = cat(3, meanRDMZ{:}); 
mRDM = squeeze(mean(mRDM, 3)); 
mRDM(logical(eye(size(mRDM, 1)))) = nan; 
imagesc(mRDM); axis square; 

 
%% correlate with category model

for subji = 1:length(meanRDM)


    rdm = meanRDM{subji}; 
    CM = kron(eye(6), ones(10));
    rdm = vectorizeRDM(rdm); 
    CM = vectorizeRDM(CM); 
    allR(subji, :) = corr(CM, rdm, 'type', 's');    
end

[h p ci t] = ttest (allR);
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);
%% Process and plot RDM during encoding
clear

f2load = 'vvs_E123_[1-8]_3-54_1_0_1_0_.1_5_1'; 
paths = load_paths_WM('vvs', 'none');
filelistSess = getFilesWM(paths.results.neuralRDMS);
load([paths.results.neuralRDMS f2load]);   

t1 = datetime; 
nSubjs = size(allNeuralRDMS, 1); 
nFreqs = size(allNeuralRDMS{1}, 3); 
nTimes = size(allNeuralRDMS{1}, 4); 

load ([paths.results.clusters 'allClustInfo_PFC_BLNET']);
idsClust = allClustInfo{1,7}.PixelIdxList{7}; 

for subji = 1:nSubjs
    neuralRDMs  = allNeuralRDMS{subji, 1}; 
    ids         = allNeuralRDMS{subji, 2};
    meanRDM{subji,:} = mean(mean(neuralRDMs(:, :,3:8, 8:12), 3), 4); 
end

%%test
% tfPlot = zeros(10, 10, 52, 46); 
% tfPlot1 = reshape(tfPlot, 10, 10, []); 
% tfPlot1(:, :, idsClust) = 1; 
% tfPlot1 = reshape(tfPlot1, 10, 10, 52, 46); 
% 
% imagesc(squeeze(tfPlot1(2,2,:,:)))


%% plot mean
meanRDM = meanRDM(2:end); 
mRDM = cat(3, meanRDM{:}); 
mRDM = squeeze(mean(mRDM, 3 )); 

mRDM(logical(eye(size(mRDM, 1)))) = nan; 
imagesc(mRDM)

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
for sessi= 6% 1:length(filelistSess) 
    
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




%% compute CCI in cluster PFC
clc
clearvars -except CCI
paths = load_paths_WM('vvs', 'none');
load ([paths.results.clusters 'allClustInfo_PFC_BLNET']);

clustinfo = allClustInfo{7}; 
for subji = 1:size(CCI, 1)
    cciSubj = squeeze(CCI(subji, :, 1:40)); 
    cciClust(subji, :) = mean(cciSubj(clustinfo.PixelIdxList{7}), 'all');
end

[h p ci ts] = ttest (cciClust)

%% plot one bar
clear data
data.data = [cciClust]; 

figure(2); set(gcf,'Position', [0 0 500 620]); 
mean_S = mean(data.data, 1);
hb = scatter([1], data.data, 100, 'k'); hold on;
%set(hb, 'lineWidth', 0.01, 'Marker', '.', 'MarkerSize',45);hold on;

h = bar (mean_S);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 3);
set(gca,'XTick',[1],'XTickLabel',{'', ''}, 'FontSize', 30, 'linew',2, 'xlim', [0 2] );
set(gca, 'ylim', [-.025 .05])
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 3);

[h p ci t] = ttest (data.data(:,1));
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth', 3);


exportgraphics(gcf, 'allM.png', 'Resolution', 300);

%% compute CCI in cluster VVS
clc
clearvars -except CCI
paths = load_paths_WM('vvs', 'none');
load ([paths.results.clusters 'all_clustinfo_VVS']);

for subji = 1:size(CCI, 1)
    cciSubj = squeeze(CCI(subji, :, :)); 
    
    clustinfo = allClustInfo{4};
    cciClust1(subji, :) = mean(cciSubj(clustinfo.PixelIdxList{14}), 'all');

    clustinfo = allClustInfo{5};
    cciClust2(subji, :) = mean(cciSubj(clustinfo.PixelIdxList{25}), 'all');

    clustinfo = allClustInfo{6};
    cciClust3(subji, :) = mean(cciSubj(clustinfo.PixelIdxList{17}), 'all');

end

[h p ci ts] = ttest (cciClust3)

%% plot one bar
clear data
data.data = [cciClust1]; 

figure(2); set(gcf,'Position', [0 0 500 620]); 
mean_S = mean(data.data, 1);
hb = scatter([1], data.data, 100, 'k'); hold on;
%set(hb, 'lineWidth', 0.01, 'Marker', '.', 'MarkerSize',45);hold on;

h = bar (mean_S);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 3);
set(gca,'XTick',[1],'XTickLabel',{'', ''}, 'FontSize', 30, 'linew',2, 'xlim', [0 2] );
set(gca, 'ylim', [-.025 .05])
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 3);

[h p ci t] = ttest (data.data(:,1));
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth', 3);


exportgraphics(gcf, 'allM.png', 'Resolution', 300);


%% COMPUTE PERMUTATIONS

clear

nPerm = 100; 
paths = load_paths_WM('vvs', 'none');
filelistSess = getFilesWM(paths.results.neuralRDMS.CAT);

for sessi= 3 %1:length(filelistSess) %this one starts at 1 and not at 3
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
for sessi= 1:2 %:length(filelistSess) %this one starts at 1 and not at 3
    
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

for sessi= 1%1:length(filelistSess) 
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
for sessi= 7 %:length(filelistSess) 
    
    disp(['File > ' num2str(sessi)]);

    clearvars -except paths filelistSess sessi
    load(filelistSess{sessi});   
    
    mICI = squeeze(mean(ICI, 'omitnan')); 
    [h p ci ts] = ttest(ICI);
    h = squeeze(h); t = squeeze(ts.tstat); 
    h(:, 1:5) = 0; % only sum p-values in clusters after the baseline
    clustinfo = bwconncomp(h);
    if ~isempty(clustinfo.PixelIdxList)
        for pixi = 1:length(clustinfo.PixelIdxList)
            allTObs(pixi, :) = sum(t(clustinfo.PixelIdxList{pixi}));
        end
    end
    [max2u id] = max(abs(allTObs));
    tObs = allTObs(id); 

    % % % ENCODING
    %h = zeros (52, 21); 
    %h(clustinfo.PixelIdxList{3}) = 1; 

    % % % MAINTENANCE
    % h = zeros (52, 46); 
    % h(clustinfo.PixelIdxList{10}) = 1; 

    

    if size(mICI, 2) == 21
        times = 1:210;
        freqs = 1:520; 
        figure; set(gcf, 'Position', [100 100 200 400])
        myCmap = colormap(brewermap([],'*Spectral'));
        colormap(myCmap)
        contourf(times, freqs, myresizem(mICI, 10), 100, 'linecolor', 'none'); hold on; colorbar
        contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 4);
        
        set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
        set(gca, 'xlim', [1 150], 'FontSize', 10); 
        set(gca, 'clim', [-.005 .005], 'FontSize', 10);
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
nPerm = 1000; 
paths = load_paths_WM('vvs', 'none');
filelistSess = getFilesWM(paths.results.neuralRDMS.IT);

t1 = datetime; 

for sessi= 1:length(filelistSess) 
    disp(['File > ' num2str(sessi)]);
    load([paths.results.neuralRDMS.IT filelistSess{sessi}]);   
    nSubjs = size(allNeuralRDMS, 1); 
    nFreqs = size(allNeuralRDMS{1}, 3); 
    nTimes = size(allNeuralRDMS{1}, 4); 

    clear ICIPerm
    for subji = 1:nSubjs
        disp(['Subji > ' num2str(subji)]);
        neuralRDMs  = allNeuralRDMS{subji, 1}; 
        ids         = allNeuralRDMS{subji, 2};
        
        parfor permi = 1:nPerm
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




%% process PERMUTATIONS


clear
nPerm = 1000; 
paths = load_paths_WM('vvs', 'none');
cd(paths.results.neuralRDMS.ITPlts);

filelistSess = dir('*mat'); filelistSess = [{filelistSess.name}'];
for sessi= 7 %:length(filelistSess) 
    
    disp(['File > ' num2str(sessi)]);

    clearvars -except paths filelistSess sessi nPerm
    load(filelistSess{sessi});   
    load(filelistSess{sessi+1});   

    % % % compute observed
    [h p ci ts] = ttest(ICI);
    h = squeeze(h); t = squeeze(ts.tstat); 
    h(:, 1:5) = 0; % only sum p-values in clusters after the baseline
    clustinfo = bwconncomp(h);
    if ~isempty(clustinfo.PixelIdxList)
        for pixi = 1:length(clustinfo.PixelIdxList)
            allTObs(pixi, :) = sum(t(clustinfo.PixelIdxList{pixi}));
        end
    end
    [max2u id] = max(abs(allTObs));
    tObs = allTObs(id); 


    for permi = 1:nPerm
        clear allTObsPerm
        ICIp = squeeze(ICIPerm(permi, :, :, :));
        [h p ci ts] = ttest(ICIp);
        h = squeeze(h); t = squeeze(ts.tstat); 
        h(:, 1:5) = 0; % only sum p-values in clusters after the baseline
        clustinfoP = bwconncomp(h);
        if ~isempty(clustinfoP.PixelIdxList)
            for pixi = 1:length(clustinfoP.PixelIdxList)
                allTObsPerm(pixi, :) = sum(t(clustinfoP.PixelIdxList{pixi}));
            end
        end
        [max2u id] = max(abs(allTObsPerm));
        tPerm(permi, :) = allTObsPerm(id); 

    end


end

allAB = tPerm(abs(tPerm) > abs(tObs));
p = 1 - ((nPerm-1) - (length (allAB)))  / nPerm




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












%% COMPUTE ONE ITEM LEVEL DNN FIT PER SUBJECT and CATEGORY, AVERAGE AND PERFORM STATS SIX CATEGORIES


clear, clc
%Network_ROI_EoM_layers_freqs_avRepet_avTimeFeatVect_freqResolv(0-1)__fitMode(0:noTrials; 1:Trials)__timeRes__win-width__mf
    
lay2u = 1; 

listF2sav = {

% % 
% 'BLNETi_vvs_E123CAT1_[8-8-56]_3-54_0_0_1_0_.1_5_1';
% 'BLNETi_vvs_E123CAT2_[8-8-56]_3-54_0_0_1_0_.1_5_1';
% 'BLNETi_vvs_E123CAT3_[8-8-56]_3-54_0_0_1_0_.1_5_1';
% 'BLNETi_vvs_E123CAT4_[8-8-56]_3-54_0_0_1_0_.1_5_1';
% 'BLNETi_vvs_E123CAT5_[8-8-56]_3-54_0_0_1_0_.1_5_1';
% 'BLNETi_vvs_E123CAT6_[8-8-56]_3-54_0_0_1_0_.1_5_1';

% 'BLNETi_pfc_E123CAT1_[8-8-56]_3-54_0_0_1_0_.1_5_1';
% 'BLNETi_pfc_E123CAT2_[8-8-56]_3-54_0_0_1_0_.1_5_1';
% 'BLNETi_pfc_E123CAT3_[8-8-56]_3-54_0_0_1_0_.1_5_1';
% 'BLNETi_pfc_E123CAT4_[8-8-56]_3-54_0_0_1_0_.1_5_1';
% 'BLNETi_pfc_E123CAT5_[8-8-56]_3-54_0_0_1_0_.1_5_1';
% 'BLNETi_pfc_E123CAT6_[8-8-56]_3-54_0_0_1_0_.1_5_1';

% 'BLNETi_vvs_M123CAT1_[8-8-56]_3-54_0_0_1_0_.1_5_1';
% 'BLNETi_vvs_M123CAT2_[8-8-56]_3-54_0_0_1_0_.1_5_1';
% 'BLNETi_vvs_M123CAT3_[8-8-56]_3-54_0_0_1_0_.1_5_1';
% 'BLNETi_vvs_M123CAT4_[8-8-56]_3-54_0_0_1_0_.1_5_1';
% 'BLNETi_vvs_M123CAT5_[8-8-56]_3-54_0_0_1_0_.1_5_1';
% 'BLNETi_vvs_M123CAT6_[8-8-56]_3-54_0_0_1_0_.1_5_1';

% 'BLNETi_pfc_M123CAT1_[8-8-56]_3-54_0_0_1_0_.1_5_1';
% 'BLNETi_pfc_M123CAT2_[8-8-56]_3-54_0_0_1_0_.1_5_1';
% 'BLNETi_pfc_M123CAT3_[8-8-56]_3-54_0_0_1_0_.1_5_1';
% 'BLNETi_pfc_M123CAT4_[8-8-56]_3-54_0_0_1_0_.1_5_1';
% 'BLNETi_pfc_M123CAT5_[8-8-56]_3-54_0_0_1_0_.1_5_1';
% 'BLNETi_pfc_M123CAT6_[8-8-56]_3-54_0_0_1_0_.1_5_1';

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


% 'Alex_vvs_M123CAT1_[1-8]_3-54_0_0_1_0_.1_5_1';
% 'Alex_vvs_M123CAT2_[1-8]_3-54_0_0_1_0_.1_5_1';
% 'Alex_vvs_M123CAT3_[1-8]_3-54_0_0_1_0_.1_5_1';
% 'Alex_vvs_M123CAT4_[1-8]_3-54_0_0_1_0_.1_5_1';
% 'Alex_vvs_M123CAT5_[1-8]_3-54_0_0_1_0_.1_5_1';
% 'Alex_vvs_M123CAT6_[1-8]_3-54_0_0_1_0_.1_5_1';

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
    clearvars -except listF2sav listi allnnFIT lay2u
        
    f2sav       = listF2sav{listi}; 
    cfg = getParams(f2sav);
    paths = load_paths_WM(cfg.brainROI, cfg.net2load);
    load([paths.results.DNNs f2sav '.mat']);
    x =  nnFit(:, 1);
    x = cellfun(@(x) x(lay2u,:,:), x, 'un', 0); 
    allnnFIT{listi,:} = x; 

end


allNF = cat(2, allnnFIT{:}); 
allNF = cellfun(@(x) squeeze(x), allNF, 'un', 0);
allSNF = cell2mat(permute(allNF, [3 4 1 2])); 
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


if strcmp(cfg.period(1), 'M')
    set(gcf, 'Position', [100 100 1200 1000])
    set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
    set(gca, 'xlim', [1 40], 'FontSize', 10); %
    plot([5 5],get(gca,'ylim'), 'k:','lineWidth', 2);
else
    set(gcf, 'Position', [100 100 150 300])
    set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
    set(gca, 'FontSize', 8, 'clim', [-5 5]);
    plot([5 5],get(gca,'ylim'), 'k:','lineWidth', 2);
    
end

exportgraphics(gcf, ['myP.png'], 'Resolution', 300); 

%% Extract activity in specific clusters from NNH


paths = load_paths_WM('vvs', 'none');

load ([paths.results.clusters 'ITM_VVS_enc']);
%load ([paths.results.clusters 'clust_PFC_ITM_ENC_px3']);


for subji = 1:size(nnH, 1)
    nnHSubj = squeeze(nnH(subji, :, :)); 
    nnHClust_vvsE1(subji, :) = mean(nnHSubj(clustinfo.PixelIdxList{6})); 
    %nnHClust_vvsE1(subji, :) = mean(nnHSubj(clustinfo.PixelIdxList{3})); 
end

% % % test cluster is correct
%tes2u = zeros(52, 15); 
%tes2u(clustinfo.PixelIdxList{6})=1; 
%figure()
%imagesc(tes2u);

[h p ci ts] = ttest(nnHClust_vvsE1);

h = squeeze(h); t = squeeze(ts.tstat); 

disp (['t = ' num2str(ts.tstat) '  ' ' p = ' num2str(p)]);

%% plot one bar
clear data
data.data = [nnHClust_vvsE1]; 

figure(2); set(gcf,'Position', [0 0 500 620]); 
mean_S = mean(data.data, 1);
hb = scatter([1], data.data, 100, 'k'); hold on;
%set(hb, 'lineWidth', 0.01, 'Marker', '.', 'MarkerSize',45);hold on;

h = bar (mean_S);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 3);
set(gca,'XTick',[1],'XTickLabel',{'', ''}, 'FontSize', 30, 'linew',2, 'xlim', [0 2] );
set(gca, 'ylim', [-.1 .15])
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 3);

[h p ci t] = ttest (data.data(:,1));
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth', 3);


exportgraphics(gcf, 'allM.png', 'Resolution', 300);

%% plot all nnH

for subji = 1:size(nnH, 1)
    
    figure()
    imagesc(squeeze(nnH(subji, :, :))); colorbar

end







%% Ttest at every time-frequency point comparing VVS and PFC fits 

%Network_ROI_ER_layers_freqs_avRepet_avTFV_fRes(0-1)_fitMode(0:noTrials; 1:Trials)_timeRes_win_mf
clear , clc

f2sav1 = 'Alex_vvs_E123_[1-8]_3-54_1_0_1_0_.1_5_1';
%f2sav1 = 'BLNETi_vvs_E123_[8-8-56]_3-54_1_0_1_0_.1_5_1';
%f2sav1 = 'CORrt_vvs_E123_[2-2-8]_3-54_1_0_1_0_.1_5_1';
cfg = getParams(f2sav1);
paths = load_paths_WM(cfg.brainROI, cfg.net2load);
load([paths.results.DNNs f2sav1 '.mat']);
nnFitVVS = nnFit; 

f2sav2 = 'Alex_pfc_E123_[1-8]_3-54_1_0_1_0_.1_5_1';
%f2sav2 = 'BLNETi_pfc_E123_[8-8-56]_3-54_1_0_1_0_.1_5_1';
%f2sav2 = 'CORrt_pfc_E123_[2-2-8]_3-54_1_0_1_0_.1_5_1';
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

%% Ttest at every time-frequency point comparing fits in ALEXNET VS BL-NET

%Network_ROI_ER_layers_freqs_avRepet_avTFV_fRes(0-1)_fitMode(0:noTrials; 1:Trials)_timeRes_win_mf
clear , clc

f2sav1 = 'BLNETi_pfc_M123_[8-8-56]_3-54_1_0_1_0_.1_5_1';
cfg = getParams(f2sav1);
paths = load_paths_WM(cfg.brainROI, cfg.net2load);
load([paths.results.DNNs f2sav1 '.mat']);
nnFitVVS = nnFit; 

f2sav2 = 'Alex_pfc_M123_[1-8]_3-54_1_0_1_0_.1_5_1';
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

for layi = 1:7 %:size(nnFit{1}, 1)
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
    
    if strcmp(cfg.brainROI, 'vvs')
        nnHVVS([18 22], :, :) = []; 
        nnHPFC([18 22], :, :) = []; 
    else
        nnHVVS([1], :, :) = []; 
        nnHPFC([1], :, :) = []; 
    end


    nnHVVS = squeeze(nnHVVS);
    nnHPFC = squeeze(nnHPFC);

    [h p ci ts] = ttest(nnHVVS, nnHPFC);
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

    if exist('allTObs1')
        [max2u id] = max(abs(allTObs1));
        tObs = allTObs1(id); 
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
        set(gca, 'xlim', [1 40], 'clim', [-3 3], 'FontSize', 10);
        plot([5 5],get(gca,'ylim'), 'k:','lineWidth', 2);
    else
        set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
        set(gca, 'FontSize', 8, 'clim', [-5 5]);
        plot([5 5],get(gca,'ylim'), 'k:','lineWidth', 2);
        
    end
    

end

%exportgraphics(gcf, [paths.results.DNNs 'myP.png'], 'Resolution', 300); 
exportgraphics(gcf, ['myP.png'], 'Resolution', 300); 

%% Ttest at every time-frequency point comparing fits in ALEXNET VS BL-NET (LAYER 7 vs 8)

%Network_ROI_ER_layers_freqs_avRepet_avTFV_fRes(0-1)_fitMode(0:noTrials; 1:Trials)_timeRes_win_mf
clear , clc

f2sav1 = 'BLNETi_pfc_M123_[8-8-56]_3-54_1_0_1_0_.1_5_1';
cfg = getParams(f2sav1);
paths = load_paths_WM(cfg.brainROI, cfg.net2load);
load([paths.results.DNNs f2sav1 '.mat']);
nnFitBLNET = nnFit; 

f2sav2 = 'Alex_pfc_M123_[1-8]_3-54_1_0_1_0_.1_5_1';
cfg = getParams(f2sav2);
paths = load_paths_WM(cfg.brainROI, cfg.net2load);
load([paths.results.DNNs f2sav2 '.mat']);
nnFitALEX = nnFit; 


tiledlayout(8,9, 'TileSpacing', 'tight', 'Padding', 'compact');
if strcmp(cfg.period(1), 'M')
    set(gcf, 'Position', [100 100 1200 1000])
else
    set(gcf, 'Position', [100 100 700 1200])
end


ax1 = nexttile;
clear nnHVVS nnHPFC
for subji = 1:length(nnFitBLNET)
   if strcmp(cfg.period(1), 'M')
     nnHBLNET(subji, : ,:) = atanh(nnFitBLNET{subji, 1}(4,:,1:40));
   else
     nnHBLNET(subji, : ,:) = atanh(nnFitBLNET{subji, 1}(4,:,1:15));
   end
end
for subji = 1:length(nnFitALEX)
   if strcmp(cfg.period(1), 'M')
     nnHALEX(subji, : ,:) = atanh(nnFitALEX{subji, 1}(4,:,1:40));
   else
     nnHALEX(subji, : ,:) = atanh(nnFitALEX{subji, 1}(4,:,1:15));
   end
end

if strcmp(cfg.brainROI, 'vvs')
    nnHBLNET([18 22], :, :) = []; 
    nnHALEX([18 22], :, :) = []; 
else
    nnHBLNET([1], :, :) = []; 
    nnHALEX([1], :, :) = []; 
end


nnHBLNET = squeeze(nnHBLNET);
nnHALEX = squeeze(nnHALEX);

[h p ci ts] = ttest(nnHBLNET, nnHALEX);
h = squeeze(h); t = squeeze(ts.tstat); 
h(:, 1:5) = 0; % only sum p-values in clusters after the baseline

freqs = 1:52; 
clustinfo = bwconncomp(h);
allClustInfo{1} = clustinfo; 

% store allTOBS
if ~isempty(clustinfo.PixelIdxList)
    for pixi = 1:length(clustinfo.PixelIdxList)
         allTObs(1, pixi, :) = sum(t(clustinfo.PixelIdxList{pixi}));
         allTObs1(pixi, :) = sum(t(clustinfo.PixelIdxList{pixi}));
    end
else
    allTObs(1, :, :) = 0;
end

if exist('allTObs1')
    [max2u id] = max(abs(allTObs1));
    tObs = allTObs1(id); 
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
    set(gca, 'xlim', [1 40], 'clim', [-3 3], 'FontSize', 10);
    plot([5 5],get(gca,'ylim'), 'k:','lineWidth', 2);
else
    set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
    set(gca, 'FontSize', 8, 'clim', [-5 5]);
    plot([5 5],get(gca,'ylim'), 'k:','lineWidth', 2);
    
end


%exportgraphics(gcf, [paths.results.DNNs 'myP.png'], 'Resolution', 300); 
exportgraphics(gcf, ['myP.png'], 'Resolution', 300); 


%% Permutations 
nPerm = 1000; 


clear max_clust_sum_perm allTObs
for permi = 1:nPerm 

    junts = cat(1, nnHBLNET, nnHALEX);
    junts = junts (randperm(size(junts, 1)), :, :);
    nnHBLNETPerm = junts(1:floor(size(junts, 1)/2), :, :); 
    nnHALEXPerm = junts(floor(size(junts, 1)/2)+1:end, :, :); 
    
    [h p ci ts] = ttest(nnHBLNETPerm, nnHALEXPerm);
    h = squeeze(h); t = squeeze(ts.tstat); 
    h(:, 1:5) = 0; % only sum p-values in clusters after the baseline

    clustinfo = bwconncomp(h);
    if ~isempty(clustinfo.PixelIdxList)
        for pixi = 1:length(clustinfo.PixelIdxList)
             allTObs(pixi, :) = sum(t(clustinfo.PixelIdxList{pixi}));
        end
    else
        allTObs(permi,:) = 0;
    end
    
    if exist('allTObs')
        [max2u id] = max(abs(allTObs));
        max_clust_sum_perm(permi,:) = allTObs(id); 
    end


end

allAB = max_clust_sum_perm(max_clust_sum_perm> tObs);
p = 1 - ((nPerm-1) - (length (allAB)))  / nPerm;

disp(['p = ' num2str(p)])



%% Extract activity in specific clusters DURING MAINTENANCE (FOR PERFORMANCE; OR BL-NET FITS)

clear nnHClust_pfc7 nnHClust_vvs4 nnHClust_vvs5 nnHClust_vvs6
load ([paths.results.clusters 'allClustInfo_PFC_BLNET']);
for subji = 1:size(nnHALEX, 1)
    nnHSubj = squeeze(nnHBLNET(subji, :, :)) - squeeze(nnHALEX(subji, :, :)) ; 
    nnHClust_pfc7(subji, :) = mean(nnHSubj(allClustInfo{7}.PixelIdxList{7})); 
end

%% plot one bar
clear data
data.data = [nnHClust_pfc7]; 

figure(2); set(gcf,'Position', [0 0 500 620]); 
mean_S = mean(data.data, 1);
hb = scatter([1], data.data, 100, 'k'); hold on;
%set(hb, 'lineWidth', 0.01, 'Marker', '.', 'MarkerSize',45);hold on;

h = bar (mean_S);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 3);
set(gca,'XTick',[1],'XTickLabel',{'', ''}, 'FontSize', 30, 'linew',2, 'xlim', [0 2] );
set(gca, 'ylim', [-.1 .15])
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 3);

[h p ci t] = ttest (data.data(:,1));
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth', 3);


exportgraphics(gcf, 'allM.png', 'Resolution', 300);



%% Ttest at every time-frequency point comparing fits in BL-NET (LAYER 7) VS CATEGORY MODEL 

%Network_ROI_ER_layers_freqs_avRepet_avTFV_fRes(0-1)_fitMode(0:noTrials; 1:Trials)_timeRes_win_mf
clear , clc

f2sav1 = 'BLNETi_vvs_M123_[8-8-56]_3-54_1_0_1_0_.1_5_1';
cfg = getParams(f2sav1);
paths = load_paths_WM(cfg.brainROI, cfg.net2load);
load([paths.results.DNNs f2sav1 '.mat']);
nnFitBLNET = nnFit; 

f2sav2 = 'CAT_vvs_M123_[1]_3-54_1_0_1_0_.1_5_1';
cfg = getParams(f2sav2);
paths = load_paths_WM(cfg.brainROI, cfg.net2load);
load([paths.results.DNNs f2sav2 '.mat']);
nnFitCAT = nnFit; 


tiledlayout(8,9, 'TileSpacing', 'tight', 'Padding', 'compact');
if strcmp(cfg.period(1), 'M')
    set(gcf, 'Position', [100 100 1200 1000])
else
    set(gcf, 'Position', [100 100 700 1200])
end


ax1 = nexttile;
clear nnHVVS nnHPFC
for subji = 1:length(nnFitBLNET)
   if strcmp(cfg.period(1), 'M')
     nnHBLNET(subji, : ,:) = atanh(nnFitBLNET{subji, 1}(4,:,1:40));
   else
     nnHBLNET(subji, : ,:) = atanh(nnFitBLNET{subji, 1}(4,:,1:15));
   end
end
for subji = 1:length(nnFitCAT)
   if strcmp(cfg.period(1), 'M')
     nnHCAT(subji, : ,:) = atanh(nnFitCAT{subji, 1}(1,:,1:40));
   else
     nnHCAT(subji, : ,:) = atanh(nnFitCAT{subji, 1}(1,:,1:15));
   end
end

if strcmp(cfg.brainROI, 'vvs')
    nnHBLNET([18 22], :, :) = []; 
    nnHCAT([18 22], :, :) = []; 
else
    nnHBLNET([1], :, :) = []; 
    nnHCAT([1], :, :) = []; 
end


nnHBLNET = squeeze(nnHBLNET);
nnHCAT = squeeze(nnHCAT);

[h p ci ts] = ttest(nnHBLNET, nnHCAT);
h = squeeze(h); t = squeeze(ts.tstat); 
h(:, 1:5) = 0; % only sum p-values in clusters after the baseline

freqs = 1:52; 
clustinfo = bwconncomp(h);
allClustInfo{1} = clustinfo; 

% store allTOBS
if ~isempty(clustinfo.PixelIdxList)
    for pixi = 1:length(clustinfo.PixelIdxList)
         allTObs(1, pixi, :) = sum(t(clustinfo.PixelIdxList{pixi}));
         allTObs1(pixi, :) = sum(t(clustinfo.PixelIdxList{pixi}));
    end
else
    allTObs(1, :, :) = 0;
end

if exist('allTObs1')
    [max2u id] = max(abs(allTObs1));
    tObs = allTObs1(id); 
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
    set(gca, 'xlim', [1 40], 'clim', [-3 3], 'FontSize', 10);
    plot([5 5],get(gca,'ylim'), 'k:','lineWidth', 2);
else
    set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
    set(gca, 'FontSize', 8, 'clim', [-5 5]);
    plot([5 5],get(gca,'ylim'), 'k:','lineWidth', 2);
    
end


%exportgraphics(gcf, [paths.results.DNNs 'myP.png'], 'Resolution', 300); 
exportgraphics(gcf, ['myP.png'], 'Resolution', 300); 


%% Permutations 
nPerm = 1000; 


clear max_clust_sum_perm allTObs
for permi = 1:nPerm 

    junts = cat(1, nnHBLNET, nnHCAT);
    junts = junts (randperm(size(junts, 1)), :, :);
    nnHBLNETPerm = junts(1:floor(size(junts, 1)/2), :, :); 
    nnHALEXPerm = junts(floor(size(junts, 1)/2)+1:end, :, :); 
    
    [h p ci ts] = ttest(nnHBLNETPerm, nnHALEXPerm);
    h = squeeze(h); t = squeeze(ts.tstat); 
    h(:, 1:5) = 0; % only sum p-values in clusters after the baseline

    clustinfo = bwconncomp(h);
    if ~isempty(clustinfo.PixelIdxList)
        for pixi = 1:length(clustinfo.PixelIdxList)
             allTObs(pixi, :) = sum(t(clustinfo.PixelIdxList{pixi}));
        end
    else
        allTObs(permi,:) = 0;
    end
    
    if exist('allTObs')
        [max2u id] = max(abs(allTObs));
        max_clust_sum_perm(permi,:) = allTObs(id); 
    end


end

allAB = max_clust_sum_perm(max_clust_sum_perm> tObs);
p = 1 - ((nPerm-1) - (length (allAB)))  / nPerm;

disp(['p = ' num2str(p)])



%% CORRELATION OF TF MAPS in BL-NET (LAYER 7) VS CATEGORY MODEL 

%Network_ROI_ER_layers_freqs_avRepet_avTFV_fRes(0-1)_fitMode(0:noTrials; 1:Trials)_timeRes_win_mf
clear , clc

f2sav1 = 'BLNETi_vvs_M123_[8-8-56]_3-54_1_0_1_0_.1_5_1';
cfg = getParams(f2sav1);
paths = load_paths_WM(cfg.brainROI, cfg.net2load);
load([paths.results.DNNs f2sav1 '.mat']);
nnFitBLNET = nnFit; 

f2sav2 = 'CAT_vvs_M123_[1]_3-54_1_0_1_0_.1_5_1';
cfg = getParams(f2sav2);
paths = load_paths_WM(cfg.brainROI, cfg.net2load);
load([paths.results.DNNs f2sav2 '.mat']);
nnFitCAT = nnFit; 

clear nnHVVS nnHPFC
for subji = 1:length(nnFitBLNET)
   if strcmp(cfg.period(1), 'M')
     nnHBLNET(subji, : ,:) = atanh(nnFitBLNET{subji, 1}(4,:,1:40));
   else
     nnHBLNET(subji, : ,:) = atanh(nnFitBLNET{subji, 1}(4,:,1:15));
   end
end
for subji = 1:length(nnFitCAT)
   if strcmp(cfg.period(1), 'M')
     nnHCAT(subji, : ,:) = atanh(nnFitCAT{subji, 1}(1,:,1:40));
   else
     nnHCAT(subji, : ,:) = atanh(nnFitCAT{subji, 1}(1,:,1:15));
   end
end

if strcmp(cfg.brainROI, 'vvs')
    nnHBLNET([18 22], :, :) = []; 
    nnHCAT([18 22], :, :) = []; 
else
    nnHBLNET([1], :, :) = []; 
    nnHCAT([1], :, :) = []; 
end


clear allR
for subji = 1:size(nnHBLNET, 1)
    
    nnHBLNETS = nnHBLNET(subji, :, :); 
    nnHCATS = nnHCAT(subji, :, :); 

    nnHBLNETS = nnHBLNETS(:); 
    nnHCATS = nnHCATS(:); 

    allR(subji, :) = corr(nnHBLNETS, nnHCATS, 'type', 's'); 

end


%% plot one bar
clear data
data.data = [allR]; 

figure(2); set(gcf,'Position', [0 0 500 620]); 
mean_S = mean(data.data, 1);
hb = scatter([1], data.data, 100, 'k'); hold on;
%set(hb, 'lineWidth', 0.01, 'Marker', '.', 'MarkerSize',45);hold on;

h = bar (mean_S);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 3);
set(gca,'XTick',[1],'XTickLabel',{'', ''}, 'FontSize', 30, 'linew',2, 'xlim', [0 2] );
set(gca, 'ylim', [0 1])
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 3);

[h p ci t] = ttest (data.data(:,1));
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth', 3);


exportgraphics(gcf, 'allM.png', 'Resolution', 300);



%% Ttest at every time-frequency point comparing correct vs incorrect in VVS 

%Network_ROI_ER_layers_freqs_avRepet_avTFV_fRes(0-1)_fitMode(0:noTrials; 1:Trials)_timeRes_win_mf
clear , clc

f2sav1 = 'Alex_vvs_E123CC_[1-8]_3-54_0_0_1_0_.1_5_1';
%f2sav1 = 'Alex_vvs_E123_[1-8]_3-54_1_0_1_0_.1_5_1';
%f2sav1 = 'BLNETi_vvs_E123_[8-8-56]_3-54_1_0_1_0_.1_5_1';
%f2sav1 = 'CORrt_vvs_E123_[2-2-8]_3-54_1_0_1_0_.1_5_1';
cfg = getParams(f2sav1);
paths = load_paths_WM(cfg.brainROI, cfg.net2load);
load([paths.results.DNNs f2sav1 '.mat']);
nnFitVVS1 = nnFit; 

f2sav2 = 'Alex_vvs_E123IC_[1-8]_3-54_0_0_1_0_.1_5_1';
%f2sav2 = 'Alex_pfc_E123_[1-8]_3-54_1_0_1_0_.1_5_1';
%f2sav2 = 'BLNETi_pfc_E123_[8-8-56]_3-54_1_0_1_0_.1_5_1';
%f2sav2 = 'CORrt_pfc_E123_[2-2-8]_3-54_1_0_1_0_.1_5_1';
cfg = getParams(f2sav2);
paths = load_paths_WM(cfg.brainROI, cfg.net2load);
load([paths.results.DNNs f2sav2 '.mat']);
nnFitVVS2 = nnFit; 

tiledlayout(8,9, 'TileSpacing', 'tight', 'Padding', 'compact');
if strcmp(cfg.period(1), 'M')
    set(gcf, 'Position', [100 100 1200 1000])
else
    set(gcf, 'Position', [100 100 700 1200])
end

for layi = 1:size(nnFit{1}, 1)
    ax1 = nexttile;
    clear nnHVVS nnHPFC
    for subji = 1:length(nnFitVVS1)
       if strcmp(cfg.period(1), 'M')
         nnHVVS1(subji, : ,:) = atanh(nnFitVVS1{subji, 1}(layi,:,1:40));
       else
         nnHVVS1(subji, : ,:) = atanh(nnFitVVS1{subji, 1}(layi,:,1:15));
       end
    end
    for subji = 1:length(nnFitVVS2)
       if strcmp(cfg.period(1), 'M')
         nnHVVS2(subji, : ,:) = atanh(nnFitVVS2{subji, 1}(layi,:,1:40));
       else
          if ~isempty(nnFitVVS2{subji, 1})
                nnHVVS2(subji, : ,:) = atanh(nnFitVVS2{subji, 1}(layi,:,1:15));
          else
                nnHVVS2(subji, : ,:) = nan; 
          end
       end
    end
    
    
    nnHVVS1 = squeeze(nnHVVS1);
    nnHVVS2 = squeeze(nnHVVS2);

    %sub2exc = [18 22 10 20]; % for the incorrect trials
    sub2exc = [11 18 22 ];
    nnHVVS1(sub2exc,:,:) = []; 
    nnHVVS2(sub2exc,:,:) = []; 

    [h p ci ts] = ttest2(nnHVVS1, nnHVVS2);
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

%% Ttest at every time-frequency point comparing correct vs incorrect in PFC 

%Network_ROI_ER_layers_freqs_avRepet_avTFV_fRes(0-1)_fitMode(0:noTrials; 1:Trials)_timeRes_win_mf
clear , clc

f2sav1 = 'Alex_pfc_E123CC_[1-8]_3-54_0_0_1_0_.1_5_1';
cfg = getParams(f2sav1);
paths = load_paths_WM(cfg.brainROI, cfg.net2load);
load([paths.results.DNNs f2sav1 '.mat']);
nnFitVVS1 = nnFit; 

f2sav2 = 'Alex_pfc_E123IC_[1-8]_3-54_0_0_1_0_.1_5_1';
cfg = getParams(f2sav2);
paths = load_paths_WM(cfg.brainROI, cfg.net2load);
load([paths.results.DNNs f2sav2 '.mat']);
nnFitVVS2 = nnFit; 

tiledlayout(8,9, 'TileSpacing', 'tight', 'Padding', 'compact');
if strcmp(cfg.period(1), 'M')
    set(gcf, 'Position', [100 100 1200 1000])
else
    set(gcf, 'Position', [100 100 700 1200])
end

for layi = 1:size(nnFit{1}, 1)
    ax1 = nexttile;
    clear nnHVVS nnHPFC
    for subji = 1:length(nnFitVVS1)
       if strcmp(cfg.period(1), 'M')
         nnHVVS1(subji, : ,:) = atanh(nnFitVVS1{subji, 1}(layi,:,1:40));
       else
         nnHVVS1(subji, : ,:) = atanh(nnFitVVS1{subji, 1}(layi,:,1:15));
       end
    end
    for subji = 1:length(nnFitVVS2)
       if strcmp(cfg.period(1), 'M')
         nnHVVS2(subji, : ,:) = atanh(nnFitVVS2{subji, 1}(layi,:,1:40));
       else
          if ~isempty(nnFitVVS2{subji, 1})
                nnHVVS2(subji, : ,:) = atanh(nnFitVVS2{subji, 1}(layi,:,1:15));
          else
                nnHVVS2(subji, : ,:) = nan; 
          end
       end
    end
    
    
    nnHVVS1 = squeeze(nnHVVS1);
    nnHVVS2 = squeeze(nnHVVS2);

    nnHVVS1([11 18 22],:,:) = []; 
    nnHVVS2([11 18 22],:,:) = []; 

    [h p ci ts] = ttest2(nnHVVS1, nnHVVS2);
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

