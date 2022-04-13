%% Calculate from epoched raw traces
clear 
tic
[ allTraces allIDs allChans] = loadTraces ('PFC');
toc


%%
%ROI_layers_freqs_freqResolv(0-1)_fitMode(0:noTrials; 1:Trials)_timeRes_win-width_mf
f2sav = 'PFC_8-16-24-32-40-48-56_13-29_0_1_all_1_1'; 
cfg = getParams(f2sav);

t1 = datetime; 

for subji = 1:length(allTraces)
   
    oneListTraces           = allTraces{subji};
    oneListIDs              = allIDs{subji};
    [oneListPow]            = extract_power_WM (oneListTraces, cfg.timeRes); % 
    oneListPow              = norm_trials(oneListPow);
    idsh                    = cell2mat(cellfun(@(x) ~strcmp(x(6), '4'), oneListIDs, 'un', 0));
    oneListIDs              = oneListIDs(idsh);
    oneListPow              = oneListPow(idsh,:,:,:);
    neuralRDMs              = createNeuralRDMs(oneListPow, cfg.freqs, cfg.win_width, cfg.mf, cfg.fR);
    networkRDMs             = createNetworkRDMs(oneListIDs, cfg.lays2load, cfg.brainROI, subji); %layers to load
    
    nnFit{subji,:}          = fitModel_WM(neuralRDMs, networkRDMs, cfg.fitMode); 
    
end

save (f2sav, 'nnFit');
t2 = datetime; 
etime(datevec(t2), datevec(t1))



%%
%ROI_layers_freqs_freqResolv(0-1)_fitMode(0:noTrials; 1:Trials)_timeRes_win-width_mf
clear 

files2process = {'PFC_8-16-24-32-40-48-56_01-54_1_1_0.1_5_1_.mat', ... 
                 'VVS_8-16-24-32-40-48-56_01-54_1_1_0.1_5_1_.mat', ... 
                 'HIP_8-16-24-32-40-48-56_01-54_1_1_0.1_5_1_.mat', ... 
                 }



for filei = 1:length(files2process)             
    
    f2sav = files2process{filei}
    f2t = strsplit(f2sav, '_'); 
    
    [ allTraces allIDs allChans] = loadTraces (f2t{1});
    
    cfg = getParams(f2sav);

    t1 = datetime; 

    for subji = 1:length(allTraces)

        oneListTraces           = allTraces{subji};
        oneListIDs              = allIDs{subji};
        [oneListPow]            = extract_power_WM (oneListTraces, cfg.timeRes); % 
        oneListPow              = norm_trials(oneListPow);
        idsh                    = cell2mat(cellfun(@(x) ~strcmp(x(6), '4'), oneListIDs, 'un', 0));
        oneListIDs              = oneListIDs(idsh);
        oneListPow              = oneListPow(idsh,:,:,:);
        neuralRDMs              = createNeuralRDMs(oneListPow, cfg.freqs, cfg.win_width, cfg.mf, cfg.fR);
        networkRDMs             = createNetworkRDMs(oneListIDs, cfg.lays2load, cfg.brainROI, subji); %layers to load

        nnFit{subji,:}          = fitModel_WM(neuralRDMs, networkRDMs, cfg.fitMode); 

    end

    save (f2sav, 'nnFit');
    
end 
t2 = datetime; 
etime(datevec(t2), datevec(t1))


%% calculate only phases

f_PHA = [3:8];
for subji = 1:length(allTraces)
   
    oneListTraces           = allTraces{subji};
    oneListIDs              = allIDs{subji};
    idsh                    = cell2mat(cellfun(@(x) ~strcmp(x(6), '4'), oneListIDs, 'un', 0));
    oneListIDs              = oneListIDs(idsh);
    oneListTraces           = oneListTraces(:, :, idsh);
    oneListPhases           = extract_phases_WM(oneListTraces, f_PHA); 
    allPhases{subji}        = oneListPhases; 
    
end

%% build distributions for all PFC electrodes
clearvars -except allPhases nnFit

f_PHA = [3:8];
nb = 10;

phase_edges=linspace(-pi, pi, nb+1);
t1 = datetime; 
for subji = 1:length(allPhases)
    oneListPhases           = allPhases{subji};
    oneListFits             = nnFit{subji}; 
    for chani = 1:size(oneListPhases, 1)
        olip = squeeze(oneListPhases(chani, :, :))'; 
        olif = squeeze(oneListFits)'; 
        for triali = 1:size(olip, 1)
            pha2use = olip(triali,:);  
            for hi=1:nb
                d2p = pha2use> phase_edges(hi) & pha2use < phase_edges(hi+1) ;
                %all_phaH{subji, hi}( :, triali) = d2p;
                j2c(d2p) = hi; %check steps
                fitPHA{subji,:}(chani, triali, hi) = mean(olif(d2p, triali), 'omitnan');
            end
        end
    end
end

t2 = datetime; 
etime(datevec(t2), datevec(t1))

%% build distributions for all VVS electrodes
clearvars -except allPhases nnFit

pfc_fits   = nnFit([2 3  5  9 10 11 12 14 15 16]);
vvs_phases = allPhases([7 9 13 18 19 20 21 23 27 28])'; 


f_PHA = [3:8];
nb = 10;
phase_edges=linspace(-pi, pi, nb+1);
        

t1 = datetime; 
for subji = 1:length(vvs_phases)
    oneListPhases           = vvs_phases{subji}(:,2000:2500,:);
    oneListFits             = pfc_fits{subji}(:,:, 2000:2500); 
    for chani = 1:size(oneListPhases, 1)
        olip = squeeze(oneListPhases(chani, :, :))'; 
        olif = squeeze(oneListFits)'; 
        for triali = 1:size(olip, 1)
            pha2use = olip(triali,:);  
            for hi=1:nb
                d2p = pha2use> phase_edges(hi) & pha2use < phase_edges(hi+1) ;
                pha2check{subji}(d2p) = hi;
                j2c(d2p) = hi; %check steps
                fitPHA{subji,:}(chani, triali, hi) = mean(olif(d2p, triali), 'omitnan');
            end
        end
    end
end

t2 = datetime; 
etime(datevec(t2), datevec(t1))

%%
d2p = pha2check{2}; 
histogram(d2p)
[count] = hist(d2p); 
sum(count)


%% calculate distance to uniform for all electrodes and trials

tic
d2n = ones(nb,1)/nb; d2n = d2n; %uniform distribution
clear allKLs
for subji =  1:length(fitPHA)
    ph2p = fitPHA{subji};
    for chani = 1:size(ph2p, 1)
        for triali = 1:size(ph2p, 2)
            d1 =  squeeze(ph2p(chani,triali, :)); 
            d1n = ( d1-min(d1) )  ./  (max(d1) - min(d1)) + eps; 
            d1n = d1n ./ sum(d1n); 
            allKLs{subji,:}(chani, triali) = kldiv ((1:nb)', d1n, d2n);
        end
    end
end

toc

%% mean KL for each subject and VVS channel
mKL = cellfun(@(x) mean(x, 2), allKLs, 'un', 0);


%% PERMUTATIONS VVS electrodes
clearvars -except allPhases nnFit

pfc_fits   = nnFit([2 3  5  9 10 11 12 14 15 16]);
vvs_phases = allPhases([7 9 13 18 19 20 21 23 27 28])'; 

nPerm = 1000; 
nb = 10;
phase_edges=linspace(-pi, pi, nb+1);
        

t1 = datetime;   
for subji = 1:length(vvs_phases)
   
   phaPerms        = zeros(nPerm, size(vvs_phases{subji}, 1), size(vvs_phases{subji}, 3), nb);
   for permi = 1:nPerm
        oneListPhases   = vvs_phases{subji};
        oneListFits     = pfc_fits{subji}; 
        for chani = 1:size(oneListPhases, 1)
            olip = squeeze(oneListPhases(chani, :, :))'; 
            olif = squeeze(oneListFits)'; 
            for triali = 1:size(olip, 1)
                pha2use = olip(triali,:);  
                %shuffle
                shiftBy =round(length(pha2use)*rand());
                pha2use_perm = circshift(pha2use, [0,shiftBy]);
                for hi=1:nb
                    d2p = pha2use_perm> phase_edges(hi) & pha2use_perm < phase_edges(hi+1) ;
                    phaPerms(permi, chani, triali, hi)    = mean(olif(d2p, triali), 'omitnan');
                end
            end
        end
   end
    fitPHA_PERM{subji} = phaPerms; 
end

save('fitPHA_PERM_0-500ms_VVS_1000p', 'fitPHA_PERM');

t2 = datetime; 
etime(datevec(t2), datevec(t1))




%% check example perm
permi = 29; 
subji  = 1;
chani  = 5;
triali = 1;


ph2p = fitPHA_PERM{subji};
d2p = squeeze(ph2p(permi, chani,triali, :))'; 
cfg.nb = 10;  
cfg.rlimi = [-.015 .025];
figure()
plot_circ_hist(d2p, cfg);
d1 =  d2p;
d1n = ( d1-min(d1) )  ./  (max(d1) - min(d1)) + eps; 
d1n = d1n ./ sum(d1n); 
d2n = ones(cfg.nbnb,1)/cfg.nb; d2n = d2n';
figure()
cfg.rlimi = [-.1 .3];
plot_circ_hist(d1n, cfg);
kl = kldiv (1:cfg.nb, d1n, d2n)
filename = [num2str(subji) '_' num2str(chani) '.png'];
%exportgraphics(gcf, filename, 'Resolution', 300);
%close all

%% calculate distance to uniform for all combinations
tic
nb = 10; 
d2n = ones(nb,1)/nb; d2n = d2n; %uniform distribution
clear allKLPERMs
for subji =  1:length(fitPHA_PERM)
    ph2p = fitPHA_PERM{subji};
    for permi = 1:size(ph2p, 1)
        for chani = 1:size(ph2p, 2)
            for triali = 1:size(ph2p, 3)
                d1 =  squeeze(ph2p(permi, chani,triali, :)); 
                d1n = ( d1-min(d1) )  ./  (max(d1) - min(d1)) + eps; 
                d1n = d1n ./ sum(d1n); 
                allKLPERMs{subji,:}(permi, chani, triali) = kldiv ((1:nb)', d1n, d2n);
            end
        end
    end
end

toc

%% z-score observed KLs wiht distribution of KL scores for each subject, trial and channel
load allKLPERMs 
load allKLs
for subji = 1:1%10
    mKLp = allKLPERMs{subji};
    mKLs = allKLs{subji};
    for chani = 1:1%size(mKLp, 2)
        for triali = 1:size(mKLp, 3)
            x = mKLs(chani, triali); 
            y = mKLp(:, chani, triali); 
            zKL = (x - mean(y) ) / std(y);
            zKLs{subji,:}(chani, triali) = zKL; 
        end
    end    
end

%% mean zKL for each subject and VVS channel
mZKL = cellfun(@(x) mean(x, 2), zKLs, 'un', 0);


%% for each electrode 1000 values
nPerm = 1000
clear mKLp
for subji =  1:length(fitPHA_PERM)
    for permi = 1:nPerm
        dperm = squeeze(allKLPERMs{subji}(permi, :, :)); 
        %mean KL for each subject and VVS channel
        mKLp{subji,:}(permi, :) = mean(dperm, 2);
    end
end

%% check rank of each channel
clear allPs
for subji = 1:length(mKL)
    obsKLp = mKL{subji}; 
    disKLp = mKLp{subji}; 
    for chani = 1:length(obsKLp)
        obsKL = obsKLp(chani);
        disKL = disKLp(:, chani); 
        
        allAb = disKL(abs(disKL) > obsKL);
        p =1 - (nPerm - (length (allAb)+1) )  /nPerm;
        allPs(chani, subji) = p; 
        %figure()
        %histogram(disKL); hold on; 
        %plot ([obsKL obsKL], get(gca, 'ylim'), 'LineWidth', 2,'Color', [.9 .3 .1] );
    end
end



%% z-score observed KLs wiht distribution of KL scores for each subject, trial and channel

load mKL

for subji = 1:10
    mKLp = allKLPERMs{subji};
    mKLs = allKLs{subji};
    for chani = 1:size(mKLp, 2)
        for triali = 1:size(mKLp, 3)
            x = mKLs(chani, triali); 
            y = mKLp(:, chani, triali); 
            zKL = (x - mean(y) ) / std(y);
            zKLs{subji,:}(chani, triali) = zKL; 
        end
    end    
end



%%
cfg.nb = 10;  
cfg.rlimi = [-.01 .03];
plot_circ_hist(d2p, cfg);




%%
nb = 10; 
u1 = d2p; 
u1n = ( u1-min(u1) )  ./  (max(u1) - min(u1)) + eps; u1n = u1n ./ sum(u1n); 
u2n = ones(nb,1)/nb; u2n = u2n';
kl = kldiv (1:nb, u1n, u2n)

%% 
figure()
bar(d2p); 
figure()
bar(u1n)













%% higher temporal resolution no time lims
clustInfo = bwconncomp(h);
clength = cell2mat(cellfun(@(x) length(x), clustInfo.PixelIdxList, 'un', 0))';
sort(clength)

h1 = zeros(size(h));
%h1(clustInfo.PixelIdxList{311}) = 1; 
h1(clustInfo.PixelIdxList{27}) = 1; 

figure()
myCmap = colormap(brewermap([],'*Spectral'));
contourf(t, 100, 'linecolor', 'none'); hold on; 
contour(h1, 1, 'Color', [0, 0, 0], 'LineWidth', 1);
set(gca, 'FontSize', 12)
exportgraphics(gcf, 'allM.png', 'Resolution', 300);


%% plot 
d2p = squeeze(oneListPhases(1,:,1))
plot(d2p)


%% only for the small matrix to be scaled
figure()
myCmap = colormap(brewermap([],'*Spectral'));
freqs = 1:520;
contourf(1:870, freqs, myresizem(t, 10), 100, 'linecolor', 'none'); hold on; 
contour(1:870, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 1); colorbar
set(gca, 'xlim', [100 600], 'clim', [-4 4])
exportgraphics(gcf, 'allM.png', 'Resolution', 300);

%% higher temporal resolution
figure()
myCmap = colormap(brewermap([],'*Spectral'));
times = -1.76:.01:6.75;
freqs = 1:52;
contourf(times, freqs, t, 100, 'linecolor', 'none'); hold on; 
contour(times, freqs, h, 1, 'Color', [0, 0, 0], 'LineWidth', 1);
plot([0 0], get(gca, 'ylim'), 'k', 'LineWidth', 2);
set(gca, 'xlim', [-1 4], 'clim', [-4 4], 'ytick', [1 30 52], 'yticklabels', {'3' '30' '150'}); colorbar
set(gca, 'FontSize', 12)
exportgraphics(gcf, 'allM.png', 'Resolution', 300);




%%



















%%
brainROI    = 'PFC'; 
lays2load   = [56];
f           = [3:54];

clear nnFit
for subji = 1:length(allTraces)
   
    oneListTraces     = allTraces{subji};
    oneListIDs          = allIDs{subji};
    [oneListPow]        = extract_power_WM (oneListTraces, 0.1); % all not possible (too large)
    oneListPow          = norm_trials(oneListPow);
    oneListPhases       = extract_phases_WM(oneListTraces, f); 
    idsh                = cell2mat(cellfun(@(x) ~strcmp(x(6), '4'), oneListIDs, 'un', 0));
    oneListIDs          = oneListIDs(idsh);
    oneListPow          = oneListPow(idsh,:,:,:);
    neuralRDMs          = createNeuralRDMs(oneListPow, f, 1, 1);
    networkRDMs         = createNetworkRDMs(oneListIDs, lays2load, brainROI, subji); %layers to load
    
    nnFit{subji,:}      = fitModel_WM(neuralRDMs, networkRDMs); %layers*frequencies*times
    
end
save ('nnFit', 'nnFit');




%% Phase extraction with EEGLAB (sanity check)
EEG.data    = oneListTraces_c;
EEG.trials  = size(oneListTraces_c, 3); 
EEG.nbchan  = 1; 
EEG.srate   = 1000;
EEG.pnts    = 9000;
EEG         = pop_eegfiltnew (EEG, f(1),f(end));

d2pp = EEG.data(1,:,1);
%plot(d2pp)

data        = squeeze(d2pp); % until 1s
d2p      = angle(hilbert(data));

figure()
plot(d2p(1,:))




%% higher temporal resolution 1-1 10ms
figure()
myCmap = colormap(brewermap([],'*Spectral'));
times = -1.76:.01:7.249
freqs = 1:52;
contourf(times, freqs, t, 100, 'linecolor', 'none'); hold on; 
contour(times, freqs, h, 1, 'Color', [0, 0, 0], 'LineWidth', 1);
plot([0 0], get(gca, 'ylim'), 'k', 'LineWidth', 2);
set(gca, 'xlim', [-1 4], 'clim', [-4 4], 'ytick', [1 30 52], 'yticklabels', {'3' '30' '150'}); colorbar
set(gca, 'FontSize', 12)
exportgraphics(gcf, 'allM.png', 'Resolution', 300);


%% 
nTimes = 901; 
mf = 1; 
win_width = 50; 
bins  =  floor ( (nTimes/mf)- win_width/mf+1 );
for timei = 1:bins %parfor possible
    timeBins(timei,:) = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
end

%% 
figure()
imagesc(squeeze(networkRDMs)); axis square ;colorbar

%%
d2p = squeeze(networkRDMs); 




%%
plot([.45 .45], get(gca, 'ylim'), 'k:', 'LineWidth', 12);
plot([.95 .95], get(gca, 'ylim'), 'k:', 'LineWidth', 12);
%plot([1.25 1.25], get(gca, 'ylim'), 'k:', 'LineWidth', 12);
plot(get(gca, 'xlim'),[300 300], 'k:', 'LineWidth', 12);
set(gca, 'xtick', [-0.055 .945 1.945 2.945 3.945 ], 'xticklabel', [0 1 2 3 4], 'xlim', [0 3.5], 'clim', [-4 4])
set(gca,'XTick',[], 'YTick', [])
set(gca,'xticklabel',[], 'FontSize', 12)




%% 
figure()
d2p = squeeze(neuralRDMs(: ,:, 25));
imagesc(d2p)




%% 
figure()
d2p = squeeze(oneListPow(1, 1, : ,:));
imagesc(d2p)




%% load raw traces in VVS

cd E:\daniel\_WM\analysis\predictor_matrix_analysis
currentF = pwd;
cd E:\daniel\_WM\analysis\out_contrasts\raw_traces\vvs\allTrials
sublist = dir('*_out_contr.mat'); sublist = {sublist.name};
disp (['measurements -> ' num2str(length(sublist))]);

for subji=1:length(sublist)
    
    load(sublist{subji});
    
    ids = cell2mat(cellfun(@(x) strcmp(x(1), '7') & ~strcmp(x(6), '4'), cfg_contrasts.oneListIds_c, 'un', 0));
    allTraces{subji,1} = cfg_contrasts.oneListTraces(:, :, ids); 
    
    allIDs{subji, 1} = cfg_contrasts.oneListIds_c(ids); 
    allChans{subji, 1} = cfg_contrasts.chanNames; 
    

end
cd (currentF)
clearvars -except allTraces allIDs allChans





%% load PFC 
load RNN_pfc_M_noAv_R_T_1-550_C_aFV_10ms_1-1win % only layer 56 was computed
%size(all_r_Times_Trials{1}) % only one layer (8th layer 8th time point)

%% align
pfc = all_r_Times_Trials([2 3  5  9 10 11 12 14 15 16])';
vvs = allTraces([7 9 13 18 19 20 21 23 27 28]);
lenVVS = cellfun(@(x) size(x), pfc, 'un', 0)

%% 

for subji = 1:10
    
    pfc_h = squeeze(pfc{subji});
    pfc_b = squeeze(mean(pfc_h(:, 18:28, :), 2)); %beta range
    
    allSB(subji, :) = mean(pfc_b);
    
    
    
    
    
    
    
end



%% 
%plot(pfc_b(2,:)); % single trial
%plot(mean(pfc_b)); % mean

m = mean(allSB);
se = std(allSB) / sqrt(10);
figure()
errorbar(m, se); % group mean




%% plot traces sanity check
%load RNN_pfc_M&C_noAv_54_Real_nT_1-49_C_aFV


d2pp = mean(all_r_Times(:, 56, 13:29, :),3, 'omitnan');
d2p = squeeze(d2pp);

m = mean(d2p, 'omitnan');
se = std(d2p, 'omitnan') / sqrt(size(d2p, 1));
figure()
errorbar(m, se); % group mean

figure()
d2p = squeeze(mean(all_r_Times(:,56,:,:), 'omitnan'));
imagesc(d2p);

%% plot each subject 
d2pp = all_r_Times(:, 56, :, :);
d2p = squeeze(d2pp(12, :, :, :));
figure()
imagesc(d2p)




%% plot all subjects 
%load RNN_pfc_M_noAv_R_T_1-49_C_aFV_100ms_5-1win
%load RNN_pfc_M_noAv_R_T_1-550_C_aFV_10ms_1-1win
load RNN_pfc_M&C_noAv_54_Real_T_1-49_C_aFV
size(all_r_Times_Trials{1})
%pfc = cellfun(@(x) squeeze(mean(x,2)), all_r_Times_Trials, 'un', 0)';
pfc = cellfun(@(x) squeeze(mean(x(7, :, :, :),2)), all_r_Times_Trials, 'un', 0)';

%pfc2 = vertcat(pfc{:});
pfc2 = [pfc{:}];
pfc3 = reshape(pfc2, [], 52, 16);
pfc3 = permute(pfc3, [3 2 1]);

%d2p = squeeze(pfc3(1, :, :));  %check that reshape is correct (n1 = all %nans)
d2p = squeeze(mean(pfc3, 'omitnan'));
figure()
imagesc(d2p);colorbar



%% 
for subji = 1:16
    d2p = squeeze(pfc3(subji, :, :));  %check that reshape is correct (n1 = all %nans) 
    figure()
    imagesc(d2p);colorbar

end




%% create similarity time-series


currentF = pwd;
cd D:\_WM\analysis\phase_RSA\SISC_EM2\3-54Hz

sublist = dir('*_rsa.mat');
sublist = {sublist.name};
disp (['measurements -> ' num2str(length(sublist))]);

for subji=1:length(sublist)
    
    load(sublist{subji});
    all{subji,1} = rsaZ; 
    if exist('allIDs')
        all_IDs{subji, :} = allIDs; 
    end

end


cd (currentF)






















%% 