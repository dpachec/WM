%% Calculate from epoched raw traces
%% first load traces
clear
%ROI__layers__freqs__avRepet__avTimeFeatVect__freqResolv(0-1)__fitMode(0:noTrials; 1:Trials)__timeRes__win-width__mf
%example f2sav = 'RNN_pfc_E_[8:8:56]_3-54_1_0_0_0_.1_5_1.mat'; 
f2sav = 'RNN_vvs_E_[8:8:56]_3-54_1_0_0_0_.1_5_1.mat'; 
cfg = getParams(f2sav);
f2t = strsplit(f2sav, '_');
region = f2t{2};
paths = load_paths_WM(region);
filelistSess = getFiles(paths.traces);

t1 = datetime; 

for sessi= 1:length(filelistSess) %this one starts at 1 and not at 3
    disp(['File > ' num2str(sessi)]);
    load([paths.traces filelistSess{sessi}]);   
   
    if strcmp(cfg.period, 'M')
        ids = cell2mat(cellfun(@(x) strcmp(x(1), '7') & ~strcmp(x(6), '4'), cfg_contrasts.oneListIds_c, 'un', 0));
    elseif strcmp(cfg.period, 'E')
        ids = cell2mat(cellfun(@(x) (strcmp(x(1), '1') | strcmp(x(1), '2') | strcmp(x(1), '3')) & ~strcmp(x(6), '4'), cfg_contrasts.oneListIds_c, 'un', 0));
    end
    oneListTraces = cfg_contrasts.oneListTraces(:,:,ids);
    cfg_contrasts.oneListIds_c    = cfg_contrasts.oneListIds_c(ids); 
    cfg_contrasts.oneListPow    = extract_power_WM (oneListTraces, cfg.timeRes); % 
    cfg_contrasts = normalize_WM(cfg_contrasts, 1, 'sess', []);
    if (cfg.avRep)
        cfg_contrasts               = average_repetitions(cfg_contrasts);
    end

    neuralRDMs                  = createNeuralRDMs(cfg_contrasts.oneListPow, cfg.freqs, cfg.win_width, cfg.mf, cfg.fR, cfg.avTFV);
    networkRDMs                 = createNetworkRDMs(cfg_contrasts.oneListIds_c, cfg.net2load, cfg.lays2load, cfg.brainROI, sessi, paths, cfg.period); %layers to load
    
    nnFit{sessi,1}              = fitModel_WM(neuralRDMs, networkRDMs, cfg.fitMode); 
    nnFit{sessi,2}              = cfg_contrasts.oneListIds_c; 
    
end

mkdir ([paths.results.DNNs]);
save([paths.results.DNNs f2sav], 'nnFit');

t2 = datetime; 
etime(datevec(t2), datevec(t1))




%% load file
%ROI__layers__freqs__avRepet__avTFV__fRes(0-1)__fitMode(0:noTrials; 1:Trials)__timeRes__win-width__mf
clear 
f2sav = 'RNN_pfc_E_[8:8:56]_3-54_1_0_0_0_.1_5_1.mat'; 
f2t = strsplit(f2sav, '_');
region = f2t{2};

paths = load_paths_WM(region);
load([paths.results.DNNs f2sav]);



%% all plot cells start with nnH

sub2exc = [];

for subji = 1:length(nnFit)
    
   %nnH(subji, : ,:) = nnFit{subji, 1}(1,:,:);
   nnH(subji, : ,:) = nnFit{subji, 1}(4,:);
   %nnH(subji, : ,:,:) = nnFit{subji};
        
end


nnH(sub2exc, :, :) = []; 
nnH = squeeze(nnH);
[h p ci ts] = ttest(nnH);
h = squeeze(h); t = squeeze(ts.tstat);



%% plot frequency resolved
d2p = squeeze(mean(nnH, 'omitnan'));
figure
freqs = 1:520; 
times = -1.75:.01:6.849; 
clustinfo = bwconncomp(h);
for pixi = 1:length(clustinfo.PixelIdxList)
   h(clustinfo.PixelIdxList{pixi}) = 0;   
end
h(clustinfo.PixelIdxList{27}) = 1;
tObs = sum(t(clustinfo.PixelIdxList{27}))

contourf(times, freqs, myresizem(t, 10), 100, 'linecolor', 'none'); hold on; %colorbar
contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 2);
set(gca, 'xlim', [-1 6])
%set(gca, 'clim', [-.025 .025])

%clustinfo = bwconncomp(h);




%% all plot cells start with nnH

sub2exc = [];

for subji = 1:length(nnFit)
    
   %nnH(subji, : ,:) = nnFit{subji, 1}(1,:,:);
   nnH(subji, : ,:) = nnFit{subji, 1}(5,:);
   %nnH(subji, : ,:,:) = nnFit{subji};
        
end


nnH(sub2exc, :, :) = []; 
nnH = squeeze(nnH);
[h p ci ts] = ttest(nnH);
h = squeeze(h); t = squeeze(ts.tstat);

%% plot bands

d2p = squeeze(mean(nnH, 'omitnan'));

times = -1.75:.1:6.849; 
figure
plot(times, d2p); hold on; 
%h(h==0) = nan; h(h==1) = .02;
%plot(h, 'lineWidth', 2)
%set(gca, 'xlim', [0 45])
%set(gca, 'clim', [-.025 .025])

[h p ci ts] = ttest(nnH);
h = squeeze(h); t = squeeze(ts.tstat);


%% plot bands all layers
clear 
f2sav = 'vvs_1-56_3-54_1_0_0_0_.1_5_1.mat'; 
f2t = strsplit(f2sav, '_');
region = f2t{1};

paths = load_paths_WM(region);
load([paths.results.DNNs f2sav]);

for subji = 1:length(nnFit)
   for layi = 1:56
      nnH(subji,layi, :) = nnFit{subji, 1}(layi,:) ;
   end
end


%% 

figure()
for layi = 1:56
    d2p = squeeze(mean(nnH(:, layi,:), 'omitnan'));
    plot(d2p); hold on; 
end
%h(h==0) = nan; h(h==1) = .02;
%plot(h, 'lineWidth', 2)


%% PLOT OBS DATA 

sub2exc = [];

for subji = 1:length(nnFit)
   
   nnH(subji, : ,:,:) = nnFit{subji};
        
end


nnH(sub2exc, :, :) = []; 
nnH = squeeze(nnH);
[h p ci ts] = ttest(nnH);
h = squeeze(h); t = squeeze(ts.tstat);

sub2exc =[];

all_r_Times = nnH; 

nSubjs =size(all_r_Times, 1);  
nLays = size(all_r_Times, 2); 
nFreqs = size(all_r_Times, 3); 
clear hL tL clustinfo
for freqi = 1:nFreqs
    allR = squeeze(all_r_Times(:, :, freqi, :));
    allR(sub2exc,:,:) = []; 
    for layi = 1:nLays
        [hL(layi, freqi, :) p ci ts] = ttest(allR(:,layi,:)); 
        tL(layi, freqi, :) = ts.tstat;
    end
end

f = 1:52;

clear max_clust_sum_obs allSTs
for layi = 1:nLays
    
    hLay = squeeze(hL(layi,f,:)); 
    tLay = squeeze(tL(layi,f,:));
    
    clear allSTs  
    clustinfo = bwconncomp(hLay);
    for pxi = 1:length(clustinfo.PixelIdxList)
        % check whether it is a combined + and - cluster
        V = tLay(clustinfo.PixelIdxList{pxi});
        Vids = clustinfo.PixelIdxList{pxi}; 
        if ~any(diff(sign(V(V~=0)))) %all time bins have tvalues of same sie
            allSTs(pxi,:) = sum(tLay(clustinfo.PixelIdxList{pxi}));
        else %remove the 
            big0 = V>0; small0 = V<0; 
            VidsS = Vids(small0); 
            ids2k = Vids(big0); 
            if sum(big0) > sum(small0)
                V(small0) = NaN; 
                hLay(VidsS) = 0; 
                clustinfo.PixelIdxList{pxi} = ids2k; 
            else
                %V(big0) = NaN; 
                %hLay(big0,:) = NaN; 
            end
            
            allSTs(pxi,:) = sum(V, 'omitnan'); 
        end
    end

    [maxh id] = max(abs(allSTs));
    max_clust_sum_obs(layi,:) = allSTs(id); 

    % % % % rem non-sig-clusters
    for ci = 1:length(clustinfo.PixelIdxList)
        modifiedCluster = abs(sum(tLay(clustinfo.PixelIdxList{ci})))+0.0001; %add this small number because removing negative values creates slighlly different valeus
       if modifiedCluster < max(abs(allSTs))   + 1    %add +1 at the very end to delete all clusters
          %disp(['layi > ' num2str(layi) ' > ' num2str(abs(sum(tLay(clustinfo.PixelIdxList{ci})))) '_' num2str(max(abs(allSTs)))]) 
          %hLay (clustinfo.PixelIdxList{ci}) =  0; 
       end
    end
    
    hL(layi, f, :) = hLay; 
    tL(layi, f, :) = tLay; 
end



figure()
layT = tiledlayout(3, 3);
layT.TileSpacing = 'compact';
layT.Padding = 'compact';
for layi = 1:nLays
    
    t = squeeze(tL(layi,:,:));
    h = squeeze(hL(layi,:,:));
    
    if size(all_r_Times, 4) < 30 & size(all_r_Times, 4) > 5
        % % % % encoding 
        if nLays > 9
            subplot(7, 8, layi)
        else
            subplot(3, 3, layi)
        end
        %times = -.5:.01:1.499;
        times = 0:.01:1.499;
        %freqs = [.1:.1:29.9 30:.5:150]; 
        freqs = 1:520;
        contourf(times, freqs, myresizem(t, 10), 100, 'linecolor', 'none'); hold on; %colorbar
        contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 2);
        plot([-0.055 -0.055], get(gca, 'ylim'), 'k:', 'LineWidth', 2);
        plot(get(gca, 'xlim'),[300 300], 'k:', 'LineWidth', 2);
        set(gca, 'xtick', [-0.055 .445 .945 1.45], 'xticklabel', [0 .5 1 1.5], 'clim', [-5 5]); colormap jet
        set(gca,'XTick',[], 'YTick', [], 'xticklabel',[], 'FontSize', 12)

    elseif size(all_r_Times, 4) == 5
        % % % % encoding 
        if nLays > 9
            subplot(7, 8, layi)
        else
            subplot(3, 3, layi)
        end
        times = 0:.01:0.499;
        %freqs = [.1:.1:29.9 30:.5:150]; 
        freqs = 1:520;
        contourf(times, freqs, myresizem(t, 10), 100, 'linecolor', 'none'); hold on; %colorbar
        contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 2);
        plot([-0.055 -0.055], get(gca, 'ylim'), 'k:', 'LineWidth', 2);
        plot(get(gca, 'xlim'),[300 300], 'k:', 'LineWidth', 2);
        set(gca, 'xtick', [-0.055 .445 .945 1.45], 'xticklabel', [0 .5 1 1.5], 'clim', [-5 5]); colormap jet
        set(gca,'XTick',[], 'YTick', [], 'xticklabel',[], 'FontSize', 12)
    elseif size(all_r_Times, 4) == 546
        % % % % encoding 
        if nLays > 9
            subplot(7, 8, layi)
        else
            subplot(3, 3, layi)
        end
        times = 0:.01:5.009;
        %freqs = [.1:.1:29.9 30:.5:150]; 
        freqs = 1:52;
        contourf(times, freqs, t, 100, 'linecolor', 'none'); hold on; %colorbar
        contour(times, freqs, h, 1, 'Color', [0, 0, 0], 'LineWidth', 2);
        %plot([-0.055 -0.055], get(gca, 'ylim'), 'k:', 'LineWidth', 2);
        %plot(get(gca, 'xlim'),[300 300], 'k:', 'LineWidth', 2);
        %set(gca, 'xtick', [-0.055 .445 .945 1.45], 'xticklabel', [0 .5 1 1.5], 'clim', [-5 5]); colormap jet
        %set(gca,'XTick',[], 'YTick', [], 'xticklabel',[], 'FontSize', 12)
    else
        if nLays > 9
            subplot(7, 8, layi)
        else
            subplot(3, 3, layi)
        end
        times = -1.75:.01:6.849; 
        freqs = 1:520;
        contourf(times, freqs, myresizem(t, 10), 100, 'linecolor', 'none'); hold on; 
        contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 1); %colorbar
        plot([-0.055 -0.055], get(gca, 'ylim'), 'k:', 'LineWidth', 2);
        plot(get(gca, 'xlim'),[300 300], 'k:', 'LineWidth', 2);
        set(gca, 'xtick', [-0.055 .945 1.945 2.945 3.945 ], 'xticklabel', [0 1 2 3 4], 'clim', [-3 3]); % 'xlim', [0 3.5],
        set(gca, 'xlim', [-1 6])
        set(gca,'XTick',[], 'YTick', [])
        set(gca,'xticklabel',[], 'FontSize', 12)
        colormap jet
    end

end
        

%f2p = [num2str(layi, '%02.f') '.png'];
%exportgraphics(gcf, f2p, 'Resolution', 300)






%% PERMUTATIONS
clear
nPerm = 1;
%ROI__layers__freqs__avRepet__avTimeFeatVect__freqResolv(0-1)__fitMode(0:noTrials; 1:Trials)__timeRes__win-width__mf
%example f2sav = 'pfc_8-16-24-32-40-48-56_13-29_0_1_500_1_1'; 
%f2sav = ['RNN_pfc_M_56_3-54_1_0_1_0_.1_5_1_p' num2str(nPerm) '.mat']; 
f2sav = ['RNN_vvs_M_1-56_3-54_1_0_1_0_.1_5_1_p' num2str(nPerm) '.mat']; 

cfg = getParams(f2sav);
f2t = strsplit(f2sav, '_');
region = f2t{2};
paths = load_paths_WM(region);
filelistSess = getFiles(paths.traces);

t1 = datetime; 

for sessi= 1:length(filelistSess) %this one starts at 1 and not at 3
    disp(['File > ' num2str(sessi)]);
    load([paths.traces filelistSess{sessi}]);   
   
    if strcmp(cfg.period, 'M')
        ids = cell2mat(cellfun(@(x) strcmp(x(1), '7') & ~strcmp(x(6), '4'), cfg_contrasts.oneListIds_c, 'un', 0));
    elseif strcmp(cfg.period, 'E')
        ids = cell2mat(cellfun(@(x) (strcmp(x(1), '1') | strcmp(x(1), '2') | strcmp(x(1), '3')) & ~strcmp(x(6), '4'), cfg_contrasts.oneListIds_c, 'un', 0));
    end
    oneListTraces = cfg_contrasts.oneListTraces(:,:,ids);
    cfg_contrasts.oneListIds_c    = cfg_contrasts.oneListIds_c(ids); 
    cfg_contrasts.oneListPow    = extract_power_WM (oneListTraces, cfg.timeRes); % 
    cfg_contrasts = normalize_WM(cfg_contrasts, 1, 'sess', []);
    if (cfg.avRep)
        cfg_contrasts               = average_repetitions(cfg_contrasts);
    end

    neuralRDMs                  = createNeuralRDMs(cfg_contrasts.oneListPow, cfg.freqs, cfg.win_width, cfg.mf, cfg.fR, cfg.avTFV);
    networkRDMs                 = createNetworkRDMs(cfg_contrasts.oneListIds_c, cfg.net2load, cfg.lays2load, cfg.brainROI, sessi, paths, cfg.period); %layers to load
    
    for permi = 1:nPerm
        sC = size(networkRDMs, 2);
        ids = randperm(sC);
        networkRDMs = networkRDMs(:, ids, ids); 
        nnFitPerm(permi, sessi,:,:, :)              = fitModel_WM(neuralRDMs, networkRDMs, cfg.fitMode); 
    end
    
end

mkdir ([paths.results.DNNs]);
save([paths.results.DNNs f2sav], 'nnFitPerm');

t2 = datetime; 
etime(datevec(t2), datevec(t1))


%% compute clusters in each permutation

nPerm = 100; 
for permi = 1:nPerm
    
    dataP = squeeze(nnFitPerm(permi, :, :,21:55));
    [h p ci ts] = ttest(dataP);
    h = squeeze(h); t = squeeze(ts.tstat);
    
    clear allSTs  
    clustinfo = bwconncomp(h);
    if length(clustinfo.PixelIdxList) > 0
        for pxi = 1:length(clustinfo.PixelIdxList)
           allSTs(permi, pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
        end
    else
        allSTs(permi, :) = 0;
    end

    [maxh id] = max(abs(allSTs(permi)));
    max_clust_sum_perm(permi,:) = allSTs(permi,id); 

end





 

%% FREQUENCY RESOLVED DNN ANALYSIS (SLOW WAY AS SANITY CHECK) 
%% 
clearvars -except act_CH act_FR 


f2sav       = 'RNN_pfc_M_Av_R_nT_1-20_C_naFV_100ms'; %'Alex_vvs_CueAll_noAv_54_Real_nT_1-49'


load_parameters_WMFRA;
loadNet_WM;
region = f2t{2};
paths = load_paths_WM(region);
filelistSess = getFiles(paths.out_contrasts_path);
nSubj = length(filelistSess); 
all_r_Times    = zeros(nSubj, nLays, nFreqs, nTimes);

tic

for subji= 1:nSubj %start with subject 2
    disp (['Subj: ' num2str(subji)])
    load([paths.out_contrasts_path filelistSess{subji}]);   
    cfg_contrasts.oneListPow = cfg_contrasts.oneListPow(:, :, :, resTime); 
    cfg_contrasts = normalize_WM(cfg_contrasts, 1, 'blo', [3 7]);
    if strcmp(avMeth,'pow')
       cfg_contrasts = average_repetitions(cfg_contrasts);
    end
    
    for freqi = 1:length(freqs2test)
        f  = freqs2test(freqi);
        
        if strcmp(f2t{9}, 'aFV')
            avTW = 1;
        else
            avTW = 0;
        end
        rdm_prev = create_rdms(cfg_contrasts, f, it, 5, 1, avTW);%
        
       
        [allS ids act_FR2 act_CH2 cM] = getIDsAndnRDMs(rdm_prev, act_CH, act_FR);
        
        for timei = 1:size(rdm_prev.rdm, 3)
            rdmp = squeeze(rdm_prev.rdm);
            rdm = squeeze(rdmp(:,:,timei));
            rdm(rdm ==1) = 1000;rdm(rdm ==0) = 2000;
            rdm = tril(rdm, -1);rdm(rdm==0) =[];
            rdm(rdm==1000) = 1;rdm(rdm==2000) = 0;
        
            for layi = 1:nLays
                if subji < subj_ch_fr
                    M =  squeeze(act_FR2(layi,:,:)); 
                else
                    M =  squeeze(act_CH2(layi,:,:)); 
                end
                M(M ==1) = 1000;M(M ==0) = 2000;
                M = tril(M, -1);M(M==0) =[];
                M(M==1000) = 1;M(M==2000) = 0;

                cMh = cM; 
                cMh(cMh ==1) = 1000;cMh(cMh ==0) = 2000;
                cMh = tril(cMh, -1);cMh(cMh==0) =[];
                cMh(cMh==1000) = 1;cMh(cMh==2000) = 0;

                
                if strcmp(f2t{8}, 'PC')
                    x = [rdm' M']; 
                    z = cMh'; 
                    allTEst_prev = partialcorr(x, z, 'type', 's');
                    allTEst = allTEst_prev(2);
                else
                    allTEst = corr(rdm', M', 'type', 's');
                end
                
                

                all_r_Times(subji, layi,freqi,timei) = allTEst;  
            end
        end
    end
    
    
    
end 

mkdir ([paths.results.DNNs]);
save([paths.results.DNNs f2sav], 'all_r_Times');

toc



%% PLOT OBS DATA 

sub2exc =[];

all_r_Times = nnH; 

nSubjs =size(all_r_Times, 1);  
nLays = size(all_r_Times, 2); 
nFreqs = size(all_r_Times, 3); 
clear hL tL clustinfo
for freqi = 1:nFreqs
    allR = squeeze(all_r_Times(:, :, freqi, :));
    allR(sub2exc,:,:) = []; 
    for layi = 1:nLays
        [hL(layi, freqi, :) p ci ts] = ttest(allR(:,layi,:)); 
        tL(layi, freqi, :) = ts.tstat;
    end
end

f = 1:52;

clear max_clust_sum_obs allSTs
for layi = 1:nLays
    
    hLay = squeeze(hL(layi,f,:)); 
    tLay = squeeze(tL(layi,f,:));
    
    clear allSTs  
    clustinfo = bwconncomp(hLay);
    for pxi = 1:length(clustinfo.PixelIdxList)
        % check whether it is a combined + and - cluster
        V = tLay(clustinfo.PixelIdxList{pxi});
        Vids = clustinfo.PixelIdxList{pxi}; 
        if ~any(diff(sign(V(V~=0)))) %all time bins have tvalues of same sie
            allSTs(pxi,:) = sum(tLay(clustinfo.PixelIdxList{pxi}));
        else %remove the 
            big0 = V>0; small0 = V<0; 
            VidsS = Vids(small0); 
            ids2k = Vids(big0); 
            if sum(big0) > sum(small0)
                V(small0) = NaN; 
                hLay(VidsS) = 0; 
                clustinfo.PixelIdxList{pxi} = ids2k; 
            else
                %V(big0) = NaN; 
                %hLay(big0,:) = NaN; 
            end
            
            allSTs(pxi,:) = sum(V, 'omitnan'); 
        end
    end

    [maxh id] = max(abs(allSTs));
    max_clust_sum_obs(layi,:) = allSTs(id); 

    % % % % rem non-sig-clusters
    for ci = 1:length(clustinfo.PixelIdxList)
        modifiedCluster = abs(sum(tLay(clustinfo.PixelIdxList{ci})))+0.0001; %add this small number because removing negative values creates slighlly different valeus
       if modifiedCluster < max(abs(allSTs))   + 1    %add +1 at the very end to delete all clusters
          %disp(['layi > ' num2str(layi) ' > ' num2str(abs(sum(tLay(clustinfo.PixelIdxList{ci})))) '_' num2str(max(abs(allSTs)))]) 
          %hLay (clustinfo.PixelIdxList{ci}) =  0; 
       end
    end
    
    hL(layi, f, :) = hLay; 
    tL(layi, f, :) = tLay; 
end

% get 3D cluster
h3D = squeeze(hL(:,f,:)); 
t3D = squeeze(tL(:,f,:));
clear allSTs   
clustinfo = bwconncomp(h3D, 6);
if length(clustinfo.PixelIdxList) > 0
    for pxi = 1:length(clustinfo.PixelIdxList)
       allSTs(pxi) = sum(t3D(clustinfo.PixelIdxList{pxi}));% it is maxi -1 cause we are adding a zero below
    end
else
    allSTs = 0;
end
[maxh id] = max(abs(allSTs));
max_clust_sum_3D = allSTs(id); 
for ci = 1:length(clustinfo.PixelIdxList)
   modifiedCluster = abs(sum(t3D(clustinfo.PixelIdxList{ci})))+0.0001; %add this small number because removing negative values creates slighlly different valeus
   if modifiedCluster < max(abs(allSTs))  % + 1    %add +1 at the very end to delete all clusters
      %disp(['layi > ' num2str(layi) ' > ' num2str(abs(sum(tLay(clustinfo.PixelIdxList{ci})))) '_' num2str(max(abs(allSTs)))]) 
      h3D (clustinfo.PixelIdxList{ci}) =  0; 
   end
end



figure()
layT = tiledlayout(3, 3);
layT.TileSpacing = 'compact';
layT.Padding = 'compact';
for layi = 1:nLays
    
    t = squeeze(tL(layi,:,:));
    h = squeeze(hL(layi,:,:));
    
    if size(all_r_Times, 4) < 30 & size(all_r_Times, 4) > 5
        % % % % encoding 
        if nLays > 9
            subplot(7, 8, layi)
        else
            subplot(3, 3, layi)
        end
        %times = -.5:.01:1.499;
        times = 0:.01:1.499;
        %freqs = [.1:.1:29.9 30:.5:150]; 
        freqs = 1:520;
        contourf(times, freqs, myresizem(t, 10), 100, 'linecolor', 'none'); hold on; %colorbar
        contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 2);
        plot([-0.055 -0.055], get(gca, 'ylim'), 'k:', 'LineWidth', 2);
        plot(get(gca, 'xlim'),[300 300], 'k:', 'LineWidth', 2);
        set(gca, 'xtick', [-0.055 .445 .945 1.45], 'xticklabel', [0 .5 1 1.5], 'clim', [-5 5]); colormap jet
        set(gca,'XTick',[], 'YTick', [], 'xticklabel',[], 'FontSize', 12)

    elseif size(all_r_Times, 4) == 5
        % % % % encoding 
        if nLays > 9
            subplot(7, 8, layi)
        else
            subplot(3, 3, layi)
        end
        times = 0:.01:0.499;
        %freqs = [.1:.1:29.9 30:.5:150]; 
        freqs = 1:520;
        contourf(times, freqs, myresizem(t, 10), 100, 'linecolor', 'none'); hold on; %colorbar
        contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 2);
        plot([-0.055 -0.055], get(gca, 'ylim'), 'k:', 'LineWidth', 2);
        plot(get(gca, 'xlim'),[300 300], 'k:', 'LineWidth', 2);
        set(gca, 'xtick', [-0.055 .445 .945 1.45], 'xticklabel', [0 .5 1 1.5], 'clim', [-5 5]); colormap jet
        set(gca,'XTick',[], 'YTick', [], 'xticklabel',[], 'FontSize', 12)
    elseif size(all_r_Times, 4) == 546
        % % % % encoding 
        if nLays > 9
            subplot(7, 8, layi)
        else
            subplot(3, 3, layi)
        end
        times = 0:.01:5.009;
        %freqs = [.1:.1:29.9 30:.5:150]; 
        freqs = 1:52;
        contourf(times, freqs, t, 100, 'linecolor', 'none'); hold on; %colorbar
        contour(times, freqs, h, 1, 'Color', [0, 0, 0], 'LineWidth', 2);
        %plot([-0.055 -0.055], get(gca, 'ylim'), 'k:', 'LineWidth', 2);
        %plot(get(gca, 'xlim'),[300 300], 'k:', 'LineWidth', 2);
        %set(gca, 'xtick', [-0.055 .445 .945 1.45], 'xticklabel', [0 .5 1 1.5], 'clim', [-5 5]); colormap jet
        %set(gca,'XTick',[], 'YTick', [], 'xticklabel',[], 'FontSize', 12)
    else
        if nLays > 9
            subplot(7, 8, layi)
        else
            subplot(3, 3, layi)
        end
        %times = 0:.01:3.49;
       % times = 0:.01:4.49;
        times = -.5:.01:3.99;
        %freqs = [.1:.1:29.9 30:.5:150]; 
        freqs = 1:520;
        contourf(times, freqs, myresizem(t, 10), 100, 'linecolor', 'none'); hold on; 
        contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 1); %colorbar
        plot([-0.055 -0.055], get(gca, 'ylim'), 'k:', 'LineWidth', 2);
        plot(get(gca, 'xlim'),[300 300], 'k:', 'LineWidth', 2);
        set(gca, 'xtick', [-0.055 .945 1.945 2.945 3.945 ], 'xticklabel', [0 1 2 3 4], 'clim', [-3 3]); % 'xlim', [0 3.5],
        %set(gca, 'xlim', [-.5 1])
        set(gca,'XTick',[], 'YTick', [])
        set(gca,'xticklabel',[], 'FontSize', 12)
        colormap jet
    end

end
        

%f2p = [num2str(layi, '%02.f') '.png'];
%exportgraphics(gcf, f2p, 'Resolution', 300)















%% FREQUENCY RESOLVED DNN ANALYSIS PERMUTATIONS

clearvars -except act_CH act_FR 


nPerm   = 1000; 
region = 'pfc';
f2sav       = 'RNN_pfc_M_Av_P_nT_1-49_C_naFV_100ms'; %'Alex_vvs_CueAll_noAv_54_Real_nT_1-49'
paths = load_paths_WM(region);
filelistSess = getFiles(paths.out_contrasts_path);

load_parameters_WMFRA;
loadNet_WM;

nSubj = length(filelistSess); 
all_r_Times_Perm    = zeros(nPerm, nSubj, nLays, nFreqs, nTimes);

tic

for subji= 1:nSubj %start with subject 2
    disp (['Subj: ' num2str(subji)])
    load([paths.out_contrasts_path filelistSess{subji}]);   
    cfg_contrasts.oneListPow = cfg_contrasts.oneListPow(:, :, :, resTime); 
    cfg_contrasts = normalize_WM(cfg_contrasts, 1, 'blo', [3 7]);
    if strcmp(avMeth,'pow')
       cfg_contrasts = average_repetitions(cfg_contrasts);
    end
    
    for freqi = 1:length(freqs2test)
        f  = freqs2test(freqi);
        
        if strcmp(f2t{9}, 'aFV')
            avTW = 1;
        else
            avTW = 0;
        end
        rdm_prev = create_rdms(cfg_contrasts, f, it, 5, 1, avTW);%
        
       
        [allS ids act_FR2 act_CH2 cM] = getIDsAndnRDMs(rdm_prev, act_CH, act_FR);

        for permi = 1:nPerm
            idp = randperm(length(ids)); %this has to be done at each permutation but not each layer
            
            for timei = 8:37 % only locked to cue
                rdmp = squeeze(rdm_prev.rdm);
                rdm = squeeze(rdmp(:,:,timei));
                rdm(rdm ==1) = 1000;rdm(rdm ==0) = 2000;
                rdm = tril(rdm, -1);rdm(rdm==0) =[];
                rdm(rdm==1000) = 1;rdm(rdm==2000) = 0;
            
                for layi = 56:56
                    if subji < subj_ch_fr
                        M =  squeeze(act_FR2(layi,:,:)); 
                    else
                        M =  squeeze(act_CH2(layi,:,:)); 
                    end
                    M = M(idp, idp);
                    M(M ==1) = 1000;M(M ==0) = 2000;
                    M = tril(M, -1);M(M==0) =[];
                    M(M==1000) = 1;M(M==2000) = 0;
    
                    cMh = cM; 
                    cMh = cMh(idp, idp);
                    cMh(cMh ==1) = 1000;cMh(cMh ==0) = 2000;
                    cMh = tril(cMh, -1);cMh(cMh==0) =[];
                    cMh(cMh==1000) = 1;cMh(cMh==2000) = 0;

                    
                    if strcmp(f2t{8}, 'PC')
                        x = [rdm' M']; 
                        z = cMh'; 
                        allTEst_prev = partialcorr(x, z, 'type', 's');
                        allTEst = allTEst_prev(2);
                    else
                        allTEst = corr(rdm', M', 'type', 's');
                    end
                    
                    
    
                    all_r_Times_Perm(permi, subji, layi,freqi,timei) = allTEst;  
                end
            end
        end
    
    end
    
end 

mkdir ([paths.results.DNNs]);
save([paths.results.DNNs f2sav], 'all_r_Times_Perm');

toc


%% Calculate correlation in cluster 

data = squeeze(all_r_Times(:, 49,:,:));
mD = squeeze(mean(data));
[h p ci ts] = ttest(data);
h = squeeze(h); t = squeeze(ts.tstat);

imagesc(t); hold on; 
contourf(h);
clustinfo = bwconncomp(h);
obsT = sum(t(clustinfo.PixelIdxList{8}));


%% calculate distribution from permutations

nPerm = 1000; 
for permi = 1:nPerm
    %dataP = squeeze(all_r_Times_Perm(permi, :, 56,:,8:17));
    dataP = squeeze(all_r_Times_Perm(permi, :, 56,:,8:37));
    [h p ci ts] = ttest(dataP);
    h = squeeze(h); t = squeeze(ts.tstat);
    
    clear allSTs  
    clustinfo = bwconncomp(h);
    if length(clustinfo.PixelIdxList) > 0
        for pxi = 1:length(clustinfo.PixelIdxList)
           allSTs(permi, pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
        end
    else
        allSTs(permi, :) = 0;
    end

    [maxh id] = max(abs(allSTs(permi)));
    max_clust_sum_perm(permi,:) = allSTs(permi,id); 


end




%% Calculate correlation in cluster for each layer / time point

data = squeeze(all_r_Times(:, 49,:,:));
mD = squeeze(mean(data));
[h p ci ts] = ttest(data);
h = squeeze(h); t = squeeze(ts.tstat);

imagesc(t); hold on; 
contourf(h);
clustinfo = bwconncomp(h);
obsT = sum(t(clustinfo.PixelIdxList{8}));


%%
c = [60.7767903066970 69.7513849905440 56.9934561621767 72.4082145601195 79.5953647172531 78.9209984633255 90.2397712325105 90.4385912362766];
c = [24 26 22 28 29 28 31 31]
figure()
plot(c)
set(gca, 'yLim', [20 35])




%% TRIAL BASED DNN ANALYSIS SLOW AS SANITY CHECK

clearvars -except act_CH act_FR 
f2sav       = 'RNN_pfc_M_noAv_R_T_1-20_C_naFV_100ms'; %'Alex_vvs_CueAll_noAv_54_Real_nT_1-49'

load_parameters_WMFRA;
loadNet_WM;
region = f2t{2};
paths = load_paths_WM(region);
filelistSess = getFiles(paths.out_contrasts_path);
nSubj = length(filelistSess); 
all_r_Times    = zeros(nSubj, nLays, nFreqs, nTimes);

tic

for subji= 1:nSubj %start with subject 2
    disp (['Subj: ' num2str(subji)])
    load([paths.out_contrasts_path filelistSess{subji}]);   
    cfg_contrasts.oneListPow = cfg_contrasts.oneListPow(:, :, :, resTime); 
    cfg_contrasts = normalize_WM(cfg_contrasts, 1, 'blo', [3 7]);
    if strcmp(avMeth,'pow')
       cfg_contrasts = average_repetitions(cfg_contrasts);
    end
    
    
    
    for freqi = 1:length(freqs2test)
        f  = freqs2test(freqi);
        if strcmp(f2t{9}, 'aFV')
            avTW = 1;
        else
            avTW = 0;
        end
        rdm_prev(freqi) = create_rdms(cfg_contrasts, f, it, 5, 1, avTW);
    end
    
    ids = rdm_prev(1).ids; % take first freq only to get ids
    idsT{subji,:} = ids; 
    ids = cellfun(@(x) strsplit(string(x)), ids, 'un', 0);
    ids0 = cellfun(@(x) x(3), ids, 'un', 0); ids0 = double(string(ids0)); 
    ids1 = char(string((ids0))); ids2 = ids1(:,[1 3]);ids3 = double(string(ids2));
    idx = find(~mod(ids3, 10)); ids4 = ids3-10; ids4(idx) = ids3(idx);

    act_FR_2 = act_FR(:,ids4,ids4);
    act_CH_2 = act_CH(:,ids4,ids4);
    
    allS = cat(4, rdm_prev.rdm);
    % no need to remove diagonal because diagonal falls in index triali and is deleted below 
    
    for layi = 1:nLays
        for triali = 1:size(allS, 1)
            if subji < subj_ch_fr
                M =  squeeze(act_FR_2(layi,:,:)); 
            else
                M =  squeeze(act_CH_2(layi,:,:)); 
            end
            M = M(triali, :);
            M(triali) = []; 
            
            for freqi = 1:size(allS, 4)  
                for timei = 1:size(allS, 3)
                    allS1 = squeeze(allS(triali, :, timei, freqi)); 
                    allS1(triali) = []; %remove same trial on diagonal (only has ones)

                    allTEst = corr(allS1', M', 'type', 's');
                    all_r_Times_Trials{subji}(layi,triali, freqi, timei) = allTEst;  
                end
            end
        end
    end
end 


toc 
save(f2sav, 'all_r_Times_Trials', 'idsT');

%% PFC cue performance effect BEHAVIOR

clear
load RNN_pfc_M_Av_R_nT_1-20_C_naFV_100ms
load RNN_pfc_M_noAv_R_T_1-20_C_naFV_100ms


sub2exc =  [];

nSubjs =size(all_r_Times, 1);  
nFreqs = size(all_r_Times, 3); 
clear hL tL clustinfo
allR = all_r_Times;
allR(sub2exc,:,:) = []; 
for layi = 1:56
    [hL(layi, :, :) p ci ts] = ttest(allR(:,layi,:,:)); 
    tL(layi, :, :) = ts.tstat;
end


%%get cluster in each layer 
clear max_clust_sum_obs allSTs clustinfo_lay
for layi = 1:56
    hLay = squeeze(hL(layi,:,:)); 
    tLay = squeeze(tL(layi,:,:));
    clustinfo_lay(layi) = bwconncomp(hLay);
end
    

%%extract mean rho values in each trial at this cluster for 1 layer 
clear mARTT 
for subji = 1:16
   for layi = 56:56
       all_r_tt = squeeze(all_r_Times_Trials{subji}(layi, :,:,:)); 
       nTrials = size(all_r_tt, 1); 
        for triali = 1:nTrials
            all_r_ttt = squeeze(all_r_tt(triali,:,:)); 
            mARTT{subji,:}(triali,:) = mean(all_r_ttt(clustinfo_lay(56).PixelIdxList{8}), 'omitnan'); 
        end   
   end
end

%% CUEPERF EFFECT in one layer

clear cR xC xI perf_cc
for subji = 1:16
   
    clear ic_i
    x = mARTT{subji};
    
    
    ids = idsT{subji};
    ids0 = cellfun(@(x) strsplit(string(x)), ids, 'un', 0);
    % 19 = correct item level ; 20 = correct category level
    ic_i = cellfun(@(x) x(19), ids0, 'un', 0); ic_i = logical(double(string(ic_i))); 
    
    xC{subji,:}= x(ic_i); 
    xI{subji,:}= x(~ic_i); 
    
end


for subji = 1:16   
    perf_cc(subji, 1) = mean(xC{subji});
    perf_cc(subji, 2) = mean(xI{subji});
end



%% 2 bar

data.data = perf_cc; 

sub2exc = [1];

data.data(sub2exc,:) =  [];

figure(2); set(gcf,'Position', [0 0 500 500]); 
mean_S = mean(data.data, 1, 'omitnan');
h = bar (mean_S);hold on;
hb = plot ([1:2], data.data); hold on;
set(hb, 'lineWidth', 3, 'Marker', '.', 'MarkerSize',40);hold on;
%h1=errorbar(mean_S, se_S,'c'); hold on;
set(h,'facecolor', 'none', 'lineWidth', 2);
%set(hb, 'Color','k','linestyle','none', 'lineWidth', 2);
set(gca,'XTick',[1:2],'XTickLabel',{'', ''}, 'FontSize', 30, 'linew',1, 'xlim', [0 3], 'ylim', [-.15 .15] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 1);
%scatter(1:2, hL+0.135, 'k', 'LineWidth', 6, 'Marker', '*')


[h p ci t] = ttest (data.data(:,1), data.data(:,2));
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth', 1);

%export_fig(2, [num2str(layi) '.png'],'-transparent', '-r100');



%% DNN ANALYSIS BANDS


clearvars -except act_CH act_FR 

f = 13:29; 
region = 'pfc';
f2sav       = 'Alex_pfc_M_Av_R_nT_1-49_C_aFV_100ms_B'; %'Alex_vvs_CueAll_noAv_54_Real_nT_1-49'
paths = load_paths_WM(region);
filelistSess = getFiles(paths.out_contrasts_path);

load_parameters_WMFRA;
loadNet_WM;

nSubj = length(filelistSess); 
all_r_Times    = zeros(nSubj, nLays, nTimes);

tic

for subji= 1:nSubj %start with subject 2
    disp (['Subj: ' num2str(subji)])
    load([paths.out_contrasts_path filelistSess{subji}]);   
    cfg_contrasts.oneListPow = cfg_contrasts.oneListPow(:, :, :, resTime); 
    cfg_contrasts = normalize_WM(cfg_contrasts, 1, 'blo', [3 7]);
    if strcmp(avMeth,'pow')
       cfg_contrasts = average_repetitions(cfg_contrasts);
    end
    

    if strcmp(f2t{9}, 'aFV')
        avTW = 1;
    else
        avTW = 0;
    end
    rdm_prev = create_rdms(cfg_contrasts, f, it, 5, 1, avTW);%
    
   
    [allS ids act_FR2 act_CH2 cM] = getIDsAndnRDMs(rdm_prev, act_CH, act_FR);
    
    for timei = 1:45
        rdmp = squeeze(rdm_prev.rdm);
        rdm = squeeze(rdmp(:,:,timei));
        rdm(rdm ==1) = 1000;rdm(rdm ==0) = 2000;
        rdm = tril(rdm, -1);rdm(rdm==0) =[];
        rdm(rdm==1000) = 1;rdm(rdm==2000) = 0;
    
        for layi = 1:7
            if subji < subj_ch_fr
                M =  squeeze(act_FR2(layi,:,:)); 
            else
                M =  squeeze(act_CH2(layi,:,:)); 
            end
            M(M ==1) = 1000;M(M ==0) = 2000;
            M = tril(M, -1);M(M==0) =[];
            M(M==1000) = 1;M(M==2000) = 0;

            cMh = cM; 
            cMh(cMh ==1) = 1000;cMh(cMh ==0) = 2000;
            cMh = tril(cMh, -1);cMh(cMh==0) =[];
            cMh(cMh==1000) = 1;cMh(cMh==2000) = 0;

            
            if strcmp(f2t{8}, 'PC')
                x = [rdm' M']; 
                z = cMh'; 
                allTEst_prev = partialcorr(x, z, 'type', 's');
                allTEst = allTEst_prev(2);
            else
                allTEst = corr(rdm', M', 'type', 's');
            end
            
            

            all_r_Times(subji, layi,timei) = allTEst;  
        end
    end
    
    
end 

mkdir ([paths.results.DNNs]);
save([paths.results.DNNs f2sav], 'all_r_Times');

toc





%% plot observed correlations per layer (BANDS)

subj2exc = [1]; 

cfg.yLim = [-.015 .035]; 
cfg.sP = 1;
cfg.Fz = 8; 
cfg.Lw = 6;
cfg.PsL = 1; 
nLays = size(all_r_Times, 2);

figure()
for layi = 1:nLays
    if nLays < 9
        subplot(3, 3, layi); 
    else
        subplot(7,8, layi); 
    end
    r_Time = squeeze(all_r_Times(:, layi, :,:)); 
    r_Time(subj2exc, : ) = []; 

    [max_clust_sum_obs(layi)] = plot_with_sigV(r_Time, cfg); hold on; 

   
end

 exportgraphics(gcf, 'figure.png', 'Resolution', 300)
 
 
%% average full maint period per layer 

subj2exc = []; 
nLays = size(all_r_Times, 2); 
nTimes = size(all_r_Times, 4); 
if nTimes > 30
    all_r_TM = squeeze(mean(all_r_Times(:, :, :, 1:35), 4));
elseif nTimes <= 10
    all_r_TM = squeeze(mean(all_r_Times, 4));
else
    all_r_TM = squeeze(mean(all_r_Times(:, :, :, 6:13), 4));
end

figure()
mART = mean(all_r_TM); 
stdART = std(all_r_TM); 
seART = stdART/ sqrt(size(all_r_TM, 1));
[h p ci t] = ttest(all_r_TM); 
hL = h; hL(hL == 0) = nan; hL(hL==1) = 0; 

shadedErrorBar(1:nLays, mART, seART, 'r', 1); hold on; 
plot(1:nLays, hL, 'LineWidth', 6)
%plot (allRTMA); hold on; 
set(gca, 'FontSize', 12)



exportgraphics(gcf, 'figure.png', 'Resolution', 300)
 


%% Bar format
data.data = [all_r_TM]; 

figure(2); set(gcf,'Position', [0 0 1200 1000]); 
mean_S = mean(data.data, 1);
hb = plot ([1:nLays], data.data); hold on;
set(hb, 'lineWidth', 3, 'Marker', '.', 'MarkerSize',10);hold on;
h = bar (mean_S);hold on;
%h1=errorbar(mean_S, se_S,'c'); hold on;
set(h,'facecolor', 'flat', 'lineWidth', 1);
set(hb, 'Color','k','linestyle','none', 'lineWidth', 2);
set(gca,'XTick',[1:nLays],'XTickLabel',{'', ''}, ...
    'FontSize', 20, 'linew',2, 'xlim', [0 nLays+1], 'ylim', [-.02 .04] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 1);
scatter(1:nLays, hL+0.025, 'k', 'LineWidth', 6, 'Marker', '*')

if nLays == 56
    c1 = [ones(1, 8) ones(1, 8)*2 ones(1, 8)*3 ones(1, 8)*4 ones(1, 8)*5 ones(1, 8)*6 ones(1, 8)*7 ones(1, 8)*8]' /7; % sorted by layer
    c1 = repmat((1:8)/8, 1, 7)'; % sorted by time point
    c2 = repmat(zeros(1,2), nLays, 1);
    cols = [c1 c2];
    h.CData = cols;
else
    c1 = (1:8)' / 8 ; 
    c2 = repmat(zeros(1,2), nLays, 1);
    cols = [c1 c2];
    h.CData = cols; 
end



[h p ci t] = ttest (data.data);
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth', 1);

export_fig(2, '_2.png','-transparent', '-r300');
close all;   


 

%% 56 Bar 
data.data = [all_r_TM]; 

figure(2); set(gcf,'Position', [0 0 1200 1000]); 
mean_S = mean(data.data, 1);
hb = plot ([1:56], data.data); hold on;
set(hb, 'lineWidth', 3, 'Marker', '.', 'MarkerSize',10);hold on;
h = bar (mean_S);hold on;
%h1=errorbar(mean_S, se_S,'c'); hold on;
set(h,'facecolor', 'flat', 'lineWidth', 1);
set(hb, 'Color','k','linestyle','none', 'lineWidth', 2);
set(gca,'XTick',[1:56],'XTickLabel',{'', ''}, ...
    'FontSize', 20, 'linew',2, 'xlim', [0 57], 'ylim', [-.045 .075] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 1);
scatter(1:56, hL+0.07, 'k', 'LineWidth', 6, 'Marker', '*')

if nLays == 56
    c1 = [ones(1, 8) ones(1, 8)*2 ones(1, 8)*3 ones(1, 8)*4 ones(1, 8)*5 ones(1, 8)*6 ones(1, 8)*7 ones(1, 8)*8]' /7; % sorted by layer
else
    c1 = (1:8)' / 8 ; 
end

c1 = repmat((1:8)/8, 1, 7)'; % sorted by time point
%c2 = repmat(zeros(1,2), nLays, 1);
cols = [c1 c2];
h.CData = cols;

[h p ci t] = ttest (data.data);
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth', 1);

export_fig(2, '_2.png','-transparent', '-r300');
close all;   


 
 
%% MEAN FOR EACH LAYER
nLays = size(all_r_Times, 2); 
subj2exc = []; 
r_Time = squeeze(mean(all_r_Times, 1)); 

% % % sort by layer
%c1 = (0:1/nLays:0.99)'; % layer and timepoint 
if nLays == 56
    c1 = [ones(1, 8) ones(1, 8)*2 ones(1, 8)*3 ones(1, 8)*4 ones(1, 8)*5 ones(1, 8)*6 ones(1, 8)*7 ones(1, 8)*8]' /7; % sorted by layer
else
    c1 = (1:8)' / 8 ; 
end

%c1 = repmat((1:8)/8, 1, 7)'; % sorted by time point
c2 = repmat(zeros(1,2), nLays, 1);
cols = [c1 c2];
colororder(cols)
plot(r_Time', 'linewidth', 2)
exportgraphics(gcf, 'figure.png', 'Resolution', 300)


 

%% check correlations with different variance

x = 5* randn(1, 1000) + 10; 
y = 5* randn(1, 1000) + 10; 

z = 5* randn(1, 1000) + 10; 
u = 40* randn(1, 1000) + 10; 

figure()
histogram(x); hold on; 
histogram(y)

figure()
histogram(u); hold on; 
histogram(z)

rho1 = corr(x', y', 'type', 's')
rho2 = corr(z', u', 'type', 's')


%% TRIAL BASED DNN ANALYSIS SLOW AS SANITY CHECK

clearvars -except act_CH act_FR 
f2sav       = 'RNN_hipp_M_noAv_R_T_1-49_C_aFV_100ms'; %'Alex_vvs_CueAll_noAv_54_Real_nT_1-49'
%f2sav       = 'RNN_pfc_M&C_noAv_54_Real_T_1-49_C_aFV';
f           = 13:29; %in case B
load_parameters_WMFRA;
loadNet_WM;

all_r_Times    = zeros(nSubj, nLays, nFreqs, nTimes);
%lois = (1:8:56) +7; 
lois = 56; 
f2t = strsplit(f2sav, '_'); 

tic
for subji=1:nSubj
    disp (['Subj: ' num2str(subji)])
    load(sublist{subji});
    cfg_contrasts.oneListPow = cfg_contrasts.oneListPow(:, :, :, resTime); 
    cfg_contrasts = normalize_WM(cfg_contrasts, 1, 'blo', [3 7]);
    if strcmp(avMeth,'pow')
       cfg_contrasts = average_repetitions(cfg_contrasts);
    end
    
    
    
    for freqi = 1:length(freqs2test)
        f  = freqs2test(freqi);
        if strcmp(f2t{9}, 'aFV')
            avTW = 1;
        else
            avTW = 0;
        end
        rdm_prev(freqi) = create_rdms(cfg_contrasts, f, it, 5, 1, avTW);
    end
    
    ids = rdm_prev(1).ids; % take first freq only to get ids
    idsT{subji,:} = ids; 
    ids = cellfun(@(x) strsplit(string(x)), ids, 'un', 0);
    ids0 = cellfun(@(x) x(3), ids, 'un', 0); ids0 = double(string(ids0)); 
    ids1 = char(string((ids0))); ids2 = ids1(:,[1 3]);ids3 = double(string(ids2));
    idx = find(~mod(ids3, 10)); ids4 = ids3-10; ids4(idx) = ids3(idx);

    act_FR_2 = act_FR(:,ids4,ids4);
    act_CH_2 = act_CH(:,ids4,ids4);
    
    allS = cat(4, rdm_prev.rdm);
    % no need to remove diagonal because diagonal falls in index triali and is deleted below 
    
    for layi = 1:length(lois) %nLays
        for triali = 1:size(allS, 1)
            if subji < subj_ch_fr
                M =  squeeze(act_FR_2(lois(layi),:,:)); 
            else
                M =  squeeze(act_CH_2(lois(layi),:,:)); 
            end
            M = M(triali, :);
            M(triali) = []; 
            
            for freqi = 1:size(allS, 4)  
                for timei = 1:size(allS, 3)      
                    allS1 = squeeze(allS(triali, :, timei, freqi)); 
                    allS1(triali) = []; %remove same trial on diagonal (only has ones)

                    allTEst = corr(allS1', M', 'type', 's');
                    all_r_Times_Trials{subji}(layi,triali, freqi, timei) = allTEst;  
                end
            end
        end
    end
end 

cd (currentF) 
toc 
save(f2sav, 'all_r_Times_Trials', 'idsT');



%% TRIAL BASED DNN ANALYSIS SLOW AS SANITY CHECK 2

clearvars -except act_CH act_FR 
f2sav       = 'RNN_pfc_M_noAv_R_T_1-49_C_aFV_100ms'; %'Alex_vvs_CueAll_noAv_54_Real_nT_1-49'
%f2sav       = 'RNN_pfc_M&C_noAv_54_Real_T_1-49_C_aFV';
f           = 13:29; %in case B
load_parameters_WMFRA;
loadNet_WM;

all_r_Times    = zeros(nSubj, nLays, nFreqs, nTimes);
%lois = (1:8:56) +7; 
lois = 56; 
f2t = strsplit(f2sav, '_'); 

tic
for subji=1:nSubj
    disp (['Subj: ' num2str(subji)])
    load(sublist{subji});
    cfg_contrasts.oneListPow = cfg_contrasts.oneListPow(:, :, :, resTime); 
    cfg_contrasts = normalize_WM(cfg_contrasts, 1, 'blo', [3 7]);
    if strcmp(avMeth,'pow')
       cfg_contrasts = average_repetitions(cfg_contrasts);
    end
    
    %build rdm_prev for each subject
    for freqi = 1:length(freqs2test)
        f  = freqs2test(freqi);
        if strcmp(f2t{9}, 'aFV')
            avTW = 1;
        else
            avTW = 0;
        end
        rdm_prev = create_rdms(cfg_contrasts, f, it, 5, 1, avTW);
    
    
        ids = rdm_prev.ids; % take first freq only to get ids
        idsT{subji,:} = ids; 
        ids = cellfun(@(x) strsplit(string(x)), ids, 'un', 0);
        ids0 = cellfun(@(x) x(3), ids, 'un', 0); ids0 = double(string(ids0)); 
        ids1 = char(string((ids0))); ids2 = ids1(:,[1 3]);ids3 = double(string(ids2));
        idx = find(~mod(ids3, 10)); ids4 = ids3-10; ids4(idx) = ids3(idx);

        act_FR_2 = act_FR(:,ids4,ids4);
        act_CH_2 = act_CH(:,ids4,ids4);
    
    
        for layi = 1:length(lois) %nLays
            for triali = 1:size(rdm_prev.rdm, 1)
                if subji < subj_ch_fr
                    M =  squeeze(act_FR_2(lois(layi),:,:)); 
                else
                    M =  squeeze(act_CH_2(lois(layi),:,:)); 
                end
                M = M(triali, :);
                M(triali) = []; 
                for timei = 1:size(rdm_prev.rdm, 3)      
                    allS = squeeze(rdm_prev.rdm); 
                    allS1 = allS(triali, :, timei); 
                    allS1(triali) = []; %remove same trial on diagonal (only has ones)

                    allTEst = corr(allS1', M', 'type', 's');
                    all_r_Times_Trials{subji}(layi,triali, freqi, timei) = allTEst;  
                end
            end
        end
    end
end 

cd (currentF) 
toc 
save(f2sav, 'all_r_Times_Trials', 'idsT');



%% plot 
%load RNN_pfc_M_noAv_R_T_1-49_C_aFV_100ms_5-1win
%load RNN_pfc_M_noAv_R_T_1-550_C_aFV_10ms_1-1win
size(all_r_Times_Trials{1})
pfc = cellfun(@(x) squeeze(mean(x,2)), all_r_Times_Trials, 'un', 0)';

%pfc2 = vertcat(pfc{:});
pfc2 = [pfc{:}];
pfc3 = reshape(pfc2, [], 52, 16);
pfc3 = permute(pfc3, [3 2 1]);

d2p = squeeze(pfc3(9, :, :));  %check that reshape is correct (n1 = all %nans)
%d2p = squeeze(mean(pfc3, 'omitnan'));
figure()
imagesc(d2p);colorbar


%% 
sub2exc = [4]
pfc4 = pfc3; 
pfc4(subj2exc, :, : ) = []; 
[h p ci t] = ttest(pfc4)
h = squeeze(h); t = squeeze(t.tstat)
imagesc(h)


%% plot steps

%d2p = squeeze( cfg_contrasts.oneListPow(1,1,:,:) );
d2p = squeeze( rdm_prev(1).rdm(:,:,1) );
figure()
imagesc(d2p); axis square




%% extract PFC cluster
clear
load RNN_pfc_M&C_noAv_54_Real_T_1-49_C_aFV
load RNN_vvs_M&C_noAv_54_Real_T_1-49_C_aFV
load RNN_pfc_M&C_noAv_54_Real_nT_1-49_C_aFV

lois = (1:8:56) +7; 

sub2exc =  [];

nSubjs =size(all_r_Times, 1);  
nFreqs = size(all_r_Times, 3); 
clear hL tL clustinfo
for freqi = 1:nFreqs
    allR = squeeze(all_r_Times(:, :, freqi, :));
    allR(sub2exc,:,:) = []; 
    for layi = 1:7
        [hL(layi, freqi, :) p ci ts] = ttest(allR(:,lois(layi),:)); 
        tL(layi, freqi, :) = ts.tstat;
    end
end

%%get cluster in each layer 
clear max_clust_sum_obs allSTs clustinfo_lay
for layi = 1:7
    
    hLay = squeeze(hL(layi,:,:)); 
    tLay = squeeze(tL(layi,:,:));
    
    clear allSTs  
    clustinfo_lay(layi) = bwconncomp(hLay);
    
end
    
%%extract mean rho values in each trial at this cluster for 1 layer 
apix = [1 2 3 4 5 7 10]; %8th time-point
%apix = [5 8 6 7 4 3 4]; %7th time-point


%% build 2 compare

clear 
load pfc_elec
load RNN_pfc_M&C_noAv_54_Real_T_1-49_C_aFV
pfc = all_r_Times_Trials([2 3  5  9 10 11 12 14 15 16])';
load RNN_vvs_M&C_noAv_54_Real_T_1-49_C_aFV
vvs = all_r_Times_Trials([7 9 13 18 19 20 21 23 27 28])';


tmp = cell2mat(cellfun(@size, chanNames_all, 'un', 0));
s2e_pfc = ~(tmp([2 3  5  9 10 11 12 14 15 16],1) > 1); 
%sub2exc = []; 
pfc(s2e_pfc) = [];vvs(s2e_pfc) = [];

ff1 = 18:28; %freq starts at 3
tt1 = 7:12; 
% ff2 = 1:6; 
% tt2 = 16:40; 

% % % use clustinfo
ff2 = 1:6; % 1: 6 is low theta 
tt2 = 16:40; 


nLays = 7;
clear cR

for subji = 1:length(pfc)
   
    x = pfc{subji};
    y = vvs{subji};
    
    figure(); set(gcf, 'Position', [100 100 1000 1000]); ploti = 1;
    for layi = 1:nLays % 7:7 %
        x1 = squeeze(mean(mean(x(layi, :, ff1, tt1), 4, 'omitnan'),3, 'omitnan'))'; 
        %x1p = squeeze(x(layi, :, :,:)); 
        %x1 = mean(x1p(:, clustinfo_lay(7).PixelIdxList{10}),2);
        for layj = 1:nLays %7:7 %
            y1 = squeeze(mean(mean(y(layj, :, ff2, tt2), 4, 'omitnan'),3, 'omitnan'))'; %16:35 max
            %y1p = squeeze(y(layj, :, :,:)); 
            %y1 = mean(y1p(:, clustinfo_lay(7).PixelIdxList{10}),2);
            
            cR(subji, layi, layj, :) = corr(x1, y1, 'type', 's');
            
            subplot(7 ,7, ploti)
            %scatterHistDiff(x1, y1);
            scatter(x1, y1, 50, '.'); hold on; 
            p = polyfit(x1,y1,1); 
            f = polyval(p,x1); 
            plot(x1,y1,'.',x1,f,'-', 'LineWidth', 2) 
            
            ploti = ploti+1;
            
            
        end
    end    
    
    exportgraphics(gcf, [num2str(subji) '.png'], 'Resolution', 300);
    close all
end

[h p ci ts] = ttest(cR, 0, 'Alpha', 0.05);
h = squeeze(h); t= squeeze(ts.tstat);


clear allSTs   
clustinfo = bwconncomp(h);
if length(clustinfo.PixelIdxList) > 0
    for pxi = 1:length(clustinfo.PixelIdxList)
       allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% it is maxi -1 cause we are adding a zero below
    end
else
    allSTs= 0;
end

[maxh id] = max(abs(allSTs));
max_clust_sum_real = allSTs(id); 



figure()
contourf(myresizem(t, 20), 40, 'linecolor', 'none'); axis equal, hold on; colorbar
contour(myresizem(h, 20), 1, 'Color', [0, 0, 0], 'LineWidth', 2);
set(gca, 'xtick', [10 30 50 70 90 110 130], 'xticklabels', {[1:7]})
set(gca, 'ytick', [10 30 50 70 90 110 130], 'yticklabels', {[1:7]})
set(gca, 'FontSize', 22) %, 'clim', [-4 4]
%set(gca, 'xlim',  [0.5 nLays+0.5], 'ylim', [.5  nLays+0.5]) %, 'clim', [0 180]


exportgraphics(gcf, 'trial_based.png', 'Resolution', 300)






%% build 2 compare hippocampus VVS

clear 
load RNN_hipp_M&C_noAv_54_Real_T_1-49_C_aFV
pfc = all_r_Times_Trials([1:16])';
load RNN_vvs_M&C_noAv_54_Real_T_1-49_C_aFV
vvs = all_r_Times_Trials([1 5 10 12:15 17 18 21 23:28])';
load hipp_elec

tmp = cell2mat(cellfun(@size, chanNames_all, 'un', 0));
s2e_pfc = ~(tmp([1:16],1) > 0); 
%sub2exc = []; 
pfc(s2e_pfc) = [];vvs(s2e_pfc) = [];



ff1 = 28:36; %freq starts at 3
tt1 = 11:40; 
% ff2 = 1:6; 
% tt2 = 16:40; 

% % % use clustinfo
ff2 = 1:6; % 1: 6 is low theta 
tt2 = 6:10; 


nLays = 7;
clear cR

for subji = 1:length(pfc)
   
    x = pfc{subji};
    y = vvs{subji};
    
    %figure(); set(gcf, 'Position', [100 100 1000 1000]); ploti = 1;
    for layi = 1:nLays % 7:7 %
        x1 = squeeze(mean(mean(x(layi, :, ff1, tt1), 4, 'omitnan'),3, 'omitnan'))'; 
        %x1p = squeeze(x(layi, :, :,:)); 
        %x1 = mean(x1p(:, clustinfo_lay(7).PixelIdxList{10}),2);
        for layj = 1:nLays %7:7 %
            y1 = squeeze(mean(mean(y(layj, :, ff2, tt2), 4, 'omitnan'),3, 'omitnan'))'; %16:35 max
            %y1p = squeeze(y(layj, :, :,:)); 
            %y1 = mean(y1p(:, clustinfo_lay(7).PixelIdxList{10}),2);
            
            cR(subji, layi, layj, :) = corr(x1, y1, 'type', 's');
            
            %subplot(7 ,7, ploti)
            %scatter(x1, y1, 50, '.'); hold on; 
            %p = polyfit(x1,y1,1); 
            %f = polyval(p,x1); 
            %plot(x1,y1,'.',x1,f,'-', 'LineWidth', 2) 
            %ploti = ploti+1;
            
            
        end
    end    
    
    %exportgraphics(gcf, [num2str(subji) '.png'], 'Resolution', 300);
    close all
end

[h p ci ts] = ttest(cR, 0, 'Alpha', 0.05);
h = squeeze(h); t= squeeze(ts.tstat);


clear allSTs   
clustinfo = bwconncomp(h);
if length(clustinfo.PixelIdxList) > 0
    for pxi = 1:length(clustinfo.PixelIdxList)
       allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% it is maxi -1 cause we are adding a zero below
    end
else
    allSTs= 0;
end

[maxh id] = max(abs(allSTs));
max_clust_sum_real = allSTs(id); 



figure()
contourf(myresizem(t, 20), 40, 'linecolor', 'none'); axis equal, hold on; colorbar
contour(myresizem(h, 20), 1, 'Color', [0, 0, 0], 'LineWidth', 2);
set(gca, 'xtick', [10 30 50 70 90 110 130], 'xticklabels', {[1:7]})
set(gca, 'ytick', [10 30 50 70 90 110 130], 'yticklabels', {[1:7]})
set(gca, 'FontSize', 22) %, 'clim', [-4 4]
%set(gca, 'xlim',  [0.5 nLays+0.5], 'ylim', [.5  nLays+0.5]) %, 'clim', [0 180]


exportgraphics(gcf, 'trial_based.png', 'Resolution', 300)


%% build 2 compare PFC Hippocampus

clear 
load RNN_pfc_M&C_noAv_54_Real_T_1-49_C_aFV
pfc = all_r_Times_Trials([5 9 12 14:16])';
load RNN_hipp_M&C_noAv_54_Real_T_1-49_C_aFV
hipp = all_r_Times_Trials([5 9 10 11 15 16])';
load hipp_elec

%cellfun(@length, pfc)
%cellfun(@length, hipp)


% tmp = cell2mat(cellfun(@size, chanNames_all, 'un', 0));
% s2e_pfc = ~(tmp([1:16],1) > 1); 
% %sub2exc = []; 
% hipp(s2e_pfc) = [];pfc(s2e_pfc) = [];



ff1 = 16:26; %freq starts at 3
tt1 = 6:40; 

ff2 = 1:6; % 1: 6 is low theta 
tt2 = 6:40; 


nLays = 7;
clear cR

for subji = 1:length(pfc)
   
    x = pfc{subji};
    y = hipp{subji};
    
    %figure(); set(gcf, 'Position', [100 100 1000 1000]); ploti = 1;
    for layi = 1:nLays % 7:7 %
        x1 = squeeze(mean(mean(x(layi, :, ff1, tt1), 4, 'omitnan'),3, 'omitnan'))'; 
        %x1p = squeeze(x(layi, :, :,:)); 
        %x1 = mean(x1p(:, clustinfo_lay(7).PixelIdxList{10}),2);
        for layj = 1:nLays %7:7 %
            y1 = squeeze(mean(mean(y(layj, :, ff2, tt2), 4, 'omitnan'),3, 'omitnan'))'; %16:35 max
            %y1p = squeeze(y(layj, :, :,:)); 
            %y1 = mean(y1p(:, clustinfo_lay(7).PixelIdxList{10}),2);
            
            cR(subji, layi, layj, :) = corr(x1, y1, 'type', 's');
            
            %subplot(7 ,7, ploti)
            %scatter(x1, y1, 50, '.'); hold on; 
            %p = polyfit(x1,y1,1); 
            %f = polyval(p,x1); 
            %plot(x1,y1,'.',x1,f,'-', 'LineWidth', 2) 
            %ploti = ploti+1;
            
            
        end
    end    
    
    %exportgraphics(gcf, [num2str(subji) '.png'], 'Resolution', 300);
    close all
end

[h p ci ts] = ttest(cR, 0, 'Alpha', 0.05);
h = squeeze(h); t= squeeze(ts.tstat);


clear allSTs   
clustinfo = bwconncomp(h);
if length(clustinfo.PixelIdxList) > 0
    for pxi = 1:length(clustinfo.PixelIdxList)
       allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% it is maxi -1 cause we are adding a zero below
    end
else
    allSTs= 0;
end

[maxh id] = max(abs(allSTs));
max_clust_sum_real = allSTs(id); 



figure()
contourf(myresizem(t, 20), 40, 'linecolor', 'none'); axis equal, hold on; colorbar
contour(myresizem(h, 20), 1, 'Color', [0, 0, 0], 'LineWidth', 2);
set(gca, 'xtick', [10 30 50 70 90 110 130], 'xticklabels', {[1:7]})
set(gca, 'ytick', [10 30 50 70 90 110 130], 'yticklabels', {[1:7]})
set(gca, 'FontSize', 22) %, 'clim', [-4 4]
%set(gca, 'xlim',  [0.5 nLays+0.5], 'ylim', [.5  nLays+0.5]) %, 'clim', [0 180]


exportgraphics(gcf, 'trial_based.png', 'Resolution', 300)





%% build 2 compare within region
load RNN_pfc_M&C_noAv_54_Real_T_1-49_C_aFV
pfc = all_r_Times_Trials';
%vvs = all_r_Times_Trials_VVS([7 9 13 18 19 20 21 23 27 28])';


ff1 = 16:26; %16:26
ff2 = 16:26; 

tt1 = 7:11; 
tt2 = 16:40; 



nLays = length(lois);
clear cR
for subji = 1:length(pfc)
   
    x = pfc{subji};
    y = pfc{subji};
    
    for layi = 1:nLays % 7:7 %
        x1 = squeeze(mean(mean(x(layi, :, ff1, tt1), 4, 'omitnan'),3, 'omitnan'))'; 
        %x1p = squeeze(x(lois(layi), :, :,:)); 
        %x1 = mean(x1p(:, clustinfo_lay(7).PixelIdxList{10}),2);
        for layj = 1:nLays %7:7 %
            y1 = squeeze(mean(mean(y(layj, :, ff2, tt2), 4, 'omitnan'),3, 'omitnan'))'; %16:35 max
            %y1p = squeeze(y(lois(layj), :, :,:)); 
            %y1 = mean(y1p(:, clustinfo_lay(7).PixelIdxList{10}),2);
            
            cR(subji, layi, layj, :) = corr(x1, y1, 'type', 's');
            
        end
    end    
    
end

[h p ci ts] = ttest(cR, 0, 'Alpha', 0.05);
h = squeeze(h); t= squeeze(ts.tstat); 


clear allSTs   
clustinfo = bwconncomp(h);
if length(clustinfo.PixelIdxList) > 0
    for pxi = 1:length(clustinfo.PixelIdxList)
       allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% it is maxi -1 cause we are adding a zero below
    end
else
    allSTs= 0;
end

[maxh id] = max(abs(allSTs));
max_clust_sum_real = allSTs(id); 

  

figure()
contourf(myresizem(t, 20), 40, 'linecolor', 'none'); axis equal, hold on; colorbar
contour(myresizem(h, 20), 1, 'Color', [0, 0, 0], 'LineWidth', 2);
set(gca, 'xtick', [10 30 50 70 90 110 130], 'xticklabels', {[1:7]})
set(gca, 'ytick', [10 30 50 70 90 110 130], 'yticklabels', {[1:7]})
set(gca, 'FontSize', 22) %, 'clim', [-4 4]
%set(gca, 'xlim',  [0.5 nLays+0.5], 'ylim', [.5  nLays+0.5]) %, 'clim', [0 180]


exportgraphics(gcf, 'trial_based.png', 'Resolution', 300)






%% one bar

data.data = cR(:,2,6); 

figure(2); set(gcf,'Position', [0 0 500 500]); 
mean_S = mean(data.data, 1, 'omitnan');
h = bar (mean_S);hold on;
hb = plot ([1], data.data); hold on;
set(hb, 'lineWidth', 3, 'Marker', '.', 'MarkerSize',40);hold on;
%h1=errorbar(mean_S, se_S,'c'); hold on;
set(h,'facecolor', 'none', 'lineWidth', 2);
%set(hb, 'Color','k','linestyle','none', 'lineWidth', 2);
set(gca,'XTick',[1],'XTickLabel',{''}, 'FontSize', 30, 'linew',1, 'xlim', [0 2], 'ylim', [-.25 .05] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 1);
%scatter(1:2, hL+0.135, 'k', 'LineWidth', 6, 'Marker', '*')


[h p ci t] = ttest (data.data);
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth', 1);

export_fig(2, [num2str(layi) '.png'],'-transparent', '-r100');





%% build2compare permutations shuffling trial labels

nPerm = 1000; 
nLays = length(lois);
clear cR max_clust_sum_perm

for permi = 1:nPerm
    for subji = 1:length(pfc)

        x = pfc{subji};
        y = vvs{subji};
        
        % apply shuffling
        idsP = randperm(size(x, 2));
        x = x(:,idsP,:,:);
        
        for layi = 1:nLays
            %x1 = squeeze(mean(mean(x(lois(layi), :, ff1, tt1), 4, 'omitnan'),3, 'omitnan'))'; 
            x1p = squeeze(x(layi, :, :,:)); 
            x1 = mean(x1p(:, clustinfo_lay(7).PixelIdxList{10}),2);
            for layj = 1:nLays
                y1 = squeeze(mean(mean(y(layj, :,ff2, tt2), 4, 'omitnan'),3, 'omitnan'))'; %16:35 max
                %y1p = squeeze(y(lois(layj), :, :,:)); 
                %y1 = mean(y1p(:, clustinfo_lay(7).PixelIdxList{10}),2);

                cR(subji, 1, 1, :) = corr(x1, y1, 'type', 's');

            end
        end  

    end

    [h p ci ts] = ttest(cR);
    h = squeeze(h); t= squeeze(ts.tstat); 
    
    clear allSTs  
    clustinfo = bwconncomp(h);
    if length(clustinfo.PixelIdxList) > 0
        for pxi = 1:length(clustinfo.PixelIdxList)
           allSTs(permi, pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
        end
    else
        allSTs(permi, :) = 0;
    end

    [maxh id] = max(abs(allSTs(permi)));
    max_clust_sum_perm(permi,:) = allSTs(permi,id); 


end

disp ('done')


%% 
t1 = max_clust_sum_real; 
%t1 = [7.6571989];
%t1 = 1.272125091552734e+02;
[counts,centers] = hist(max_clust_sum_perm, 14);

figure()
h = bar(centers, counts); hold on;
h.FaceColor = 'w';
h.EdgeColor = 'k';
h.LineWidth =2;
set(gca, 'FontSize', 20, 'ylim', [0 1000] );
plot ([t1 t1], get(gca, 'ylim'), 'LineWidth', 3,'Color', [.1 .3 .6] );
filename = ['histogram.png'];
export_fig(1, filename,'-transparent', '-r300');
close all;


%% 
t1 = max_clust_sum_real; 
%t1 = [7.6571989];
%t1 = 1.272125091552734e+02;
[counts,centers] = hist(max_clust_sum_perm, 14);

figure()
h = bar(centers, counts); hold on;
h.FaceColor = 'w';
h.EdgeColor = 'k';
h.LineWidth =2;
set(gca, 'FontSize', 20, 'ylim', [0 1000] );
plot ([t1 t1], get(gca, 'ylim'), 'LineWidth', 3,'Color', [.1 .3 .6] );
filename = ['histogram.png'];
export_fig(1, filename,'-transparent', '-r300');
close all;


%% get rank (for stats)

clear p
for layi = 1:nLays
    mcsR = max_clust_sum_real; 
    mcsP = max_clust_sum_perm; 
    allAb = mcsP(abs(mcsP) > abs(mcsR));
    p = 1 - ((nPerm-1) - (length (allAb)))  / nPerm;
end

p




%% check EM2 relationship 
% do first the plot (withohut excluding subjects) to get all_cond1 all_cond2
all_cond1 = all_cond1([7 9 13 18 19 20 21 23 27 28])';
all_cond2 = all_cond2([7 9 13 18 19 20 21 23 27 28])';














%% calculate only in one layer 

L_pfc = 1;
L_vvs = 2;

ff1 = 16:26; 
ff2 = 16:26; 
tt1 = 7:11; 
tt2 = 16:40; 

clear cR
for subji = 1:length(pfc)
   
    x = pfc{subji};
    y = vvs{subji};
    x1 = squeeze(mean(mean(x(L_pfc, :, ff1, tt1), 4, 'omitnan'),3, 'omitnan'))'; 
    y1 = squeeze(mean(mean(y(L_vvs, :, ff2, tt2), 4, 'omitnan'),3, 'omitnan'))'; %16:35 max
    
    figure(); 
    p = polyfit(x1,y1,1); 
    f = polyval(p,x1); 
    plot(x1,y1,'o',x1,f,'-')   
    %scatter(x1,y1);hold on
    cR(subji, :) = corr(x1, y1, 'type', 's');
    
end

[h p ci ts] = ttest(cR, 0, 'Alpha', 0.05);
h = squeeze(h); t= squeeze(ts.tstat); 

max_clust_sum_real = t; 



%% permutations only in one layer 

nPerm = 1000; 
nLays = 7;
clear cR max_clust_sum_perm

for permi = 1:nPerm
    for subji = 1:10
        x = pfc{subji};
        y = vvs{subji};
        
        % apply shuffling
        idsP = randperm(size(x, 2));
        x = x(:,idsP,:,:);
        
        
        x1 = squeeze(mean(mean(x(L_pfc, :, ff1, tt1), 4, 'omitnan'),3, 'omitnan'))'; 
        y1 = squeeze(mean(mean(y(L_vvs, :,ff2, tt2), 4, 'omitnan'),3, 'omitnan'))'; %16:35 max
        cR(subji, 1, 1, :) = corr(x1, y1, 'type', 's');

    end

    [h p ci ts] = ttest(cR);
    h = squeeze(h); t= squeeze(ts.tstat); 
    
    clear allSTs  
    max_clust_sum_perm(permi,:) = t;% 

end

disp ('done')

%% 
t1 = max_clust_sum_real; 
[counts,centers] = hist(max_clust_sum_perm, 14);

figure()
h = bar(centers, counts); hold on;
h.FaceColor = 'w';
h.EdgeColor = 'k';
h.LineWidth =2;
set(gca, 'FontSize', 20, 'ylim', [0 1000] );
plot ([t1 t1], get(gca, 'ylim'), 'LineWidth', 3,'Color', [.1 .3 .6] );
filename = ['histogram.png'];
%export_fig(1, filename,'-transparent', '-r300');
%close all;


%% get rank (for stats)

clear p
for layi = 1:nLays
    mcsR = max_clust_sum_real; 
    mcsP = max_clust_sum_perm; 
    allAb = mcsP(abs(mcsP) > abs(mcsR));
    p = 1 - ((nPerm-1) - (length (allAb)))  / nPerm;
end

p

%% build 2 compare SLIDING
   
pfc = all_r_Times_Trials_PFC([2 3  5  9 10 11 12 14 15 16])';
vvs = all_r_Times_Trials_VVS([7 9 13 18 19 20 21 23 27 28])';


%%compare all2all
%lois         = [1 9 17 25 33 41 49]; 
lois         = [8 16 24 32 40 48 54]; 


ff1 = 16:26; 
tt1 = 7:11; 
ff2 = 16:26; 
tt2 = 6; 

nLays = length(lois);
clear cR
for subji = 1:length(pfc)
   
    x = pfc{subji};
    y = vvs{subji};
    
    for layi = 1:nLays
        x1 = squeeze(mean(mean(x(lois(layi), :, ff1,tt1), 4, 'omitnan'),3, 'omitnan'))'; 
        %x1p = squeeze(x(lois(layi), :, :,:)); 
        %x1 = mean(x1p(:, clustinfo_lay(7).PixelIdxList{10}),2);
        for layj = 1:nLays
            for timei = 1:35
                y1 = squeeze(mean(mean(y(lois(layj), :, ff2,timei: timei+tt2), 4, 'omitnan'),3, 'omitnan'))'; %16:35 max
                %y1p = squeeze(y(lois(layj), :, :,:)); 
                %y1 = mean(y1p(:, clustinfo_lay(7).PixelIdxList{10}),2);

                cR(subji, timei, layi, layj, :) = corr(x1, y1, 'type', 's');
            end
            
        end
    end    
    
end


%%
clear h t
for timei = 1:35
   
cR1 = squeeze(cR(:,timei, :, :));
[h2 p2 ci2 ts] = ttest(cR1);
h(timei,:,:) = squeeze(h2); t(timei, :,:)= squeeze(ts.tstat); 

t1 = squeeze(t(timei,:,:));
h1 = squeeze(h(timei,:,:));

figure()
contourf(myresizem(t1, 20), 40, 'linecolor', 'none'); axis equal, hold on; colorbar
contour(myresizem(h1, 20), 1, 'Color', [0, 0, 0], 'LineWidth', 2);
set(gca, 'xtick', [10 30 50 70 90 110 130], 'xticklabels', {[1:7]})
set(gca, 'ytick', [10 30 50 70 90 110 130], 'yticklabels', {[1:7]})
set(gca, 'FontSize', 22, 'clim', [-4 4]) %
%set(gca, 'xlim',  [0.5 nLays+0.5], 'ylim', [.5  nLays+0.5]) %, 'clim', [0 180]

end





%% plot observed 3D clusters

h3D2 = h; 
%h3D2 = permute(h3D2, [3 2 1]); 

xslice = [0:51];                               % define the cross sections to view
yslice = 0;
zslice = 0;

x = 1:1:7;
y = 1:1:35;
z = 1:1:7;
[X,Y,Z] = meshgrid(x,y,z);

figure(); set(gcf, 'Position', [100 100 400 400]); hold on; 

%slice(X, Y, Z, h3D2, xslice, yslice, zslice)    % display the slices
%contourslice(X, Y, Z, h3D2, xslice, yslice, zslice, 1)    % display the slices

% % % red 3D shape
 fv = isosurface(X, Y, Z, h3D2); 
 p1 = patch(fv);
 p1.FaceColor = 'red';
 p1.EdgeColor = [0.5 0.5 0.5];
 daspect([1 1 1]); 



% % % leave this for the grey plane at 0
 x = 0:.01:54;
 y =  ones(size(x))*6;
 z = ones(size(x))*8;
 col = y;  % This is the color, vary with x in this case.
 p2 = patch([zeros(size(x));z],[y;y],[x;x],[col;col]);
 set(p2,'edgealpha',0.01)

xlim([1 7])
ylim([1 35])
zlim([1 7])
view(34,24)

%set(gca, 'ytick', [6 54], 'yticklabel', {'0' '4'} )
%set(gca, 'ztick', [1 45], 'zticklabel', {'1' '150'}, 'FontSize', 18, 'ylim', [1 25])

exportgraphics(gcf, 'clust3D.png', 'Resolution', 150); 








%% 
%d2p = squeeze(cR(1,:,:,:));
d2p = squeeze(mean(cR, 1));
figure()
subplot(121)
imagesc(d2p); axis equal; colorbar
set(gca, 'xlim', [1 56], 'ylim', [1 56], 'clim', [-.075 .075])



%% here
load RNN_pfc_M&C_noAv_54_Real_T_2extractCluster
load RNN_pfc_M&C_noAv_54_Real_nT_2extractCluster

%lois         = [1 9 17 25 33 41 49]; 
lois         = [8 16 24 32 40 48 56]; %(8th time point)
%lois         = [7 15 23 31 39 47 55]; %(7th time point)
sub2exc =  [];

nSubjs =size(all_r_Times, 1);  
nLays = 7; 
nFreqs = size(all_r_Times, 3); 
clear hL tL clustinfo
for freqi = 1:nFreqs
    allR = squeeze(all_r_Times(:, :, freqi, :));
    allR(sub2exc,:,:) = []; 
    for layi = 1:7
        [hL(layi, freqi, :) p ci ts] = ttest(allR(:,lois(layi),:)); 
        tL(layi, freqi, :) = ts.tstat;
    end
end

%% get cluster in each layer 
f = 1:54;

clear max_clust_sum_obs allSTs clustinfo_lay
for layi = 1:nLays
    
    hLay = squeeze(hL(layi,f,:)); 
    tLay = squeeze(tL(layi,f,:));
    
    clear allSTs  
    clustinfo_lay(layi) = bwconncomp(hLay);
    
end
    

%% extract mean rho values in each trial at this cluster for 1 layer 
apix = [1 1 6 6 4 2 4]; %8th time-point
%apix = [5 8 6 7 4 3 4]; %7th time-point

clear mARTT mARTTCI
for subji = 1:16
   for layi = 1:7
       all_r_tt = squeeze(all_r_Times_Trials{subji}(layi, :,:,:)); 
       nTrials = size(all_r_tt, 1); 
        for triali = 1:nTrials
            all_r_ttt = squeeze(all_r_tt(triali,:,:)); 
            mARTT(triali) = mean(all_r_ttt(clustinfo_lay(layi).PixelIdxList{apix(layi)}), 'omitnan'); 
        end 
        
        ids = idsT{subji};
        ids0 = cellfun(@(x) strsplit(string(x)), ids, 'un', 0);
        ic_i = cellfun(@(x) x(19), ids0, 'un', 0); ic_i = logical(double(string(ic_i)));  % 19 = correct item level ; 20 = correct category level
        mARTTCI(subji,layi, 1) = mean(mARTT(ic_i), 'omitnan');
        mARTTCI(subji,layi, 2) = mean(mARTT(~ic_i), 'omitnan');
   end
end



%% 

for layi = 1:7

    data.data = squeeze(mARTTCI(:,layi,:)); 

    figure(2); set(gcf,'Position', [0 0 500 500]); 
    mean_S = mean(data.data, 1, 'omitnan');
    h = bar (mean_S);hold on;
    hb = plot ([1:2], data.data); hold on;
    set(hb, 'lineWidth', 3, 'Marker', '.', 'MarkerSize',40);hold on;
    %h1=errorbar(mean_S, se_S,'c'); hold on;
    set(h,'facecolor', 'none', 'lineWidth', 2);
    set(hb, 'Color','k','linestyle','none', 'lineWidth', 2);
    set(gca,'XTick',[1:2],'XTickLabel',{'', ''}, 'FontSize', 30, 'linew',1, 'xlim', [0 3], 'ylim', [-.1 .1] );
    plot(get(gca,'xlim'), [0 0],'k','lineWidth', 1);
    %scatter(1:2, hL+0.135, 'k', 'LineWidth', 6, 'Marker', '*')


    [h p ci t] = ttest (data.data);
    disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

    %[p,h,stats] = signrank(data.data(:,1), data.data(:,2));
    %disp (['W = ' num2str(stats.signedrank) '  ' ' p = ' num2str(p)]);


    set(gca, 'LineWidth', 1);

    export_fig(2, [num2str(layi) '.png'],'-transparent', '-r100');
    close all;   
    
end










%%

figure()
layT = tiledlayout(8, 7);
layT.TileSpacing = 'compact';
layT.Padding = 'compact';
for layi = 1:nLays
    
    t = squeeze(tL(layi,:,:));
    h = squeeze(hL(layi,:,:));
    
    if size(all_r_Times, 4) < 30 & size(all_r_Times, 4) > 5
        % % % % encoding 
        if nLays > 9
            subplot(7, 8, layi)
        else
            subplot(3, 3, layi)
        end
        %times = -.5:.01:1.499;
        times = 0:.01:1.499;
        %freqs = [.1:.1:29.9 30:.5:150]; 
        freqs = 1:540;
        contourf(times, freqs, myresizem(t, 10), 100, 'linecolor', 'none'); hold on; %colorbar
        contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 2);
        plot([-0.055 -0.055], get(gca, 'ylim'), 'k:', 'LineWidth', 2);
        plot(get(gca, 'xlim'),[300 300], 'k:', 'LineWidth', 2);
        set(gca, 'xtick', [-0.055 .445 .945 1.45], 'xticklabel', [0 .5 1 1.5], 'clim', [-5 5]); colormap jet
        set(gca,'XTick',[], 'YTick', [], 'xticklabel',[], 'FontSize', 12)

    elseif size(all_r_Times, 4) == 5
        % % % % encoding 
        if nLays > 9
            subplot(7, 8, layi)
        else
            subplot(3, 3, layi)
        end
        times = 0:.01:0.499;
        %freqs = [.1:.1:29.9 30:.5:150]; 
        freqs = 1:540;
        contourf(times, freqs, myresizem(t, 10), 100, 'linecolor', 'none'); hold on; %colorbar
        contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 2);
        plot([-0.055 -0.055], get(gca, 'ylim'), 'k:', 'LineWidth', 2);
        plot(get(gca, 'xlim'),[300 300], 'k:', 'LineWidth', 2);
        set(gca, 'xtick', [-0.055 .445 .945 1.45], 'xticklabel', [0 .5 1 1.5], 'clim', [-5 5]); colormap jet
        set(gca,'XTick',[], 'YTick', [], 'xticklabel',[], 'FontSize', 12)
    else
        if nLays > 9
            subplot(7, 8, layi)
        else
            subplot(3, 3, layi)
        end
        %times = 0:.01:3.59;
        %times = 0:.01:4.49;
        %freqs = [.1:.1:29.9 30:.5:150]; 
        freqs = 1:540;
        contourf(times, freqs, myresizem(t, 10), 100, 'linecolor', 'none'); hold on; 
        contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 1); %colorbar
        plot([-0.055 -0.055], get(gca, 'ylim'), 'k:', 'LineWidth', 2);
        plot(get(gca, 'xlim'),[300 300], 'k:', 'LineWidth', 2);
        set(gca, 'xtick', [-0.055 .945 1.945 2.945 3.945 ], 'xticklabel', [0 1 2 3 4], 'xlim', [0 3.5], 'clim', [-3 3])
        set(gca,'XTick',[], 'YTick', [])
        set(gca,'xticklabel',[], 'FontSize', 12)
        colormap jet
    end

end
        

f2p = [num2str(layi, '%02.f') '.png'];
exportgraphics(gcf, f2p, 'Resolution', 300)






%% DNN ANALYSIS BAND PERMUTATIONS
clearvars -except act_CH act_FR 
nPerm               = 1000; 
f2sav               = 'RNN_vvs_Maint_noAv_B_PERM'
f                   = 1:54; %in case B

load_parameters_WMFRA;
loadNet_WM;
all_r_Times_perm    = zeros(nPerm, nSubj, nLays, nFreqs, nTimes);

tic
for subji=1:nSubj
    disp (['Subj: ' num2str(subji)])
    load(sublist{subji});
    cfg_contrasts.oneListPow = cfg_contrasts.oneListPow(:, :, :, resTime); 
    cfg_contrasts = normalize_WM(cfg_contrasts, 1, 'blo', [3 7]);
    if strcmp(avMeth,'pow')
       cfg_contrasts = average_repetitions(cfg_contrasts);
    end
    
    
    f2t = strsplit(f2sav, '_'); 
    if strcmp(f2t{5}, '54') | strcmp(f2t{5}, '150') 
        for freqi = 1:length(freqs2test)
            f  = freqs2test(freqi);
            rdm_prev(freqi) = create_rdms(cfg_contrasts, f, it, 5, 1);
        end
    else
        rdm_prev = create_rdms(cfg_contrasts, f, it, 5, 1);
        if subji == 1 % change name only once
            f2sav = [f2sav(1:end-4) num2str(f(1)) '-' num2str(f(end)) 'Hz_PERM'];
        end
    end
    
    ids = rdm_prev(1).ids; % take first freq only to get ids
    ids = cellfun(@(x) strsplit(string(x)), ids, 'un', 0);
    ids0 = cellfun(@(x) x(3), ids, 'un', 0); ids0 = double(string(ids0)); 
    ids1 = char(string((ids0))); ids2 = ids1(:,[1 3]);ids3 = double(string(ids2));
    idx = find(~mod(ids3, 10)); ids4 = ids3-10; ids4(idx) = ids3(idx);


    
    allS = cat(4, rdm_prev.rdm);
    allS2 = allS;allS2(allS2==1) = 1000; allS2(allS2==0) = 2000;
    allS3 = arrayfun(@(j)  (arrayfun(@(i)tril(squeeze(allS2(:,:,i,j)), -1), 1:size(allS2,3), 'un', 0)), 1:size(allS2,4), 'un', 0);
    allS4 = vertcat(allS3{:});
    allS7 = cellfun(@(x) x(x~=0), allS4, 'un', 0);     
    allS8 = cat(2, allS7 {:});
    allS8 (allS8 ==1000) = 1;allS8 (allS8 ==2000) = 0;

    for permi = 1:nPerm
        parfor layi = 1:nLays
            if subji < subj_ch_fr
                M =  squeeze(act_FR(layi,ids4,ids4)); 
            else
                M =  squeeze(act_CH(layi,ids4,ids4)); 
            end
            idp = randperm(length(ids0));
            Mp = M(idp, idp);
            Mp(Mp ==1) = 1000;Mp(Mp ==0) = 2000;
            Mp = tril(Mp, -1);Mp(Mp==0) =[];
            Mp(Mp==1000) = 1;Mp(Mp==2000) = 0;
            allTEst = corr(allS8, Mp', 'type', 's');
            all_r_Times_perm(permi, subji, layi,:, : ) = reshape(allTEst, nFreqs, nTimes );  
        end
    end
 
end 

for permi = 1:nPerm
    for layi = 1:nLays
        d2u = squeeze(all_r_Times_perm(permi,:, layi, :, :));
        [hs p ci ts] =   ttest(d2u);
        hPerm(permi, layi,:,:) = squeeze(hs); 
        tPerm(permi, layi,:,:) = squeeze(ts.tstat); 
    end
end     

cd (currentF) 
toc 
save(f2sav, 'hPerm', 'tPerm', 'all_r_Times_perm');






%%  Exclude subjects from all_r_Times_perm file 
subj2exc = [1]; 
for permi = 1:nPerm
    for layi = 1:nLays
        d2u = squeeze(all_r_Times_perm(permi,:, layi, :, :));
        d2u(subj2exc, :, :) = []; 
        [hs p ci ts] =   ttest(d2u);
        hPerm(permi, layi,:,:) = squeeze(hs); 
        tPerm(permi, layi,:,:) = squeeze(ts.tstat); 
    end
end  

%% check 1 permutation (onlyt t and h)

perm = 37; 

subj2exc = []; 
cfg.yLim = [-.015 .025]; 
cfg.sP = 1;
cfg.Fz = 8; 
cfg.Lw = 4; %zero = not drawn 


figure()
for layi = 1:nLays
    subplot(3, 3, layi); 
    r_Time = squeeze(all_r_Times_perm(perm, :, layi, :,:)); 
    r_Time(subj2exc, : ) = []; 

    [max_clust_sum_obs_check(layi)] = plot_with_sigV(r_Time, cfg); hold on;    

end

 %exportgraphics(gcf, 'figure.png', 'Resolution', 300)




%% count clusters
nPerm = size(tPerm, 1); 

clear max_clust_sum 
for permi = 1:nPerm
clear allSTs   
    for layi = 1:nLays

        d2pT = squeeze(tPerm(permi, layi,:,:));
        d2pH = squeeze(hPerm(permi, layi,:,:));

        clustinfo = bwconncomp(d2pH);
        if length(clustinfo.PixelIdxList) > 0
            for pxi = 1:length(clustinfo.PixelIdxList)
               allSTs(layi, permi, pxi) = sum(d2pT(clustinfo.PixelIdxList{pxi}));% it is maxi -1 cause we are adding a zero below
            end
        else
            allSTs(layi, permi, :) = 0;
        end

        [maxh id] = max(abs(allSTs(layi, permi, :)));
        max_clust_sum_perm(permi,layi, :) = allSTs(layi, permi,id); 
    end
    
end

%%

for layi = 1:nLays
    
    t1 = max_clust_sum_obs(layi); 
    [counts,centers] = hist(max_clust_sum_perm(:,layi), 14);
    
    figure()
    h = bar(centers, counts); hold on;
    h.FaceColor = 'w';
    h.EdgeColor = 'k';
    h.LineWidth =2;
    set(gca, 'FontSize', 20, 'ylim', [0 300] );
    plot ([t1 t1], get(gca, 'ylim'), 'LineWidth', 3,'Color', [.1 .3 .6] );
    filename = ['histogram.png'];
    %export_fig(1, filename,'-transparent', '-r300');
    %close all;
end


%% get rank (for stats)

clear p
for layi = 1:nLays
    mcsR = max_clust_sum_obs(layi); 
    mcsP = max_clust_sum_perm(:,layi); 
    allAb = mcsP(abs(mcsP) > abs(mcsR));
    p(layi, :) = 1 - ((nPerm-1) - (length (allAb)))  / nPerm;
end

p




%% PFC cue performance effect behavioral

clear
load RNN_pfc_M_Av_R_nT_1-49_C_naFV



sub2exc =  [];

nSubjs =size(all_r_Times, 1);  
nFreqs = size(all_r_Times, 3); 
clear hL tL clustinfo
allR = all_r_Times;
allR(sub2exc,:,:) = []; 
for layi = 1:56
    [hL(layi, :, :) p ci ts] = ttest(allR(:,layi,:,:)); 
    tL(layi, :, :) = ts.tstat;
end


%%get cluster in each layer 
clear max_clust_sum_obs allSTs clustinfo_lay
for layi = 1:56
    hLay = squeeze(hL(layi,:,:)); 
    tLay = squeeze(tL(layi,:,:));
    clustinfo_lay(layi) = bwconncomp(hLay);
end
    

%%extract mean rho values in each trial at this cluster for 1 layer 
clear mARTT 
for subji = 1:16
   for layi = 7:7
       all_r_tt = squeeze(all_r_Times_Trials{subji}(layi, :,:,:)); 
       nTrials = size(all_r_tt, 1); 
        for triali = 1:nTrials
            all_r_ttt = squeeze(all_r_tt(triali,:,:)); 
            mARTT{subji,:}(triali,:) = mean(all_r_ttt(clustinfo_lay(56).PixelIdxList{10}), 'omitnan'); 
        end   
   end
end

%% CUEPERF EFFECT in one layer

clear cR xC xI perf_cc
for subji = 1:16
   
    clear ic_i
    x = mARTT{subji};
    
    
    ids = idsT{subji};
    ids0 = cellfun(@(x) strsplit(string(x)), ids, 'un', 0);
    % 19 = correct item level ; 20 = correct category level
    ic_i = cellfun(@(x) x(20), ids0, 'un', 0); ic_i = logical(double(string(ic_i))); 
    
    xC{subji,:}= x(ic_i); 
    xI{subji,:}= x(~ic_i); 
    
end


for subji = 1:16   
    perf_cc(subji, 1) = mean(xC{subji});
    perf_cc(subji, 2) = mean(xI{subji});
end



%% 2 bar

data.data = perf_cc; 

sub2exc = [1 7 13 16];

data.data(sub2exc,:) =  [];

figure(2); set(gcf,'Position', [0 0 500 500]); 
mean_S = mean(data.data, 1, 'omitnan');
h = bar (mean_S);hold on;
hb = plot ([1:2], data.data); hold on;
set(hb, 'lineWidth', 3, 'Marker', '.', 'MarkerSize',40);hold on;
%h1=errorbar(mean_S, se_S,'c'); hold on;
set(h,'facecolor', 'none', 'lineWidth', 2);
%set(hb, 'Color','k','linestyle','none', 'lineWidth', 2);
set(gca,'XTick',[1:2],'XTickLabel',{'', ''}, 'FontSize', 30, 'linew',1, 'xlim', [0 3], 'ylim', [-.15 .15] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 1);
%scatter(1:2, hL+0.135, 'k', 'LineWidth', 6, 'Marker', '*')


[h p ci t] = ttest (data.data(:,1), data.data(:,2));
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth', 1);

export_fig(2, [num2str(layi) '.png'],'-transparent', '-r100');




%% 

for layi = 1:7

    data.data = squeeze(perf_cc(2:end,layi,:)); 

    figure(2); set(gcf,'Position', [0 0 500 500]); 
    mean_S = mean(data.data, 1, 'omitnan');
    h = bar (mean_S);hold on;
    hb = plot ([1:2], data.data); hold on;
    set(hb, 'lineWidth', 3, 'Marker', '.', 'MarkerSize',40);hold on;
    %h1=errorbar(mean_S, se_S,'c'); hold on;
    set(h,'facecolor', 'none', 'lineWidth', 2);
    set(hb, 'Color','k','linestyle','none', 'lineWidth', 2);
    set(gca,'XTick',[1:2],'XTickLabel',{'', ''}, 'FontSize', 30, 'linew',1, 'xlim', [0 3], 'ylim', [-.02 .02] );
    plot(get(gca,'xlim'), [0 0],'k','lineWidth', 1);
    %scatter(1:2, hL+0.135, 'k', 'LineWidth', 6, 'Marker', '*')


    [h p ci t] = ttest (data.data);
    disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

    set(gca, 'LineWidth', 1);

    export_fig(2, [num2str(layi) '.png'],'-transparent', '-r100');
    close all;   
    
end








%% CREATE RDM _ PREV with all frequencies

clearvars 
freqs2test = [1:54]';

win_width           = 5; 
mf                  = 1; 
meanInTime          = 1; 
meanInFreq          = 0; 
takeElec            = 0; 
takeFreq            = 0;
TG                  = 1; %temporal generalization
contr2save          = {'ALL_EE'}; %{};use 'ALL_EE' for rdms
bline               = [3 7];
acrossTrials        = 1;
zScType             = 'blo'; %'sess' %'blo' % 'none' % 'none' all trials irrespective of sesson or block
avMeth              = 'pow'; %'pow' 'corr' 'none'



tic
sublist = dir('*contr.mat');
sublist = {sublist.name};
disp (['measurements -> ' num2str(length(sublist))]);


% % % % correlate at each frequency
for freqi = 1:length(freqs2test)
      disp (['fq: ' num2str(freqi)])
      for subji=1:length(sublist)
            %disp(['File > ' num2str(subji)]);
            load(sublist{subji});
        
            
            %cfg_contrasts = averageBands_WM(cfg_contrasts);
            
            %cut time
            %cfg_contrasts.oneListPow = cfg_contrasts.oneListPow(:, :, :, 51:300); 
            
            cfg_contrasts = normalize_WM(cfg_contrasts, acrossTrials, zScType, bline);

            if strcmp(avMeth,'pow')
                cfg_contrasts = average_repetitions(cfg_contrasts);
            end
            
            cfg_contrasts.contr2save = contr2save;

            f  = freqs2test(freqi);
            rdm_prev = create_rdms(cfg_contrasts, f, {'7'}, win_width, mf, meanInTime);
            
            
            if strcmp(avMeth,'corr')
                rdm_prev = average_repetitions_corr(rdm_prev);
            end

            all_rdm_prev(subji, freqi, :, :) = rdm_prev;
            
            

      end
end

filename = 'rdms_1-54_Maint_Av';
save(filename, 'all_rdm_prev', '-v7.3');

toc


%% FREQUENCY RESOLVED DNN ANALYSIS (SLOW WAY AS SANITY CHECK) 

clearvars -except act_CH act_FR 
f2sav       = 'RNN_pfc_M_noAv_R_nT_1-49_C_aFV_100ms'; %'Alex_vvs_CueAll_noAv_54_Real_nT_1-49'
load_parameters_WMFRA;
loadNet_WM;

all_r_Times    = zeros(nSubj, nLays, nFreqs, nTimes);

tic

for subji= 1:nSubj %start with subject 2
    disp (['Subj: ' num2str(subji)])
    load(sublist{subji});
    cfg_contrasts.oneListPow = cfg_contrasts.oneListPow(:, :, :, resTime); 
    cfg_contrasts = normalize_WM(cfg_contrasts, 1, 'blo', [3 7]);
    if strcmp(avMeth,'pow')
       cfg_contrasts = average_repetitions(cfg_contrasts);
    end
    
    for freqi = 1:length(freqs2test)
        f  = freqs2test(freqi);
        
        if strcmp(f2t{9}, 'aFV')
            avTW = 1;
        else
            avTW = 0;
        end
        rdm_prev = create_rdms(cfg_contrasts, f, it, 1, 1, avTW);%
        
       
        [allS ids act_FR2 act_CH2 cM] = getIDsAndnRDMs(rdm_prev, act_CH, act_FR);
        
        for timei = 1:45
            rdmp = squeeze(rdm_prev.rdm);
            rdm = squeeze(rdmp(:,:,timei));
            rdm(rdm ==1) = 1000;rdm(rdm ==0) = 2000;
            rdm = tril(rdm, -1);rdm(rdm==0) =[];
            rdm(rdm==1000) = 1;rdm(rdm==2000) = 0;
        
            for layi = 56:56
                if subji < subj_ch_fr
                    M =  squeeze(act_FR2(layi,:,:)); 
                else
                    M =  squeeze(act_CH2(layi,:,:)); 
                end
                M(M ==1) = 1000;M(M ==0) = 2000;
                M = tril(M, -1);M(M==0) =[];
                M(M==1000) = 1;M(M==2000) = 0;

                cMh = cM; 
                cMh(cMh ==1) = 1000;cMh(cMh ==0) = 2000;
                cMh = tril(cMh, -1);cMh(cMh==0) =[];
                cMh(cMh==1000) = 1;cMh(cMh==2000) = 0;

                
                if strcmp(f2t{8}, 'PC')
                    x = [rdm' M']; 
                    z = cMh'; 
                    allTEst_prev = partialcorr(x, z, 'type', 's');
                    allTEst = allTEst_prev(2);
                else
                    allTEst = corr(rdm', M', 'type', 's');
                end
                
                

                all_r_Times(subji, layi,freqi,timei) = allTEst;  
            end
        end
    end
    
    
    
end 

cd (currentF) 
toc 
save(f2sav, 'all_r_Times');

%% plot traces sanity check
d2pp = mean(all_r_Times(:, 56, 13:29, :),3, 'omitnan');
d2p = squeeze(d2pp);

m = mean(d2p, 'omitnan');
se = std(d2p, 'omitnan') / sqrt(size(d2p, 1));
figure()
errorbar(m, se); % group mean


d2p = squeeze(mean(all_r_Times(:,56,:,:), 'omitnan'));
figure()
imagesc(d2p);

%% check
load hospital;
hospital.SexID = grp2idx(hospital.Sex);
x = [hospital.Age hospital.BloodPressure];
z = [hospital.SexID hospital.Smoker];
[rho,pval] = partialcorr(x,z)

%% FREQUENCY RESOLVED DNN ANALYSIS

clearvars -except act_CH act_FR  
f2sav       = 'RNN_pfc_M_Av_R_nT_1-550_C_aFV_10ms'; %'Alex_vvs_CueAll_noAv_54_Real_nT_1-49'
load_parameters_WMFRA;
loadNet_WM;

all_r_Times    = zeros(nSubj, nLays, nFreqs, nTimes);

tic
f2u = strsplit(f2sav, '_');
for subji=1:nSubj
    disp (['Subj: ' num2str(subji)])
    load(sublist{subji});
    cfg_contrasts.oneListPow = cfg_contrasts.oneListPow(:, :, :, resTime); 
    cfg_contrasts = normalize_WM(cfg_contrasts, 1, 'blo', [3 7]);
    if strcmp(avMeth,'pow')
       cfg_contrasts = average_repetitions(cfg_contrasts);
    end
    
    for freqi = 1:length(freqs2test)
        f  = freqs2test(freqi);
        if strcmp(f2u{9}, 'aFV')
            avTW = 1;
        else
            avTW = 0;
        end
        rdm_prev(freqi) = create_rdms(cfg_contrasts, f, it, 1, 1, avTW);%
    end
    
    [allS ids act_FR2 act_CH2 cM] = getIDsAndnRDMs(rdm_prev, act_CH, act_FR);
    
    parfor layi = 1:nLays
        if subji < subj_ch_fr
            M =  squeeze(act_FR2(layi,:,:)); 
        else
            M =  squeeze(act_CH2(layi,:,:)); 
        end
        
        M(M ==1) = 1000;M(M ==0) = 2000;
        M = tril(M, -1);M(M==0) =[];
        M(M==1000) = 1;M(M==2000) = 0;
        
        cMh = cM; 
        cMh(cMh ==1) = 1000;cMh(cMh ==0) = 2000;
        cMh = tril(cMh, -1);cMh(cMh==0) =[];
        cMh(cMh==1000) = 1;cMh(cMh==2000) = 0;

        if strcmp(f2u{8}, 'PC')
            allTEst = partialcorri(allS, M', cMh', 'type', 's');    
        else
           allTEst = corr(allS, M', 'type', 's'); 
        end
        
        all_r_Times(subji, layi,:, : ) = reshape(allTEst, nFreqs, nTimes );  
    end
    
    
    
end 

cd (currentF) 

save(f2sav, 'all_r_Times');
toc 


%% FREQUENCY RESOLVED DNN ANALYSIS PARTIAL CORRELATION ALL LAYERS

clearvars -except act_CH act_FR 
f2sav       = 'RNN_pfc_M&C_noAv_54_Real_nT_1-49_PCAL_aFV'; %'Alex_vvs_CueAll_noAv_54_Real_nT_1-49'
load_parameters_WMFRA;
loadNet_WM;

all_r_Times    = zeros(nSubj, nLays, nFreqs, nTimes);

tic
f2u = strsplit(f2sav, '_');
for subji=1:nSubj
    disp (['Subj: ' num2str(subji)])
    load(sublist{subji});
    cfg_contrasts.oneListPow = cfg_contrasts.oneListPow(:, :, :, resTime); 
    cfg_contrasts = normalize_WM(cfg_contrasts, 1, 'blo', [3 7]);
    if strcmp(avMeth,'pow')
       cfg_contrasts = average_repetitions(cfg_contrasts);
    end
    
    for freqi = 1:length(freqs2test)
        f  = freqs2test(freqi);
        if strcmp(f2u{10}, 'aFV')
            avTW = 1;
        else
            avTW = 0;
        end
        rdm_prev(freqi) = create_rdms(cfg_contrasts, f, it, 5, 1, avTW);%
    end
    
    [allS ids act_FR2 act_CH2 cM] = getIDsAndnRDMs(rdm_prev, act_CH, act_FR);
    lois = (1:8:56) + 7;
    
    for layi = 56:56
        if subji < subj_ch_fr
            M1 = squeeze(act_FR2(lois(1),:,:)); 
            M2 = squeeze(act_FR2(lois(2),:,:)); 
            M3 = squeeze(act_FR2(lois(3),:,:)); 
            M4 = squeeze(act_FR2(lois(4),:,:)); 
            M5 = squeeze(act_FR2(lois(5),:,:)); 
            M6 = squeeze(act_FR2(lois(6),:,:)); 
            M7 = squeeze(act_FR2(lois(7),:,:)); 
        else
            M1 = squeeze(act_CH2(lois(1),:,:)); 
            M2 = squeeze(act_CH2(lois(2),:,:)); 
            M3 = squeeze(act_CH2(lois(3),:,:)); 
            M4 = squeeze(act_CH2(lois(4),:,:)); 
            M5 = squeeze(act_CH2(lois(5),:,:)); 
            M6 = squeeze(act_CH2(lois(6),:,:)); 
            M7 = squeeze(act_CH2(lois(7),:,:)); 
        end
        
        M1(M1 ==1) = 1000;M1(M1 ==0) = 2000;
        M1 = tril(M1, -1);M1(M1==0) =[];
        M1(M1==1000) = 1;M1(M1==2000) = 0;
        
        M2(M2 ==1) = 1000;M2(M2 ==0) = 2000;
        M2 = tril(M2, -1);M2(M2==0) =[];
        M2(M2==1000) = 1;M2(M2==2000) = 0;
        
        M3(M3 ==1) = 1000;M3(M3 ==0) = 2000;
        M3 = tril(M3, -1);M3(M3==0) =[];
        M3(M3==1000) = 1;M3(M3==2000) = 0;
        
        M4(M4 ==1) = 1000;M4(M4 ==0) = 2000;
        M4 = tril(M4, -1);M4(M4==0) =[];
        M4(M4==1000) = 1;M4(M4==2000) = 0;
        
        M5(M5 ==1) = 1000;M5(M5 ==0) = 2000;
        M5 = tril(M5, -1);M5(M5==0) =[];
        M5(M5==1000) = 1;M5(M5==2000) = 0;
        
        M6(M6 ==1) = 1000;M6(M6 ==0) = 2000;
        M6 = tril(M6, -1);M6(M6==0) =[];
        M6(M6==1000) = 1;M6(M6==2000) = 0;
        
        M7(M7 ==1) = 1000;M7(M7 ==0) = 2000;
        M7 = tril(M7, -1);M7(M7==0) =[];
        M7(M7==1000) = 1;M7(M7==2000) = 0;
        
        cMh = [M1; M2; M3; M4; M5; M6]; 
        
        

        %allTEst = corr(allS, M', 'type', 's');
        allTEst = partialcorri(allS, M7', cMh', 'type', 's');
        
        all_r_Times(subji, layi,:, : ) = reshape(allTEst, nFreqs, nTimes );  
    end
    
    
    
end 

cd (currentF) 
toc 
save(f2sav, 'all_r_Times');




%% PLOT OBS DATA (all layers)

%sub2exc =  [1 4 6:8 13 ]; %[1 3 4 7 10 13 16];
sub2exc =[];

nSubjs =size(all_r_Times, 1);  
nLays = size(all_r_Times, 2); 
nFreqs = size(all_r_Times, 3); 
clear hL tL clustinfo
for freqi = 1:nFreqs
    allR = squeeze(all_r_Times(:, :, freqi, :));
    allR(sub2exc,:,:) = []; 
    for layi = 1:nLays
        [hL(layi, freqi, :) p ci ts] = ttest(allR(:,layi,:)); 
        tL(layi, freqi, :) = ts.tstat;
    end
end

f = 1:52;

clear max_clust_sum_obs allSTs
for layi = 1:nLays
    
    hLay = squeeze(hL(layi,f,:)); 
    tLay = squeeze(tL(layi,f,:));
    
    clear allSTs  
    clustinfo = bwconncomp(hLay);
    for pxi = 1:length(clustinfo.PixelIdxList)
        % check whether it is a combined + and - cluster
        V = tLay(clustinfo.PixelIdxList{pxi});
        Vids = clustinfo.PixelIdxList{pxi}; 
        if ~any(diff(sign(V(V~=0)))) %all time bins have tvalues of same sie
            allSTs(pxi,:) = sum(tLay(clustinfo.PixelIdxList{pxi}));
        else %remove the 
            big0 = V>0; small0 = V<0; 
            VidsS = Vids(small0); 
            ids2k = Vids(big0); 
            if sum(big0) > sum(small0)
                V(small0) = NaN; 
                hLay(VidsS) = 0; 
                clustinfo.PixelIdxList{pxi} = ids2k; 
            else
                %V(big0) = NaN; 
                %hLay(big0,:) = NaN; 
            end
            
            allSTs(pxi,:) = sum(V, 'omitnan'); 
        end
    end

    [maxh id] = max(abs(allSTs));
    max_clust_sum_obs(layi,:) = allSTs(id); 

    % % % % rem non-sig-clusters
    for ci = 1:length(clustinfo.PixelIdxList)
        modifiedCluster = abs(sum(tLay(clustinfo.PixelIdxList{ci})))+0.0001; %add this small number because removing negative values creates slighlly different valeus
       if modifiedCluster < max(abs(allSTs))   + 1    %add +1 at the very end to delete all clusters
          %disp(['layi > ' num2str(layi) ' > ' num2str(abs(sum(tLay(clustinfo.PixelIdxList{ci})))) '_' num2str(max(abs(allSTs)))]) 
          %hLay (clustinfo.PixelIdxList{ci}) =  0; 
       end
    end
    
    hL(layi, f, :) = hLay; 
    tL(layi, f, :) = tLay; 
end

% get 3D cluster
h3D = squeeze(hL(:,f,:)); 
t3D = squeeze(tL(:,f,:));
clear allSTs   
clustinfo = bwconncomp(h3D, 6);
if length(clustinfo.PixelIdxList) > 0
    for pxi = 1:length(clustinfo.PixelIdxList)
       allSTs(pxi) = sum(t3D(clustinfo.PixelIdxList{pxi}));% it is maxi -1 cause we are adding a zero below
    end
else
    allSTs = 0;
end
[maxh id] = max(abs(allSTs));
max_clust_sum_3D = allSTs(id); 
for ci = 1:length(clustinfo.PixelIdxList)
   modifiedCluster = abs(sum(t3D(clustinfo.PixelIdxList{ci})))+0.0001; %add this small number because removing negative values creates slighlly different valeus
   if modifiedCluster < max(abs(allSTs))  % + 1    %add +1 at the very end to delete all clusters
      %disp(['layi > ' num2str(layi) ' > ' num2str(abs(sum(tLay(clustinfo.PixelIdxList{ci})))) '_' num2str(max(abs(allSTs)))]) 
      h3D (clustinfo.PixelIdxList{ci}) =  0; 
   end
end



figure()
layT = tiledlayout(8, 7);
layT.TileSpacing = 'compact';
layT.Padding = 'compact';
for layi = 56:56 %1:nLays
    
    t = squeeze(tL(layi,:,:));
    h = squeeze(hL(layi,:,:));
    
    if size(all_r_Times, 4) < 30 & size(all_r_Times, 4) > 5
        % % % % encoding 
        if nLays > 9
            subplot(7, 8, layi)
        else
            subplot(3, 3, layi)
        end
        %times = -.5:.01:1.499;
        times = 0:.01:1.499;
        %freqs = [.1:.1:29.9 30:.5:150]; 
        freqs = 1:520;
        contourf(times, freqs, myresizem(t, 10), 100, 'linecolor', 'none'); hold on; %colorbar
        contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 2);
        plot([-0.055 -0.055], get(gca, 'ylim'), 'k:', 'LineWidth', 2);
        plot(get(gca, 'xlim'),[300 300], 'k:', 'LineWidth', 2);
        set(gca, 'xtick', [-0.055 .445 .945 1.45], 'xticklabel', [0 .5 1 1.5], 'clim', [-5 5]); colormap jet
        set(gca,'XTick',[], 'YTick', [], 'xticklabel',[], 'FontSize', 12)

    elseif size(all_r_Times, 4) == 5
        % % % % encoding 
        if nLays > 9
            subplot(7, 8, layi)
        else
            subplot(3, 3, layi)
        end
        times = 0:.01:0.499;
        %freqs = [.1:.1:29.9 30:.5:150]; 
        freqs = 1:520;
        contourf(times, freqs, myresizem(t, 10), 100, 'linecolor', 'none'); hold on; %colorbar
        contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 2);
        plot([-0.055 -0.055], get(gca, 'ylim'), 'k:', 'LineWidth', 2);
        plot(get(gca, 'xlim'),[300 300], 'k:', 'LineWidth', 2);
        set(gca, 'xtick', [-0.055 .445 .945 1.45], 'xticklabel', [0 .5 1 1.5], 'clim', [-5 5]); colormap jet
        set(gca,'XTick',[], 'YTick', [], 'xticklabel',[], 'FontSize', 12)
    elseif size(all_r_Times, 4) == 546
        % % % % encoding 
        if nLays > 9
            subplot(7, 8, layi)
        else
            subplot(3, 3, layi)
        end
        times = 0:.01:5.009;
        %freqs = [.1:.1:29.9 30:.5:150]; 
        freqs = 1:52;
        contourf(times, freqs, t, 100, 'linecolor', 'none'); hold on; %colorbar
        contour(times, freqs, h, 1, 'Color', [0, 0, 0], 'LineWidth', 2);
        %plot([-0.055 -0.055], get(gca, 'ylim'), 'k:', 'LineWidth', 2);
        %plot(get(gca, 'xlim'),[300 300], 'k:', 'LineWidth', 2);
        %set(gca, 'xtick', [-0.055 .445 .945 1.45], 'xticklabel', [0 .5 1 1.5], 'clim', [-5 5]); colormap jet
        %set(gca,'XTick',[], 'YTick', [], 'xticklabel',[], 'FontSize', 12)
    else
        if nLays > 9
            subplot(7, 8, layi)
        else
            subplot(3, 3, layi)
        end
        %times = 0:.01:3.49;
       % times = 0:.01:4.49;
        times = -.5:.01:3.99;
        %freqs = [.1:.1:29.9 30:.5:150]; 
        freqs = 1:520;
        contourf(times, freqs, myresizem(t, 10), 100, 'linecolor', 'none'); hold on; 
        contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 1); %colorbar
        plot([-0.055 -0.055], get(gca, 'ylim'), 'k:', 'LineWidth', 2);
        plot(get(gca, 'xlim'),[300 300], 'k:', 'LineWidth', 2);
        set(gca, 'xtick', [-0.055 .945 1.945 2.945 3.945 ], 'xticklabel', [0 1 2 3 4], 'clim', [-3 3]); % 'xlim', [0 3.5],
        set(gca,'XTick',[], 'YTick', [])
        set(gca,'xticklabel',[], 'FontSize', 12)
        colormap jet
    end

end
        

%f2p = [num2str(layi, '%02.f') '.png'];
%exportgraphics(gcf, f2p, 'Resolution', 300)



%% beta effect time window extracted from EM2 contrast

load RNN_vvs_Cue_noAv_54_Real_nT_1-49

foi = 13:27;
toi = 16:27; 

clear e2 h t
for subji = 1:28
    for layi = 1:56
        e2(subji, layi,:,:) = squeeze(mean(all_r_Times(subji, layi,foi,toi), 'all'));
        
    end
end


mtDiff = mean(e2,'omitnan');
stdDiff = std(e2,'omitnan');
nsy = size(e2, 1) -1; 
seDiff = stdDiff / sqrt(nsy);

lois = (1:8:56) + 7;
figure()
%plot(mtDiff(lois), 'LineWidth', 2)
%shadedErrorBar(1:7, mtDiff(lois), seDiff(lois), 'r', 1); hold on; 
%plot(mtDiff, 'LineWidth', 2)
shadedErrorBar(1:56, mtDiff, seDiff, 'r', 1); hold on; 
set(gca, 'FontSize', 20)%, 'ylim', [-.3 .8]
exportgraphics(gcf, 'layers.png', 'Resolution', 300)


[h p ci ts] = ttest(e2(:,lois));
%[h p ci ts] = ttest(tDiff);
p



exportgraphics(gcf, 'layers.png', 'Resolution', 300)










%% mean beta effect all subjects
cd D:\_WM\analysis\predictor_matrix_analyses\pfc\beta_effect
load RNN_pfc_Maint_noAv_54_
load Cat_pfc_Maint_noAv_54_

foi = 13:29;
toi = 6:11; 

for subji = 1:16
    for layi = 1:56

        e1 = squeeze(all_r_Times_catM(subji, :,:));
        e2 = squeeze(all_r_Times_RNN(subji, layi,:,:));
        tDiff(subji, layi) = mean(e2(foi, toi), 'all') - mean(e1(foi, toi), 'all') ;

    end
end


mtDiff = mean(tDiff,'omitnan');
stdDiff = std(tDiff,'omitnan');
seDiff = stdDiff / sqrt(15);

lois = (1:8:56) + 7;
figure()
%plot(mtDiff(lois), 'LineWidth', 2)
%shadedErrorBar(1:7, mtDiff(lois), seDiff(lois), 'r', 1); hold on; 
%plot(mtDiff, 'LineWidth', 2)
shadedErrorBar(1:56, mtDiff, seDiff, 'r', 1); hold on; 
set(gca, 'FontSize', 20)%, 'ylim', [-.3 .8]
exportgraphics(gcf, 'layers.png', 'Resolution', 300)


[h p ci ts] = ttest(tDiff(:,lois));
%[h p ci ts] = ttest(tDiff);
p



exportgraphics(gcf, 'layers.png', 'Resolution', 300)











%% mean beta effect only for t
% Category model > 1.3238
mean(t(13:29, 6:11), 'all')

for layi = 1:56
   
    t = squeeze(tL(layi,:,:));
    tDiff(layi) = mean(t(13:29, 6:11), 'all') - 1.3238;
     
end

lois = (1:8:56) + 7
figure()
plot(tDiff(lois), 'LineWidth', 2)
set(gca, 'FontSize', 20, 'ylim', [-.3 .8])
exportgraphics(gcf, 'layers.png', 'Resolution', 300)








%% plot observed 3D clusters

h3D2 = h3D; 
h3D2 = permute(h3D2, [3 1 2]); 

xslice = [0:51];                               % define the cross sections to view
yslice = 0;
zslice = 0;

x = 1:1:8;
y = 1:1:45;
z = 1:1:52;
[X,Y,Z] = meshgrid(x,y,z);

figure(); set(gcf, 'Position', [100 100 400 400]); hold on; 

%slice(X, Y, Z, h3D2, xslice, yslice, zslice)    % display the slices
%contourslice(X, Y, Z, h3D2, xslice, yslice, zslice, 1)    % display the slices

% % % red 3D shape
 fv = isosurface(X, Y, Z, h3D2); 
 p1 = patch(fv);
 p1.FaceColor = 'red';
 p1.EdgeColor = [0.5 0.5 0.5];
 daspect([1 1 1]); 



% % % leave this for the grey plane at 0
 x = 0:.01:54;
 y =  ones(size(x))*6;
 z = ones(size(x))*8;
 col = y;  % This is the color, vary with x in this case.
 p2 = patch([zeros(size(x));z],[y;y],[x;x],[col;col]);
 set(p2,'edgealpha',0.01)

xlim([1 8])
ylim([1 54])
zlim([1 45])
view(34,24)

set(gca, 'ytick', [6 54], 'yticklabel', {'0' '4'} )
set(gca, 'ztick', [1 45], 'zticklabel', {'1' '150'}, 'FontSize', 18, 'ylim', [1 25])

exportgraphics(gcf, 'clust3D.png', 'Resolution', 150); 


%% plot observed 3D clusters RNN

for timei = 0:7

lois = (1:8:56) +timei

h3D2 = h3D; 
h3D2 = permute(h3D2, [3 1 2]); 

xslice = [0:7];                               % define the cross sections to view
yslice = 0;
zslice = 0;

x = 1:1:7;
y = 1:1:45;
z = 1:1:52;
[X,Y,Z] = meshgrid(x,y,z);

figure(); set(gcf, 'Position', [100 100 400 400]); hold on; 

%slice(X, Y, Z, h3D2, xslice, yslice, zslice)    % display the slices
%contourslice(X, Y, Z, h3D2, xslice, yslice, zslice, 1)    % display the slices

% % % red 3D shape
 fv = isosurface(X, Y, Z, h3D2(:,lois,:)); 
 p1 = patch(fv);
 p1.FaceColor = 'red';
 p1.EdgeColor = [0.5 0.5 0.5];
 daspect([1 1 1]); 



% % % leave this for the grey plane at 0
 x = 0:.01:52;
 y =  ones(size(x))*6;
 z = ones(size(x))*8;
 col = y;  % This is the color, vary with x in this case.
 p2 = patch([zeros(size(x));z],[y;y],[x;x],[col;col]);
 set(p2,'edgealpha',0.01)

xlim([1 8])
ylim([1 54])
zlim([1 45])
view(34,24)

set(gca, 'ytick', [6 54], 'yticklabel', {'0' '4'} )
set(gca, 'ztick', [1 45], 'zticklabel', {'1' '150'}, 'FontSize', 18, 'ylim', [1 25])

filename = [num2str(timei, '%02.f') 'clust3D.png']
exportgraphics(gcf, filename, 'Resolution', 150); 

end

%% plot only last time point one line

figure(); set(gcf, 'Position', [100 100 1500 1000]);
lois  = [8 16 24 32 40 48 54]; 

myCmap = colormap(brewermap([],'YlOrRd'));
for layi = 1:7
    subplot(1,7,layi)
    %times = 0:.01:3.49;
    times = 0:.01:4.49;
    %freqs = [.1:.1:29.9 30:.5:150]; 
    freqs = 1:540;
    contourf(times, freqs, myresizem(t, 10), 100, 'linecolor', 'none'); hold on; 
    contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 1); %colorbar
    plot([10 10], get(gca, 'ylim'), 'k:', 'LineWidth', 2);
    %plot(get(gca, 'xlim'),[300 300], 'k:', 'LineWidth', 2);
    set(gca, 'xtick', [-0.055 .945 1.945 2.945 3.945 ], 'xticklabel', [0 1 2 3 4], 'xlim', [0 3.5], 'clim', [-7 7])
    set(gca,'XTick',[], 'YTick', [])
    set(gca,'xticklabel',[], 'FontSize', 12)
    
end

ha=get(gcf,'children');
n = .1; 
set(ha(1),'position',[ 6/9 0 n n ])
set(ha(2),'position',[ 5/9 0 n n ])
set(ha(3),'position',[ 4/9 0 n n])
set(ha(4),'position',[ 3/9 0 n n ])
set(ha(5),'position',[ 2/9 0 n n ])
set(ha(6),'position',[ 1/9 0 n n])
set(ha(7),'position',[ 0 0 n n ])
exportgraphics(gcf, 'allM.png', 'Resolution', 300);




%% plot locked to cue

allT_cue = squeeze(mean(all_r_Times(:, 8, :, 6:10), 4, 'omitnan'));
plot(mean(allT_cue, 'omitnan')); hold on; 
size(allT_cue)
[h p ci ts] = ttest(allT_cue)
bL = h; bL(bL==0) = nan; bL(bL == 1) = 0; 





%% Plot in 3 different time periods

nLays = 56; 

%allTC= squeeze(mean(all_r_Times(:, :, :, :), 4, 'omitnan'));
allTC= squeeze(mean(all_r_Times, 4, 'omitnan'));
size(allTC)
[h p ci ts] = ttest(allTC);
t = squeeze(ts.tstat); h = squeeze(h);


figure()
subplot(211)
imagesc(flipud(squeeze(mean(allTC, 'omitnan')))); colorbar; hold on; 
plot (get(gca, 'ylim'), [8+.5 8+.5], 'LineWidth', 1,'Color', [0 0 0 ] );
plot (get(gca, 'ylim'), [16+.5 16+.5], 'LineWidth', 1,'Color', [0 0 0 ] );
plot (get(gca, 'ylim'), [24+.5 24+.5], 'LineWidth', 1,'Color', [0 0 0 ] );
plot (get(gca, 'ylim'), [32+.5 32+.5], 'LineWidth', 1,'Color', [0 0 0 ] );
plot (get(gca, 'ylim'), [40+.5 40+.5], 'LineWidth', 1,'Color', [0 0 0 ] );
plot (get(gca, 'ylim'), [48+.5 48+.5], 'LineWidth', 1,'Color', [0 0 0 ] );
set(gca, 'xtick', [], 'xticklabels', [], 'ytick', [], 'yticklabels', [])
plot ([30 30],get(gca, 'xlim'),  'k:', 'LineWidth', 2);

subplot(212)
imagesc(flipud(myresizem(t, 10))); colorbar; hold on;
contour(flipud(myresizem(h, 10)), 1, 'Color', [0, 0, 0], 'LineWidth', 1);
plot (get(gca, 'ylim'), [80 80], 'LineWidth', 1,'Color', [0 0 0 ] );
plot (get(gca, 'ylim'), [160 160], 'LineWidth', 1,'Color', [0 0 0 ] );
plot (get(gca, 'ylim'), [240 240], 'LineWidth', 1,'Color', [0 0 0 ] );
plot (get(gca, 'ylim'), [320 320], 'LineWidth', 1,'Color', [0 0 0 ] );
plot (get(gca, 'ylim'), [400 400], 'LineWidth', 1,'Color', [0 0 0 ] );
plot (get(gca, 'ylim'), [480 480], 'LineWidth', 1,'Color', [0 0 0 ] );
set(gca, 'xtick', [], 'xticklabels', [], 'ytick', [], 'yticklabels', [])
plot ([300 300],get(gca, 'xlim'),  'k:', 'LineWidth', 2);









%% Plot in 3 different time periods vertical

nLays = 56; 

allTC= squeeze(mean(all_r_Times, 4, 'omitnan'));
[h p ci ts] = ttest(allTC);
t = squeeze(ts.tstat); h = squeeze(h);

mAllTC = squeeze(mean(allTC, 'omitnan'))';


figure()
subplot(211)
imagesc(flipud(myresizem(mAllTC, 10))); colorbar; hold on; 
plot ([80 80], get(gca, 'xlim'), 'LineWidth', 1,'Color', [0 0 0 ] );
plot ([160 160],get(gca, 'xlim'),  'LineWidth', 1,'Color', [0 0 0 ] );
plot ([240 240], get(gca, 'xlim'), 'LineWidth', 1,'Color', [0 0 0 ] );
plot ([320 320], get(gca, 'xlim'), 'LineWidth', 1,'Color', [0 0 0 ] );
plot ([400 400], get(gca, 'xlim'), 'LineWidth', 1,'Color', [0 0 0 ] );
plot ([480 480], get(gca, 'xlim'), 'LineWidth', 1,'Color', [0 0 0 ] );
set(gca, 'xtick', [], 'xticklabels', [], 'ytick', [], 'yticklabels', [])
plot ([0 560],  [540-300 540-300],'k:', 'LineWidth', 2);axis equal
set(gca, 'xlim', [0 560])

subplot(212)
imagesc(flipud(myresizem(t', 10))); colorbar; hold on;
contour(flipud(myresizem(h', 10)), 1, 'Color', [0, 0, 0], 'LineWidth', 1);
plot ([80 80], get(gca, 'xlim'), 'LineWidth', 1,'Color', [0 0 0 ] );
plot ([160 160],get(gca, 'xlim'),  'LineWidth', 1,'Color', [0 0 0 ] );
plot ([240 240], get(gca, 'xlim'), 'LineWidth', 1,'Color', [0 0 0 ] );
plot ([320 320], get(gca, 'xlim'), 'LineWidth', 1,'Color', [0 0 0 ] );
plot ([400 400], get(gca, 'xlim'), 'LineWidth', 1,'Color', [0 0 0 ] );
plot ([480 480], get(gca, 'xlim'), 'LineWidth', 1,'Color', [0 0 0 ] );
set(gca, 'xtick', [], 'xticklabels', [], 'ytick', [], 'yticklabels', [])
plot ([0 560],  [540-300 540-300],'k:', 'LineWidth', 2);axis equal







%% Plot in 3 different time periods vertical Alexnet

allTC= squeeze(mean(all_r_Times, 4, 'omitnan'));
[h p ci ts] = ttest(allTC);
t = squeeze(ts.tstat); h = squeeze(h);
%h = zeros(8, 54); %delete outline

clustinfo = bwconncomp(h);
c2r = [ 1 2 4:7  ]; % cue pfc% cue pfc
%c2r = [1:8]; % enc pfc
%c2r = [1:9]% maint pfc
%c2r = [1 3:6]; % Maint VVS
%c2r = [1:7]; % Cue VVS
%c2r = [1]; % Enc VVS
for i = 1:length(c2r)
    h(clustinfo.PixelIdxList{c2r(i)})= 0;
end

mAllTC = squeeze(mean(allTC, 'omitnan'))';



figure(); set(gcf, 'Position', [0 0 10 1300])
%subplot(211)
%imagesc(flipud(myresizem(mAllTC, 10))); hold on; 
%set(gca, 'xtick', [], 'xticklabels', [], 'ytick', [], 'yticklabels', [])
%plot ([0 560],  [540-300 540-300],'k:', 'LineWidth', 3);
%set(gca, 'xlim', [0 80])
%set(gca, 'xtick', [], 'xticklabels', [], 'ytick', [], 'yticklabels', [])

subplot(212)
imagesc(flipud(myresizem(t', 10))); hold on;
contour(flipud(myresizem(h', 10)), 1, 'Color', [0, 0, 0], 'LineWidth', 3);
plot ([0 560],  [540-300 540-300],'k:', 'LineWidth', 5);
set(gca, 'xlim', [0 80])

set(gca, 'xtick', [], 'xticklabels', [], 'ytick', [], 'yticklabels', [], 'clim', [-7 7])

myCmap = colormap(brewermap([],'YlOrRd'));
colormap(myCmap)
exportgraphics(gcf, 'vertical_plots.png', 'Resolution', 300); 

%% Plot in 3 different time periods vertical RNN

%lois = [8 16 24 32 40 48 56]; 
lois = (1:8:56) + 7


allTC= squeeze(mean(all_r_Times(:, lois,:,:), 4, 'omitnan'));
[h p ci ts] = ttest(allTC);
t = squeeze(ts.tstat); h = squeeze(h);
%h = zeros(8, 54); %delete outline

clustinfo = bwconncomp(h);
%c2r = [ 1 2 4:7  ]; % cue pfc
%c2r = [1:7]; % enc pfc
%c2r = [1:3]% maint pfc
%c2r = [2:7]; % Maint VVS
%c2r = [1:7]; % Cue VVS
%c2r = [1]; % Enc VVS
c2r = []; % Empty
for i = 1:length(c2r)
    h(clustinfo.PixelIdxList{c2r(i)})= 0;
end

mAllTC = squeeze(mean(allTC, 'omitnan'))';



figure(); set(gcf, 'Position', [0 0 10 1300])
%subplot(211)
%imagesc(flipud(myresizem(mAllTC, 10))); hold on; 
%set(gca, 'xtick', [], 'xticklabels', [], 'ytick', [], 'yticklabels', [])
%plot ([0 560],  [540-300 540-300],'k:', 'LineWidth', 3);
%set(gca, 'xlim', [0 80])
%set(gca, 'xtick', [], 'xticklabels', [], 'ytick', [], 'yticklabels', [])

subplot(212)
imagesc(flipud(myresizem(t', 10))); hold on;
contour(flipud(myresizem(h', 10)), 1, 'Color', [0, 0, 0], 'LineWidth', 3);
plot ([0 560],  [540-300 540-300],'k:', 'LineWidth', 5);
set(gca, 'xlim', [0 70])

set(gca, 'xtick', [], 'xticklabels', [], 'ytick', [], 'yticklabels', [], 'clim', [-7 7])

myCmap = colormap(brewermap([],'YlOrRd'));
colormap(myCmap)
exportgraphics(gcf, 'vertical_plots.png', 'Resolution', 300); 


%% FREQUENCY RESOLVED DNN ANALYSIS PERMUTATIONS

clearvars -except act_CH act_FR 

%lois = (1:8:56) +0; 
lois = 56; % use for alexnet 
nPerm =     100; 
f2sav       = 'RNN_pfc_M&C_noAv_54_Perm_nT_1-49_C_naFV'
load_parameters_WMFRA;
loadNet_WM;


all_r_Times_perm    = zeros(nPerm, nSubj, nLays, nFreqs, nTimes);

tic

f2u = strsplit(f2sav, '_');
for subji=1:nSubj
    disp (['Subj: ' num2str(subji)])
    load(sublist{subji});
    cfg_contrasts.oneListPow = cfg_contrasts.oneListPow(:, :, :, resTime); 
    cfg_contrasts = normalize_WM(cfg_contrasts, 1, 'blo', [3 7]);
    if strcmp(avMeth,'pow')
       cfg_contrasts = average_repetitions(cfg_contrasts);
    end
    
    for freqi = 1:length(freqs2test)
        f  = freqs2test(freqi);
        if strcmp(f2u{10}, 'aFV')
            avTW = 1;
        else
            avTW = 0;
        end
        rdm_prev(freqi) = create_rdms(cfg_contrasts, f, it, 5, 1, avTW);
    end
    
   [allS ids act_FR2 act_CH2 cM] = getIDsAndnRDMs(rdm_prev, act_CH, act_FR);
    
    for permi = 1:nPerm
        idp = randperm(length(ids)); %this has to be done at each permutation but not each layer
        parfor layi = 1:length(lois)
            loi = lois(layi);
            if subji < subj_ch_fr
                M =  squeeze(act_FR2(loi,:,:)); 
            else
                M =  squeeze(act_CH2(loi,:,:)); 
            end
            Mp = M(idp, idp);
            Mp(Mp ==1) = 1000;Mp(Mp ==0) = 2000;
            Mp = tril(Mp, -1);Mp(Mp==0) =[];
            Mp(Mp==1000) = 1;Mp(Mp==2000) = 0;
            
            cMh = cM; 
            cMh = cMh(idp, idp);
            cMh(cMh ==1) = 1000;cMh(cMh ==0) = 2000;
            cMh = tril(cMh, -1);cMh(cMh==0) =[];
            cMh(cMh==1000) = 1;cMh(cMh==2000) = 0;

            
            if strcmp(f2u{9}, 'PC')
                allTEst = partialcorri(allS, Mp', cMh', 'type', 's');    
            else
               allTEst = corr(allS, Mp', 'type', 's'); 
            end
            
            all_r_Times_perm(permi, subji, layi,:, : ) = reshape(allTEst, nFreqs, nTimes );  
        end
    end   
end 

cd (currentF) 

f2sav = [f2sav '_' num2str(nPerm) 'p']; 
save(f2sav, 'all_r_Times_perm', '-v7.3');

toc





%% get obs data
clear allSTs max_clust_sum_obs
clustinfo = bwconncomp(h);
if length(clustinfo.PixelIdxList) > 0
    for pxi = 1:length(clustinfo.PixelIdxList)
       allSTs(pxi,:)= sum(t(clustinfo.PixelIdxList{pxi}));% it is maxi -1 cause we are adding a zero below
    end
else
    allSTs = 0;
end

[maxh id] = max(abs(allSTs));
max_clust_sum_obs = allSTs(id); 

%%  count permutations

nPerm = 1000;
tR = 1:35; %time range 

sub2exc =  [1 3 4 7 10 13 16];

clear max_clust_sum_perm
for permi = 1:nPerm

allTC= squeeze(all_r_Times_perm(permi, :, 1, :, tR));

allTC(sub2exc, :, :) = []; 

[h p ci ts] = ttest(allTC);
t = squeeze(ts.tstat); h = squeeze(h);

clear allSTs 
clustinfo = bwconncomp(h);
if length(clustinfo.PixelIdxList) > 0
    for pxi = 1:length(clustinfo.PixelIdxList)
       allSTs(pxi,:)= sum(t(clustinfo.PixelIdxList{pxi}));% it is maxi -1 cause we are adding a zero below
    end
else
    allSTs = 0;
end

[maxh id] = max(abs(allSTs));
max_clust_sum_perm(permi, :) = allSTs(id); 



end



%% %%
t1 = max_clust_sum_obs; 
%t1 = [120.384014285890] % pfc beta effect, 56 > last layer and time point for no average data 
%t1 = 1.117001518186384e+02 % pfc bea effect, 56 > last layer and time point for average data
%t1 = 59.713604023815556 % pfc beta effect, 56 > last layer and time point for no average data with partially correcting for category model
%t1 = 51.086190103391460 %hippocampus, encoding average 
%t1 = 59.007; %pfc beta effect last layer and time point when not averaging FV in time

[counts,centers] = hist(max_clust_sum_perm, 14);

figure()
h = bar(centers, counts); hold on;
h.FaceColor = 'w';
h.EdgeColor = 'k';
h.LineWidth =2;
set(gca, 'FontSize', 20, 'ylim', [0 350] );
plot ([t1 t1], get(gca, 'ylim'), 'LineWidth', 3,'Color', [.1 .3 .6] );
filename = ['histogram.png'];
%export_fig(1, filename,'-transparent', '-r300');
%close all;

%% get rank for stats

%mcsR = max_clust_sum_obs; 
mcsR = t1; 
mcsP = max_clust_sum_perm; 
allAb = mcsP(abs(mcsP) > abs(mcsR));
p = 1 - ((nPerm-1) - (length (allAb)))  / nPerm







%% FREQUENCY RESOLVED DNN ANALYSIS FOR 7 LAYERS ONLY and not across time

clearvars -except act_CH act_FR 

%lois         = [1 9 17 25 33 41 49]; 
lois         = [8 16 24 32 40 48 56]; 
%lois           = 1:8; %use this for alexnet
f2sav       = 'RNN_vvs_CueAll_noAv_54_Real_nT_6-14_AvRDM'
load_parameters_WMFRA;
loadNet_WM;

all_r_Times    = zeros(nSubj, length(lois), nFreqs);

tic

for subji=1:nSubj
    disp (['Subj: ' num2str(subji)])
    load(sublist{subji});
    cfg_contrasts.oneListPow = cfg_contrasts.oneListPow(:, :, :, resTime); 
    cfg_contrasts = normalize_WM(cfg_contrasts, 1, 'blo', [3 7]);
    if strcmp(avMeth,'pow')
       cfg_contrasts = average_repetitions(cfg_contrasts);
    end
    
    for freqi = 1:length(freqs2test)
        f  = freqs2test(freqi);
        rdm_prev_prev = create_rdms(cfg_contrasts, f, it, 5, 1);
        % average RDMS across time
        rdm_prev_prev.rdm = squeeze(mean(rdm_prev_prev.rdm,3));
        rdm_prev(freqi) = rdm_prev_prev; 
    end
    
   [allS ids{subji} act_FR2 act_CH2] = getIDsAndnRDMs(rdm_prev, act_CH, act_FR);
    
    parfor layi = 1:length(lois)
        loi = lois(layi);
        if subji < subj_ch_fr
            M =  squeeze(act_FR2(loi,:,:)); 
        else
            M =  squeeze(act_CH2(loi,:,:)); 
        end
        M(M ==1) = 1000;M(M ==0) = 2000;
        M = tril(M, -1);M(M==0) =[];
        M(M==1000) = 1;M(M==2000) = 0;

        allTEst = corr(allS, M', 'type', 's');

        all_r_Times(subji, layi,:) = allTEst;  
    end
 
end 

cd (currentF) 
toc 
save(f2sav, 'all_r_Times', 'ids');






%% Plot in 3 different time periods vertical Alexnet

allTC= squeeze(mean(all_r_Times, 4, 'omitnan'));
[h p ci ts] = ttest(allTC);
t = squeeze(ts.tstat); h = squeeze(h);
%h = zeros(8, 54); %delete outline

clustinfo = bwconncomp(h);
%c2r = [1:4 6:11]; % cue pfc
%c2r = [1:8]; % enc pfc
%c2r = [2:7]% maint pfc
%c2r = [1 3:6]; % Maint VVS
%c2r = [1:7]; % Cue VVS
%c2r = [1:2 4:7]; % Maint VVS
c2r = [1 3:7]; 
for i = 1:length(c2r)
    h(clustinfo.PixelIdxList{c2r(i)})= 0;
end

mAllTC = squeeze(mean(allTC, 'omitnan'))';



figure(); set(gcf, 'Position', [0 0 10 1300])
%subplot(211)
%imagesc(flipud(myresizem(mAllTC, 10))); hold on; 
%set(gca, 'xtick', [], 'xticklabels', [], 'ytick', [], 'yticklabels', [])
%plot ([0 560],  [540-300 540-300],'k:', 'LineWidth', 3);
%set(gca, 'xlim', [0 80])
%set(gca, 'xtick', [], 'xticklabels', [], 'ytick', [], 'yticklabels', [])

subplot(212)
imagesc(flipud(myresizem(t', 10))); hold on;
contour(flipud(myresizem(h', 10)), 1, 'Color', [0, 0, 0], 'LineWidth', 3);
plot ([0 560],  [540-300 540-300],'k:', 'LineWidth', 5);
set(gca, 'xlim', [0 80])

set(gca, 'xtick', [], 'xticklabels', [], 'ytick', [], 'yticklabels', [], 'clim', [-7 7])

myCmap = colormap(brewermap([],'YlOrRd'));
colormap(myCmap)
exportgraphics(gcf, 'vertical_plots.png', 'Resolution', 300); 

%% 
clear allSTs max_clust_sum_obs
clustinfo = bwconncomp(h);
if length(clustinfo.PixelIdxList) > 0
    for pxi = 1:length(clustinfo.PixelIdxList)
       allSTs(pxi,:)= sum(t(clustinfo.PixelIdxList{pxi}));% it is maxi -1 cause we are adding a zero below
    end
else
    allSTs = 0;
end

[maxh id] = max(abs(allSTs));
max_clust_sum_obs = allSTs(id); 





%% Plot in 3 different time periods vertical RNN


lois = 1:7


allTC= squeeze(mean(all_r_Times(:, lois,:,:), 4, 'omitnan'));
[h p ci ts] = ttest(allTC);
t = squeeze(ts.tstat); h = squeeze(h);
%h = zeros(8, 54); %delete outline

clustinfo = bwconncomp(h);
%c2r = [ 2:5  ]; % cue pfc
%c2r = [1:7]; % enc pfc
%c2r = [1:3]% maint pfc
%c2r = [1 3:7]; % Maint VVS
%c2r = [1:7]; % Cue VVS
%c2r = [1]; % Enc VVS
c2r = [1 2]; % Empty
for i = 1:length(c2r)
    h(clustinfo.PixelIdxList{c2r(i)})= 0;
end

mAllTC = squeeze(mean(allTC, 'omitnan'))';



figure(); set(gcf, 'Position', [0 0 10 1300])
%subplot(211)
%imagesc(flipud(myresizem(mAllTC, 10))); hold on; 
%set(gca, 'xtick', [], 'xticklabels', [], 'ytick', [], 'yticklabels', [])
%plot ([0 560],  [540-300 540-300],'k:', 'LineWidth', 3);
%set(gca, 'xlim', [0 80])
%set(gca, 'xtick', [], 'xticklabels', [], 'ytick', [], 'yticklabels', [])

subplot(212)
imagesc(flipud(myresizem(t', 10))); hold on;
contour(flipud(myresizem(h', 10)), 1, 'Color', [0, 0, 0], 'LineWidth', 3);
plot ([0 560],  [540-300 540-300],'k:', 'LineWidth', 5);
set(gca, 'xlim', [0 70])

set(gca, 'xtick', [], 'xticklabels', [], 'ytick', [], 'yticklabels', [], 'clim', [-7 7])

myCmap = colormap(brewermap([],'YlOrRd'));
colormap(myCmap)
exportgraphics(gcf, 'vertical_plots.png', 'Resolution', 300); 

%% Plot in 3 different time periods vertical separately for first and last time points in each layer

%allTC= squeeze(mean(all_r_Times, 4, 'omitnan'));
allTC= all_r_Times;
[h p ci ts] = ttest(allTC);
t = squeeze(ts.tstat); h = squeeze(h);

mAllTC = squeeze(mean(allTC, 'omitnan'))';



figure(); set(gcf, 'Position', [100 100 100 500])
subplot(211)
imagesc(flipud(myresizem(mAllTC, 10))); colorbar; hold on; 
%set(gca, 'xtick', [], 'xticklabels', [], 'ytick', [], 'yticklabels', [])
plot ([0 560],  [540-300 540-300],'k:', 'LineWidth', 2);
set(gca, 'xlim', [0 70])

subplot(212)
imagesc(flipud(myresizem(t', 10))); colorbar; hold on;
contour(flipud(myresizem(h', 10)), 1, 'Color', [0, 0, 0], 'LineWidth', 1);
plot ([0 560],  [540-300 540-300],'k:', 'LineWidth', 2);
set(gca, 'xlim', [0 70])




exportgraphics(gcf, 'vertical_plots.png', 'Resolution', 300); 



%% FREQUENCY RESOLVED DNN ANALYSIS FOR 7 LAYERS ONLY and not across time PERMUTATIONS IN LOOP 

clearvars -except act_CH act_FR 

nPerm         = 100;
%lois         = [1 9 17 25 33 41 49]; 
lois         = [8 16 24 32 40 48 56]; 
%lois          = 1:8; %use this for alexnet


listF2sav = {
             %'Alex_pfc_Enc_Av_54_Perm_nT_6-17'  'Alex_vvs_Enc_Av_54_Perm_nT_6-17' ...   
             %'Alex_pfc_Cue_Av_54_Perm_nT_6-14' ...
             %'Alex_vvs_Cue_Av_54_Perm_nT_6-14'    ...
             %'Alex_pfc_Maint_Av_54_Perm_nT_11-49' 
             %'Alex_vvs_Maint_Av_54_Perm_nT_11-49' ...
             %'RNN_pfc_Enc_Av_54_Perm_nT_6-17'  'RNN_vvs_Enc_Av_54_Perm_nT_6-17' ...   
             %'RNN_pfc_Cue_Av_54_Perm_nT_6-14' 
             %'RNN_vvs_Cue_Av_54_Perm_nT_6-14'    ...
             %'RNN_pfc_Maint_Av_54_Perm_nT_11-49' 'RNN_vvs_Maint_Av_54_Perm_nT_11-49' 
             'RNN_vvs_MaintAll_noAv_54_Perm_nT_11-49_AvRDM'
             };
    

for listi = 1:length(listF2sav)
    clearvars -except act_CH act_FR listF2sav listi nPerm lois
    
    f2sav       = listF2sav{listi}
    load_parameters_WMFRA;
    loadNet_WM;
    all_r_Times_perm    = zeros(nPerm, nSubj, length(lois), nFreqs);
    
    for subji=1:nSubj
        disp (['Subj: ' num2str(subji)])
        load(sublist{subji});
        cfg_contrasts.oneListPow = cfg_contrasts.oneListPow(:, :, :, resTime); 
        cfg_contrasts = normalize_WM(cfg_contrasts, 1, 'blo', [3 7]);
        if strcmp(avMeth,'pow')
           cfg_contrasts = average_repetitions(cfg_contrasts);
        end

        for freqi = 1:length(freqs2test)
            f  = freqs2test(freqi);
            rdm_prev_prev = create_rdms(cfg_contrasts, f, it, 5, 1);
            % average RDMS across time
            rdm_prev_prev.rdm = squeeze(mean(rdm_prev_prev.rdm,3));
            rdm_prev(freqi) = rdm_prev_prev; 
        end

       [allS ids act_FR2 act_CH2] = getIDsAndnRDMs(rdm_prev, act_CH, act_FR);

       for permi = 1:nPerm
           idp = randperm(length(ids)); %this has to be done at each permutation but not each layer
            parfor layi = 1:length(lois)
                loi = lois(layi);
                if subji < subj_ch_fr
                    M =  squeeze(act_FR2(loi,:,:)); 
                else
                    M =  squeeze(act_CH2(loi,:,:)); 
                end
                Mp = M(idp, idp);
                Mp(Mp ==1) = 1000;Mp(Mp ==0) = 2000;
                Mp = tril(Mp, -1);Mp(Mp==0) =[];
                Mp(Mp==1000) = 1;Mp(Mp==2000) = 0;

                allTEst = corr(allS, Mp', 'type', 's');

                all_r_Times_perm(permi, subji, layi,:) = allTEst;  
            end
       end

    end 

    cd (currentF) 
    toc 
    save(f2sav, 'all_r_Times_perm', 'ids');

end






%%  count permutations

nPerm = 100;

clear max_clust_sum_perm
for permi = 1:nPerm

allTC= squeeze(mean(all_r_Times_perm(permi, :, :, :, :), 5, 'omitnan'));
[h p ci ts] = ttest(allTC);
t = squeeze(ts.tstat); h = squeeze(h);

clear allSTs 
clustinfo = bwconncomp(h);
if length(clustinfo.PixelIdxList) > 0
    for pxi = 1:length(clustinfo.PixelIdxList)
       allSTs(pxi,:)= sum(t(clustinfo.PixelIdxList{pxi}));% it is maxi -1 cause we are adding a zero below
    end
else
    allSTs = 0;
end

[maxh id] = max(abs(allSTs));
max_clust_sum_perm(permi, :) = allSTs(id); 



end


%% %%
%[67.8299432512560;73.1548118121768;55.2789515185440] % 1st 
%[60.747837472955810;65.541581939972790;45.307019796192414] % 8th
%t1 = max_clust_sum_obs; 
t1 = 77.516198717576510; %pfc cue 8th
%t1 = [58.0412563673653]; %pfc cue Alex
%t1 = [82.3523945322749]; %vvs maint alex no av
t1 = [63.7020746868958]
[counts,centers] = hist(max_clust_sum_perm, 14);

figure()
h = bar(centers, counts); hold on;
h.FaceColor = 'w';
h.EdgeColor = 'k';
h.LineWidth =2;
set(gca, 'FontSize', 20, 'ylim', [0 350] );
plot ([t1 t1], get(gca, 'ylim'), 'LineWidth', 3,'Color', [.1 .3 .6] );
filename = ['histogram.png'];
%export_fig(1, filename,'-transparent', '-r300');
%close all;

%% get rank for stats

%mcsR = max_clust_sum_obs; 
mcsR = t1; 
mcsP = max_clust_sum_perm; 
allAb = mcsP(abs(mcsP) > abs(mcsR));
p = 1 - ((nPerm-1) - (length (allAb)))  / nPerm


%% Plot wihout non significant clusters

allTC= squeeze(mean(all_r_Times, 4, 'omitnan'));
[h p ci ts] = ttest(allTC);
t = squeeze(ts.tstat); h = squeeze(h);
clustinfo = bwconncomp(h);
%c2r = [1 3 5 6] % 1st vvs
%c2r = [3 4] % 8th vvs
%c2r = [2 3] % 1st and 8th (coincidence) pfc
c2r = 1:length(clustinfo.PixelIdxList); % delete all
for i=1:length(c2r)
    
    h(clustinfo.PixelIdxList{c2r(i)})= 0;
end


figure(); set(gcf, 'Position', [100 100 100 500])
subplot(212)
imagesc(flipud(myresizem(t', 10))); colorbar; hold on;
contour(flipud(myresizem(h', 10)), 1, 'Color', [0, 0, 0], 'LineWidth', 1);
plot ([0 560],  [540-300 540-300],'k:', 'LineWidth', 2);
set(gca, 'xlim', [0 70])



exportgraphics(gcf, 'vertical_plots.png', 'Resolution', 300); 


%% Blue plots > average full maint period per layer for 1 or 8 layer

lois         = [1 9 17 25 33 41 49]; 
%lois         = [8 16 24 32 40 48 54]; 

subj2exc = []; 
nLays = length(lois); 
nTimes = size(all_r_Times, 4); 
all_r_TM = squeeze(mean(all_r_Times(:, lois, :, :), 4));

figure()
mART = squeeze(mean(all_r_TM, 'omitnan')); 
stdART = std(all_r_TM, 'omitnan'); 
seART = stdART/ sqrt(size(all_r_TM, 1));
[h p ci t] = ttest(all_r_TM); 
hL = h; hL(hL == 0) = nan; hL(hL==1) = 0; 

shadedErrorBar(1:nLays, mART, seART, 'r', 1); hold on; 
plot(1:nLays, hL, 'LineWidth', 6)
%plot (allRTMA); hold on; 
set(gca, 'FontSize', 12)



%exportgraphics(gcf, 'figure.png', 'Resolution', 300)

%%Bar format
data.data = [all_r_TM]; 

figure(2); set(gcf,'Position', [0 0 1200 1000]); 
mean_S = mean(data.data, 1);
h = bar (mean_S);hold on;
hb = plot ([1:nLays], data.data); hold on;
set(hb, 'lineWidth', 3, 'Marker', '.', 'MarkerSize',20);hold on;
%h1=errorbar(mean_S, se_S,'c'); hold on;
set(h,'facecolor', 'flat', 'lineWidth', 1);
set(hb, 'Color','k','linestyle','none', 'lineWidth', 2);
set(gca,'XTick',[1:nLays],'XTickLabel',{'', ''}, ...
    'FontSize', 30, 'linew',1, 'xlim', [0 nLays+1], 'ylim', [-.05 .15] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 1);
scatter(1:nLays, hL+0.135, 'k', 'LineWidth', 6, 'Marker', '*')


[h p ci t] = ttest (data.data);
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth', 1);

export_fig(2, '_2.png','-transparent', '-r300');
close all;   


%% 
clear allSTs max_clust_sum_obs
clustinfo = bwconncomp(h);
if length(clustinfo.PixelIdxList) > 0
    for pxi = 1:length(clustinfo.PixelIdxList)
       allSTs(pxi,:)= sum(t(clustinfo.PixelIdxList{pxi}));% it is maxi -1 cause we are adding a zero below
    end
else
    allSTs = 0;
end

[maxh id] = max(abs(allSTs));
max_clust_sum_obs = allSTs(id); 






%% FREQUENCY RESOLVED DNN ANALYSIS TRIALS

clearvars -except act_CH act_FR 
f2sav       = 'RNN_vvs_M&C_noAv_54_Real_T_1-49_C_aFV'; %'Alex_vvs_CueAll_noAv_54_Real_nT_1-49'
load_parameters_WMFRA;
loadNet_WM;

all_r_Times    = zeros(nSubj, nLays, nFreqs, nTimes);

tic
f2u = strsplit(f2sav, '_');
for subji=1:nSubj
    disp (['Subj: ' num2str(subji)])
    load(sublist{subji});
    cfg_contrasts.oneListPow = cfg_contrasts.oneListPow(:, :, :, resTime); 
    cfg_contrasts = normalize_WM(cfg_contrasts, 1, 'blo', [3 7]);
    if strcmp(avMeth,'pow')
       cfg_contrasts = average_repetitions(cfg_contrasts);
    end
    
    for freqi = 1:length(freqs2test)
        f  = freqs2test(freqi);
        if strcmp(f2u{10}, 'aFV')
            avTW = 1;
        else
            avTW = 0;
        end
        rdm_prev(freqi) = create_rdms(cfg_contrasts, f, it, 5, 1, avTW);%
        
    end
    
    ids = rdm_prev(1).ids; % take first freq only to get ids
    ids = cellfun(@(x) strsplit(string(x)), ids, 'un', 0);
    ids0 = cellfun(@(x) x(3), ids, 'un', 0); ids0 = double(string(ids0)); 
    ids1 = char(string((ids0))); ids2 = ids1(:,[1 3]);ids3 = double(string(ids2));
    idx = find(~mod(ids3, 10)); ids4 = ids3-10; ids4(idx) = ids3(idx);

    
    act_FR_2 = act_FR(:,ids4,ids4);
    act_CH_2 = act_CH(:,ids4,ids4);
    
    for triali = 1:size(act_FR_2, 2)
            
        allS = cat(4, rdm_prev.rdm);
        allS1 = squeeze(allS(triali, :, :)); 
        allS2 = allS1;allS2(allS2==1) = NaN;allS2 = squeeze(allS2);
        allS2(triali, :) = []; %remove first row (only has nans, cause comes from diagonal)
        

        for layi = 56:56 %1:nLays
            if subji < subj_ch_fr
                M =  squeeze(act_FR_2(layi,triali,:)); 
            else
                M =  squeeze(act_CH_2(layi,triali,:)); 
            end
            M(M ==1) = 1000;M(M ==0) = 2000;
            M = tril(M, -1);M(M==0) =[];
            M(M==1000) = 1;M(M==2000) = 0;

            allTEst = corr(allS2, M, 'type', 's');

            all_r_Times_Trials{subji}(layi,triali, : ,: ) = reshape(allTEst, nFreqs, nTimes );  
        end
    end
    
end 

cd (currentF) 
toc 
save(f2sav, 'all_r_Times_Trials');









%% plot and get cluster







%% build 2 compare resolved

v1 = all_r_Times_Trials_PFC([1:3 5 9:16])';
v2 = all_r_Times_Trials_VVS([6 7 9 13 18:23 27:28])';
vvs_pfc = [v1 v2];


%% compare resolved

clear cR
for subji = 1:12
   
    x = vvs_pfc{subji,1};
    y = vvs_pfc{subji,2};
    
    for layi = 1:56
        for freqi = 1:54
            x1 = squeeze(mean(x(layi, :, freqi, :), 4, 'omitnan')); 
            y1 = squeeze(mean(y(layi, :, freqi, :), 4, 'omitnan')); 
            cR(subji, layi, freqi, :) = corr(x1', y1', 'type', 's');
        end    
    end
end

%[h p ci ts] = ttest(cR);
%h = h'; t= ts.tstat'; 
%% 
d2p = squeeze(mean(cR, 'omitnan'));
[h p ci ts] = ttest(cR);
h= squeeze(h);

figure()
imagesc(d2p)
figure()
imagesc(h)







%%  Ttest against mean baseline and not against zero
subj2exc = []; 
nPerm = size(all_r_Times, 1);


for layi = 1:nLays
    d2u = squeeze(all_r_Times(:, layi, :, :));
    d2u(subj2exc, :, :) = []; 
    mBas = mean(d2u(:, :, 1:5), 3); 
    mBas = mean(mBas,2);

    clear hs p ts
    for freqi = 1:54
        for timei = 1:45
            [hs(freqi, timei) p ci tsprev] =   ttest(d2u(:, freqi, timei), mBas);  
            ts(freqi, timei) = tsprev.tstat;
        end
    end

    h3D(layi,:,:) = squeeze(hs); 
    t3D(layi,:,:) = squeeze(ts); 
end


% get 3D cluster
clear allSTs   
clustinfo = bwconncomp(h3D, 6);
if length(clustinfo.PixelIdxList) > 0
    for pxi = 1:length(clustinfo.PixelIdxList)
       allSTs(pxi) = sum(t3D(clustinfo.PixelIdxList{pxi}));% it is maxi -1 cause we are adding a zero below
    end
else
    allSTs(permi, :) = 0;
end
[maxh id] = max(abs(allSTs));
max_clust_sum_3D = allSTs(id); 
for ci = 1:length(clustinfo.PixelIdxList)
   modifiedCluster = abs(sum(t3D(clustinfo.PixelIdxList{ci})))+0.0001; %add this small number because removing negative values creates slighlly different valeus
   if modifiedCluster < max(abs(allSTs))  % + 1    %add +1 at the very end to delete all clusters
      %disp(['layi > ' num2str(layi) ' > ' num2str(abs(sum(tLay(clustinfo.PixelIdxList{ci})))) '_' num2str(max(abs(allSTs)))]) 
      h3D (clustinfo.PixelIdxList{ci}) =  0; 
   end
end






































%%  T-test against zero
subj2exc = [1]; 
for permi = 1:nPerm
    for layi = 1:nLays
        d2u = squeeze(all_r_Times_perm(permi,:, layi, :, :));
        d2u(subj2exc, :, :) = []; 
        [hs p ci ts] =   ttest(d2u);
        hPerm(permi, layi,:,:) = squeeze(hs); 
        tPerm(permi, layi,:,:) = squeeze(ts.tstat); 
    end
end  

%%  Ttest against mean baseline and not against zero
subj2exc = [1]; 
nPerm = size(all_r_Times_perm, 1);

for permi = 1:nPerm
    for layi = 51:51 %1:nLays
        d2u = squeeze(all_r_Times_perm(permi,:, layi, :, :));
        d2u(subj2exc, :, :) = []; 
        mBas = mean(d2u(:, :, 1:5), 3); 
        mBas = mean(mBas,2);
        
        clear hs p ts
        for freqi = 1:54
            for timei = 1:40
                [hs(freqi, timei) p ci tsprev] =   ttest(d2u(:, freqi, timei), mBas);  
                ts(freqi, timei) = tsprev.tstat;
            end
        end
        
        hPerm(permi, layi,:,:) = squeeze(hs); 
        tPerm(permi, layi,:,:) = squeeze(ts); 
    end
   
end  

%% check 1 permutation (onlyt t and h)

perm = 1; 
nLays = size(all_r_Times_perm, 3);
for layi = 1:nLays
    d2pT = squeeze(tPerm(perm, layi, :,:));
    d2pH = squeeze(hPerm(perm, layi, :,:));

    figure()
    imagesc(myresizem(d2pT, 10)); colorbar; hold on;
    contour(myresizem(d2pH, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 1);
end


%% count clusters
nPerm = size(tPerm, 1); 
f = 1:54;

nLays = size(all_r_Times_perm, 3);
clear max_clust_sum 
for permi = 1:nPerm
clear allSTs   
    for layi = 51:51 %1:nLays

        d2pT = squeeze(tPerm(permi, layi,f,:));
        d2pH = squeeze(hPerm(permi, layi,f,:));

        clustinfo = bwconncomp(d2pH);
        if length(clustinfo.PixelIdxList) > 0
            for pxi = 1:length(clustinfo.PixelIdxList)
               allSTs(layi, permi, pxi) = sum(d2pT(clustinfo.PixelIdxList{pxi}));% it is maxi -1 cause we are adding a zero below
            end
        else
            allSTs(layi, permi, :) = 0;
        end

        [maxh id] = max(abs(allSTs(layi, permi, :)));
        max_clust_sum_perm(permi,layi, :) = allSTs(layi, permi,id); 
    end
    
end

%%

for layi = 1:nLays
    
    t1 = 20%max_clust_sum_obs(layi); 
    [counts,centers] = hist(max_clust_sum_perm(:,layi), 14);
    
    figure()
    h = bar(centers, counts); hold on;
    h.FaceColor = 'w';
    h.EdgeColor = 'k';
    h.LineWidth =2;
    set(gca, 'FontSize', 20, 'ylim', [0 50] );
    plot ([t1 t1], get(gca, 'ylim'), 'LineWidth', 3,'Color', [.1 .3 .6] );
    filename = ['histogram.png'];
    %export_fig(1, filename,'-transparent', '-r300');
    %close all;
end


%% check for the correlations with negative mean

clear fr max_clust_sum_perm allSTs
for permi = 1:100
    for subji = 1:10
        x = rand(50, 50, 40);
        y = rand(50, 50, 40);

        for timei = 1:40
            x1 = x(:, :, timei); x1 = x1(:); 
            y1 = y(:, :, timei); y1 = y1(:); 
            fr(subji, timei, :) = corr(x1, y1, 'type', 's');
        end  
    end
    [h p ci ts] = ttest(fr);
    t = ts.tstat; 
    clustinfo = bwconncomp(h);
    if length(clustinfo.PixelIdxList) > 0
        for pxi = 1:length(clustinfo.PixelIdxList)
           allSTs(permi, pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% it is maxi -1 cause we are adding a zero below
        end
    else
        allSTs(permi, :) = 0;
    end

    [maxh id] = max(abs(allSTs(permi, :)));
    max_clust_sum_perm(permi,:) = allSTs(permi,id); 
    
end


%%

t1 = 0%max_clust_sum_obs(layi); 
[counts,centers] = hist(max_clust_sum_perm, 24);

figure()
h = bar(centers, counts); hold on;
h.FaceColor = 'w';
h.EdgeColor = 'k';
h.LineWidth =2;
set(gca, 'FontSize', 20, 'ylim', [0 50] );
plot ([t1 t1], get(gca, 'ylim'), 'LineWidth', 3,'Color', [.1 .3 .6] );
filename = ['histogram.png'];
%export_fig(1, filename,'-transparent', '-r300');
    %close all;


%% get rank (for stats)

clear p
for layi = 1:8
    clear allAb
    mcsR = max_clust_sum_obs(layi); 
    mcsP = max_clust_sum_perm(:,layi); 
    allAb = mcsP(abs(mcsP) > abs(mcsR));
    p(layi, :) = 1 - ((nPerm-1) - (length (allAb)))  / nPerm;
end

p

% % % %  5 % % % 5 %5 % % 5 % 



%% 3D clusters 

nPerm = size(tPerm, 1); 
f = 1:54;

clear max_clust_sum 
for permi = 1:nPerm
    clear allSTs   
    d2pT = squeeze(tPerm(permi, :,f,:));
    d2pH = squeeze(hPerm(permi, :,f,:));

    clustinfo = bwconncomp(d2pH);
    if length(clustinfo.PixelIdxList) > 0
        for pxi = 1:length(clustinfo.PixelIdxList)
           allSTs(permi, pxi) = sum(d2pT(clustinfo.PixelIdxList{pxi}));% it is maxi -1 cause we are adding a zero below
        end
    else
        allSTs(permi, :) = 0;
    end

    [maxh id] = max(abs(allSTs(permi, :)));
    max_clust_sum_perm3D(permi,:) = allSTs(permi,id); 
    
end


%% 

t1 = max_clust_sum_3D; 
[counts,centers] = hist(max_clust_sum_perm3D, 14);

figure()
h = bar(centers, counts); hold on;
h.FaceColor = 'w';
h.EdgeColor = 'k';
h.LineWidth =2;
set(gca, 'FontSize', 20, 'ylim', [0 300] );
plot ([t1 t1], get(gca, 'ylim'), 'LineWidth', 3,'Color', [.1 .3 .6] );
filename = ['histogram.png'];
%export_fig(1, filename,'-transparent', '-r300');
%close all;







%% STIMULUS REPRESENTATION IN DNNs
% % % 
clear
f2sav       = 'RNN';
loadNet_WM;



















%% RNN all RDMS

figure()
for layi = 1:56
   subplot (7, 8, layi)
   d2p = squeeze(act_CH(layi, :,:)); 
   imagesc(d2p); axis square; 
   set(gca,'cLim', [0 1])
   set(gca,'XTick',[], 'YTick', [], 'xticklabel',[])
   title(num2str(layi))
end

set(gcf, 'Position', [100 100 700 700]); 
ha=get(gcf,'children');
n = .1;
count = 56;
rowc = 8; 
for rowi = 7:-1:1
    for coli = 8:-1:1
        set(ha(count),'position',[0+rowi/9 0+coli/9 n n ])
        count = count-1;
        rowc = rowc-1; 
    end
end

exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);

%% RNN all RDMS

figure()
for layi = 1:56
   subplot (7, 8, layi)
   d2p = squeeze(act_CH(layi, :,:)); 
   imagesc(d2p); axis square; 
   set(gca,'cLim', [0 1])
   set(gca,'XTick',[], 'YTick', [], 'xticklabel',[])
   %title(num2str(layi)) % check order
end

set(gcf, 'Position', [100 100 700 700]); 
ha=get(gcf,'children');
n = .1;
count = 56;
for rowi = 1:7
    for coli = 1:8
        set(ha(count),'position',[0+coli/9 0+rowi/9 n n ])
        count = count+-1;
        
    end
end

exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);


%% RNN all MDS
cols = brewermap(6, 'Dark2') % % Scheme|'Accent'|'Dark2'|'Paired'|'Pastel1'|'Pastel2'|'Set1'|'Set2'|'Set3'|
cols = repelem(cols, 10, 1);

figure()
for layi = 1:56
   subplot (7, 8, layi)
   d2p = squeeze(act_CH(layi, :,:)); 
   d2p = 1- d2p;
   [rdmMDS] = cmdscale(d2p);
   scatter(rdmMDS(:,1),rdmMDS(:,2),350, cols, '.'); axis square
   set(gca,'XTick',[], 'YTick', [], 'xticklabel',[])
   box on
end

set(gcf, 'Position', [100 100 700 700]); 
ha=get(gcf,'children');
n = .1;
count = 56;
for rowi = 1:7
    for coli = 1:8
        set(ha(count),'position',[0+rowi/9 0+coli/9 n n ])
        count = count-1;
    end
end

exportgraphics(gcf, 'matrixRNN.png', 'Resolution', 300);

%% plot only last time point

figure(); set(gcf, 'Position', [100 100 500 500]);
lois  = [8 16 24 32 40 48 54]; 
for layi = 1:7
   subplot (3, 3, layi)
   d2p = squeeze(act_CH(lois(layi), :,:)); 
   imagesc(d2p); axis square; 
   set(gca,'cLim', [0 1])
   set(gca,'XTick',[], 'YTick', [], 'xticklabel',[])
end
 ha=get(gcf,'children');
 n = .25
 set(ha(1),'position',[0 0 n n ])
 set(ha(2),'position',[.55 .27 n n ])
 set(ha(3),'position',[.275 .27 n n])
 set(ha(4),'position',[.0 .27 n n ])
 set(ha(5),'position',[.55 .54 n n ])
 set(ha(6),'position',[.275 .54 n n])
 set(ha(7),'position',[.0 .54 n n ])
exportgraphics(gcf, 'allM.png', 'Resolution', 300);

%% plot only last time point one line

figure(); set(gcf, 'Position', [100 100 1500 500]);
lois  = [8 16 24 32 40 48 54]; 
for layi = 1:7
   subplot (1, 7, layi)
   d2p = squeeze(act_CH(lois(layi), :,:)); 
   imagesc(d2p); axis square; 
   set(gca,'cLim', [0 1])
   set(gca,'XTick',[], 'YTick', [], 'xticklabel',[])
end
 ha=get(gcf,'children');
 n = .25; 
 set(ha(1),'position',[ 6/9 0 n n ])
 set(ha(2),'position',[ 5/9 0 n n ])
 set(ha(3),'position',[ 4/9 0 n n])
 set(ha(4),'position',[ 3/9 0 n n ])
 set(ha(5),'position',[ 2/9 0 n n ])
 set(ha(6),'position',[ 1/9 0 n n])
 set(ha(7),'position',[ 0 0 n n ])
exportgraphics(gcf, 'allM.png', 'Resolution', 300);




%% plot only last time point one line vertical

figure(); set(gcf, 'Position', [100 100 700 700]);
lois  = [8 16 24 32 40 48 54]; 
for layi = 1:7
   subplot (1, 7, layi)
   d2p = squeeze(act_CH(lois(layi), :,:)); 
   imagesc(d2p); axis square; 
   set(gca,'cLim', [0 1])
   set(gca,'XTick',[], 'YTick', [], 'xticklabel',[])
end
 ha=get(gcf,'children');
 n = .1; 
 set(ha(1),'position',[0 6/9 n n ])
 set(ha(2),'position',[0 5/9 n n ])
 set(ha(3),'position',[0 4/9 n n])
 set(ha(4),'position',[0 3/9 n n ])
 set(ha(5),'position',[0 2/9 n n ])
 set(ha(6),'position',[0 1/9 n n])
 set(ha(7),'position',[0 .0 n n ])
exportgraphics(gcf, 'allM.png', 'Resolution', 300);


%% plot only last time point one line ALEXNET

figure(); set(gcf, 'Position', [100 100 1500 500]);
for layi = 1:8
   subplot (1, 8, layi)
   d2p = squeeze(act_CH(layi, :,:)); 
   imagesc(d2p); axis square; 
   set(gca,'cLim', [0 1])
   set(gca,'XTick',[], 'YTick', [], 'xticklabel',[])
end
 ha=get(gcf,'children');
 n = .25; 
 set(ha(1),'position',[.09 *7  0 n n ])
 set(ha(2),'position',[.09 *6  0 n n ])
 set(ha(3),'position',[.09 * 5 0 n n ])
 set(ha(4),'position',[.09 * 4 0 n n])
 set(ha(5),'position',[.09 * 3 0 n n ])
 set(ha(6),'position',[.09 * 2 0 n n ])
 set(ha(7),'position',[.09 0 n n])
 set(ha(8),'position',[.0 0 n n ])
exportgraphics(gcf, 'allM.png', 'Resolution', 300);



%% plot MDS only last time point one line
cols = brewermap(6, 'Dark2') % % Scheme|'Accent'|'Dark2'|'Paired'|'Pastel1'|'Pastel2'|'Set1'|'Set2'|'Set3'|
cols = repelem(cols, 10, 1);


figure(); set(gcf, 'Position', [100 100 1500 500]);
lois  = [8 16 24 32 40 48 54]; 
for layi = 1:7
   subplot (1, 7, layi)
   d2p = squeeze(act_CH(lois(layi), :,:)); 
   d2p = 1- d2p;
   [rdmMDS] = cmdscale(d2p);
   scatter(rdmMDS(:,1),rdmMDS(:,2),350, cols, '.'); axis square
   set(gca,'XTick',[], 'YTick', [], 'xticklabel',[])
   box on
end
 ha=get(gcf,'children');
 n = .25; 
 set(ha(1),'position',[.09 *6  0.1 n n ])
 set(ha(2),'position',[.09 * 5 0.1 n n ])
 set(ha(3),'position',[.09 * 4 0.1 n n])
 set(ha(4),'position',[.09 * 3 0.1 n n ])
 set(ha(5),'position',[.09 * 2 0.1 n n ])
 set(ha(6),'position',[.09 0.1 n n])
 set(ha(7),'position',[.0 0.1 n n ])
 
 
exportgraphics(gcf, 'allM.png', 'Resolution', 300);

%% plot MDS only last time point one line ALEXNET
cols = brewermap(6, 'Dark2') % % Scheme|'Accent'|'Dark2'|'Paired'|'Pastel1'|'Pastel2'|'Set1'|'Set2'|'Set3'|
cols = repelem(cols, 10, 1);


figure(); set(gcf, 'Position', [100 100 1500 500]);
for layi = 1:8
   subplot (1, 8, layi)
   d2p = squeeze(act_CH(layi, :,:)); 
   d2p = 1- d2p;
   [rdmMDS] = cmdscale(d2p);
   scatter(rdmMDS(:,1),rdmMDS(:,2),350, cols, '.'); axis square
   set(gca,'XTick',[], 'YTick', [], 'xticklabel',[])
   box on
end
 ha=get(gcf,'children');
 n = .25; 
 set(ha(1),'position',[.09 *7  0 n n ])
 set(ha(2),'position',[.09 *6  0 n n ])
 set(ha(3),'position',[.09 * 5 0 n n ])
 set(ha(4),'position',[.09 * 4 0 n n])
 set(ha(5),'position',[.09 * 3 0 n n ])
 set(ha(6),'position',[.09 * 2 0 n n ])
 set(ha(7),'position',[.09 0 n n])
 set(ha(8),'position',[.0 0 n n ])
 
 
exportgraphics(gcf, 'allM.png', 'Resolution', 300);

%% 
figure
plot (1:10,1:10, 'Color', [0.901960784313726,0.670588235294118,0.00784313725490196], 'linewidth', 10)

%% MDS all layers (black and yellow circles plot)

act = act_CH; 
act1 = arrayfun(@(i)tril(squeeze(act(i,:,:)), -1), 1:size(act,1), 'un', 0);
act2 = cat(3, act1{:}); act2(act2==0) = nan; act2 = permute(act2, [3, 1, 2]);
act2 = reshape(act2, 56, []); act3 = act2(:,all(~isnan(act2)));  


allM = corr(act3', 'type', 'p');allMS = allM.^2;

% % % matrix
%imagesc(allMS); axis square; colorbar

% % % mds
%c1 = [ones(1, 8) ones(1, 8)*2 ones(1, 8)*3 ones(1, 8)*4 ones(1, 8)*5 ones(1, 8)*6 ones(1, 8)*7]' /7; % sorted by layer
%c2 = [ones(1, 8) ones(1, 8)*2 ones(1, 8)*3 ones(1, 8)*4 ones(1, 8)*5 ones(1, 8)*6 ones(1, 8)*7]' /7; % sorted by layer
c1 = repmat((0:7)/7, 1, 7)'; % sorted by time point
c2 = repmat((0:7)/7, 1, 7)'; % sorted by time point
c3 = repmat(zeros(1), 56, 1);

cols = [c1 c2 c3];



d2p = 1- allM;
[rdmMDS] = cmdscale(d2p);
figure()
plot(rdmMDS(1:8,1),rdmMDS(1:8,2),'k', 'linewidth', 2);hold on; 
plot(rdmMDS(9:16,1),rdmMDS(9:16,2),'k', 'linewidth', 2)
plot(rdmMDS(17:24,1),rdmMDS(17:24,2),'k', 'linewidth', 2)
plot(rdmMDS(25:32,1),rdmMDS(25:32,2),'k', 'linewidth', 2)
plot(rdmMDS(33:40,1),rdmMDS(33:40,2),'k', 'linewidth', 2)
plot(rdmMDS(41:48,1),rdmMDS(41:48,2),'k', 'linewidth', 2)
plot(rdmMDS(49:56,1),rdmMDS(49:56,2),'k', 'linewidth', 2)
scatter(rdmMDS(:,1),rdmMDS(:,2),3500,cols, '.'); 
set(gca,'FontSize', 26);



%scatter3(rdmMDS(:,1),rdmMDS(:,2),rdmMDS(:,3),500, '.'); hold on; 
%plot(rdmMDS(:,1),rdmMDS(:,2)); hold on; 
%brewermap_view({gca})

%exportgraphics(gcf, 'allM.png', 'Resolution', 300);

%% MDS all layers (black and yellow circles plot) only for first and last time points

act = act_CH; 
act1 = arrayfun(@(i)tril(squeeze(act(i,:,:)), -1), 1:size(act,1), 'un', 0);
act2 = cat(3, act1{:}); act2(act2==0) = nan; act2 = permute(act2, [3, 1, 2]);
act2 = reshape(act2, 56, []); act3 = act2(:,all(~isnan(act2)));  


allM = corr(act3', 'type', 'p');%allMS = allM.^2;

% % % matrix
%imagesc(allMS); axis square; colorbar

% % % mds
%c1 = [ones(1, 8) ones(1, 8)*2 ones(1, 8)*3 ones(1, 8)*4 ones(1, 8)*5 ones(1, 8)*6 ones(1, 8)*7]' /7; % sorted by layer
%c2 = [ones(1, 8) ones(1, 8)*2 ones(1, 8)*3 ones(1, 8)*4 ones(1, 8)*5 ones(1, 8)*6 ones(1, 8)*7]' /7; % sorted by layer
clear c1 c2 c3 cols
c1 = (1:7)/7'; % sorted by time point
%c2 = repmat(zeros(1), 7, 1)';

c2 = repmat(zeros(1), 7, 1)';

c3 = repmat(ones(1), 7, 1)';

cols = [c1 ; c2 ; c3]';
cols = repelem(cols, 2,1);



d2p = 1- allM;
[rdmMDS] = cmdscale(d2p);
figure()
lois         = [1 8 9 16 17 24 25 32 33 40 41 48 49 54]; 
plot(rdmMDS([1 8],1),rdmMDS([1 8],2),'k', 'linewidth', 2);hold on; 
plot(rdmMDS([9 16],1),rdmMDS([9 16],2),'k', 'linewidth', 2);hold on; 
plot(rdmMDS([17 24],1),rdmMDS([17 24],2),'k', 'linewidth', 2);hold on; 
plot(rdmMDS([25 32],1),rdmMDS([25 32],2),'k', 'linewidth', 2);hold on; 
plot(rdmMDS([33 40],1),rdmMDS([33 40],2),'k', 'linewidth', 2);hold on; 
plot(rdmMDS([41 48],1),rdmMDS([41 48],2),'k', 'linewidth', 2);hold on; 
plot(rdmMDS([49 54],1),rdmMDS([49 54],2),'k', 'linewidth', 2);hold on; 



fmt = {'o'; '^'};
fmt = repelem(fmt, 1, 7);

for i = 1:14
    scatter(rdmMDS(lois(i),1),rdmMDS(lois(i),2),200,cols(i, :),fmt{i}, 'Markerfacecolor', cols(i,:)); 
    set(gca,'FontSize', 20, 'xlim', [-.4 .4], 'ylim', [-.2 .5]); 
end



%scatter3(rdmMDS(:,1),rdmMDS(:,2),rdmMDS(:,3),500, '.'); hold on; 
%plot(rdmMDS(:,1),rdmMDS(:,2)); hold on; 
%brewermap_view({gca})

exportgraphics(gcf, 'allM.png', 'Resolution', 300);


%% Representational consistency all layers / time points

act = act_CH; 
act1 = arrayfun(@(i)tril(squeeze(act(i,:,:)), -1), 1:size(act,1), 'un', 0);
act2 = cat(3, act1{:}); act2(act2==0) = nan; act2 = permute(act2, [3, 1, 2]);
act2 = reshape(act2, 56, []); act3 = act2(:,all(~isnan(act2)));  


allM = corr(act3', 'type', 's');allMS = allM.^2;

% % % matrix
figure()
imagesc(allMS); axis square; colorbar
set(gca, 'FontSize', 15)


%% Representational consistency only last time point

lois         = [1 9 17 25 33 41 49]; 
%lois         = [8 16 24 32 40 48 54]; 


act = act_CH; 
act1 = arrayfun(@(i)tril(squeeze(act(i,:,:)), -1), 1:size(act,1), 'un', 0);
act2 = cat(3, act1{:}); act2(act2==0) = nan; act2 = permute(act2, [3, 1, 2]);
act2 = reshape(act2, 56, []); act3 = act2(:,all(~isnan(act2)));  


allM = corr(act3', 'type', 's');
allMS = allM(lois, lois);

% % % matrix
figure()
imagesc(allMS); axis square; colorbar
set(gca, 'FontSize', 20)
set(gca, 'xtick', [1:7], 'xticklabel', {[1:7]}, 'ytick', [1:7], 'yticklabel', {[1:7]},'clim', [.5 1])
exportgraphics(gcf, 'allM.png', 'Resolution', 300);



%% Representational consistency ALEXNET

act = act_CH; 
act1 = arrayfun(@(i)tril(squeeze(act(i,:,:)), -1), 1:size(act,1), 'un', 0);
act2 = cat(3, act1{:}); act2(act2==0) = nan; act2 = permute(act2, [3, 1, 2]);
act2 = reshape(act2, 8, []); act3 = act2(:,all(~isnan(act2)));  


allM = corr(act3', 'type', 's');allMS = allM.^2;

% % % matrix
figure()
imagesc(allMS); axis square; colorbar
set(gca, 'FontSize', 20)
set(gca, 'xtick', [1:8], 'xticklabel', {[1:8]}, 'ytick', [1:8], 'yticklabel', {[1:8]},'clim', [0.3 1])
exportgraphics(gcf, 'allM.png', 'Resolution', 300);



%% Changes in representational consistency 

act = act_CH; 
act1 = arrayfun(@(i)tril(squeeze(act(i,:,:)), -1), 1:size(act,1), 'un', 0);
act2 = cat(3, act1{:}); act2(act2==0) = nan; act2 = permute(act2, [3, 1, 2]);
act2 = reshape(act2, 56, []); act3 = act2(:,all(~isnan(act2)));  

c1 = (1:7)/7'; % sorted by time point
c2 = repmat(zeros(1), 7, 1)';
c3 = repmat(ones(1), 7, 1)';
cols = [c1 ; c2 ; c3]';


clear allM2
for tlyi = 1:55

    allM2(tlyi,:) = corr(act3(tlyi,:)',act3(tlyi+1,:)', 'type', 's');allMS = allM;
    
end

figure()
lw = 3;

plot(abs(allM2(1:7)), 'Linewidth', lw, 'Color', [cols(1,:)]); hold on; 
plot(abs(allM2(9:15)), 'Linewidth', lw, 'Color', [cols(2,:)])
plot(abs(allM2(17:23)), 'Linewidth', lw, 'Color', [cols(3,:)])
plot(abs(allM2(25:31)), 'Linewidth', lw, 'Color', [cols(4,:)])
plot(abs(allM2(33:39)), 'Linewidth', lw, 'Color', [cols(5,:)])
plot(abs(allM2(41:47)), 'Linewidth', lw, 'Color', [cols(6,:)])
plot(abs(allM2(49:55)), 'Linewidth', lw, 'Color', [cols(7,:)])
%legend({'Layer 1' 'Layer 2' 'Layer 3' 'Layer 4' 'Layer 5' 'Layer 6' 'Layer 7' })
set(gca, 'ylim', [0.85 1], 'xlim', [0 8], 'xtick', [1:7], 'xticklabels', {'1-2' '2-3' '3-4' '4-5' '5-6' '6-7' '7-8'}, 'FontSize', 20)
exportgraphics(gcf, 'allM.png', 'Resolution', 300);


%% example pytorch DQN

layers = importONNXLayers('test_model.onnx','OutputLayerType', 'classification', 'ImportWeights', true)
layers = removeLayers(layers,'Transpose_0')
layers = removeLayers(layers,'Input_input')
inputlayer = imageInputLayer([84 84 3],'Name','input', 'Normalization', 'none')
layers = addLayers(layers,inputlayer)
layers = connectLayers(layers,'input', 'Conv_1')
analyzeNetwork(layers)
%%
net = assembleNetwork(layers);

%% 
analyzeNetwork(net)

%%
im = rand(84, 84, 3);

act1 = activations(net,im,'Conv_1');


%% START HERE WITH THE ALL TRIALS ANALYSIS
% 1) bands
clearvars -except act_CH act_FR 
f2sav       = 'Alex_pfc_CueAll_noAv_54_Real_nT_1-49'
f           = 3:54; %in case B
load_parameters_WMFRA;
loadNet_WM;

all_r_Times    = zeros(nSubj, nLays, nFreqs, nTimes);

tic
for subji=1:nSubj
    disp (['Subj: ' num2str(subji)])
    load(sublist{subji});
    cfg_contrasts.oneListPow = cfg_contrasts.oneListPow(:, :, :, resTime); 
    cfg_contrasts = normalize_WM(cfg_contrasts, 1, 'blo', [3 7]);
    if strcmp(avMeth,'pow')
       cfg_contrasts = average_repetitions(cfg_contrasts);
    end
    
    f2t = strsplit(f2sav, '_'); 
    if strcmp(f2t{5}, '54') | strcmp(f2t{5}, '150') 
        for freqi = 1:length(freqs2test)
            f  = freqs2test(freqi);
            rdm_prev(freqi) = create_rdms(cfg_contrasts, f, it, 5, 1);
        end
    elseif strcmp(f2t{5}, 'B') 
        rdm_prev = create_rdms(cfg_contrasts, f, it, 5, 1);
        if subji == 1 % change name only once
            f2sav = [f2sav num2str(f(1)) '-' num2str(f(end)) 'Hz'];
        end
    end
    
    [allS ids act_FR2 act_CH2] = getIDsAndnRDMs(rdm_prev, act_CH, act_FR); 
   
    
    for layi = 1:nLays
        if subji < subj_ch_fr
            M =  squeeze(act_FR2(layi,:,:)); 
        else
            M =  squeeze(act_CH2(layi,:,:)); 
        end
        M(M ==1) = 1000;M(M ==0) = 2000;
        M = tril(M, -1);M(M==0) =[];
        M(M==1000) = 1;M(M==2000) = 0;

        allTEst = corr(allS, M', 'type', 's');

        all_r_Times(subji, layi,:, : ) = reshape(allTEst, nFreqs, nTimes );  
    end
 
end 

cd (currentF) 
toc 
save(f2sav, 'all_r_Times');


%% PLOT OBS DATA (all layers) - case with mix negative and positive clusters

sub2exc =  [];

nSubjs =size(all_r_Times, 1);  
nLays = size(all_r_Times, 2); 
nFreqs = size(all_r_Times, 3); 
clear hL tL clustinfo
for freqi = 1:nFreqs
    allR = squeeze(all_r_Times(:, :, freqi, :));
    allR(sub2exc,:,:) = []; 
    for layi = 1:nLays
        [hL(layi, freqi, :) p ci ts] = ttest(allR(:,layi,:)); 
        tL(layi, freqi, :) = ts.tstat;
    end
end

f = 1:52; % from 3 to 54

clear max_clust_sum_obs allSTs
for layi = 1:nLays
    
    hLay = squeeze(hL(layi,f,:)); 
    tLay = squeeze(tL(layi,f,:));
    
    clear allSTs  
    clustinfo = bwconncomp(hLay);
    for pxi = 1:length(clustinfo.PixelIdxList)
        % check whether it is a combined + and - cluster
        V = tLay(clustinfo.PixelIdxList{pxi});
        Vids = clustinfo.PixelIdxList{pxi}; 
        if ~any(diff(sign(V(V~=0)))) %all time bins have tvalues of same sie
            allSTs(pxi,:) = sum(tLay(clustinfo.PixelIdxList{pxi}));
        else %remove the 
            big0 = V>0; small0 = V<0; 
            VidsS = Vids(small0); 
            ids2k = Vids(big0); 
            if sum(big0) > sum(small0)
                V(small0) = NaN; 
                hLay(VidsS) = 0; 
                clustinfo.PixelIdxList{pxi} = ids2k; 
            else
                %V(big0) = NaN; 
                %hLay(big0,:) = NaN; 
            end
            
            allSTs(pxi,:) = sum(V, 'omitnan'); 
        end
    end

    [maxh id] = max(abs(allSTs));
    max_clust_sum_obs(layi,:) = allSTs(id); 

    % % % % rem non-sig-clusters
    for ci = 1:length(clustinfo.PixelIdxList)
        modifiedCluster = abs(sum(tLay(clustinfo.PixelIdxList{ci})))+0.0001; %add this small number because removing negative values creates slighlly different valeus
       if modifiedCluster < max(abs(allSTs))   %+ 1    %add +1 at the very end to delete all clusters
          %disp(['layi > ' num2str(layi) ' > ' num2str(abs(sum(tLay(clustinfo.PixelIdxList{ci})))) '_' num2str(max(abs(allSTs)))]) 
          hLay (clustinfo.PixelIdxList{ci}) =  0; 
       end
    end
    
    hL(layi, f, :) = hLay; 
    tL(layi, f, :) = tLay; 
end



figure()
layT = tiledlayout(8, 7);
layT.TileSpacing = 'compact';
layT.Padding = 'compact';
for layi = 56:56 %1:nLays
    
    t = squeeze(tL(layi,:,:));
    h = squeeze(hL(layi,:,:));
    
    if size(all_r_Times, 4) < 30 & size(all_r_Times, 4) > 5
        % % % % encoding 
        if nLays > 9
            subplot(7, 8, layi)
        else
            subplot(3, 3, layi)
        end
        %times = -.5:.01:1.499;
        times = 0:.01:1.499;
        %freqs = [.1:.1:29.9 30:.5:150]; 
        freqs = 1:520;
        contourf(times, freqs, myresizem(t, 10), 100, 'linecolor', 'none'); hold on; %colorbar
        contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 2);
        plot([-0.055 -0.055], get(gca, 'ylim'), 'k:', 'LineWidth', 2);
        plot(get(gca, 'xlim'),[300 300], 'k:', 'LineWidth', 2);
        set(gca, 'xtick', [-0.055 .445 .945 1.45], 'xticklabel', [0 .5 1 1.5], 'clim', [-5 5]); colormap jet
        set(gca,'XTick',[], 'YTick', [], 'xticklabel',[], 'FontSize', 12)

    elseif size(all_r_Times, 4) == 5
        % % % % encoding 
        if nLays > 9
            subplot(7, 8, layi)
        else
            subplot(3, 3, layi)
        end
        times = 0:.01:0.499;
        %freqs = [.1:.1:29.9 30:.5:150]; 
        freqs = 1:520;
        contourf(times, freqs, myresizem(t, 10), 100, 'linecolor', 'none'); hold on; %colorbar
        contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 2);
        plot([-0.055 -0.055], get(gca, 'ylim'), 'k:', 'LineWidth', 2);
        plot(get(gca, 'xlim'),[300 300], 'k:', 'LineWidth', 2);
        set(gca, 'xtick', [-0.055 .445 .945 1.45], 'xticklabel', [0 .5 1 1.5], 'clim', [-5 5]); colormap jet
        set(gca,'XTick',[], 'YTick', [], 'xticklabel',[], 'FontSize', 12)
    else
        if nLays > 9
            subplot(7, 8, layi)
        else
            subplot(3, 3, layi)
        end
        %times = 0:.01:3.49;
        times = 0:.01:4.49;
        %freqs = [.1:.1:29.9 30:.5:150]; 
        freqs = 1:520;
        contourf(times, freqs, myresizem(t, 10), 100, 'linecolor', 'none'); hold on; 
        contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 1); %colorbar
        plot([-0.055 -0.055], get(gca, 'ylim'), 'k:', 'LineWidth', 2);
        plot(get(gca, 'xlim'),[300 300], 'k:', 'LineWidth', 2);
        set(gca, 'xtick', [-0.055 .945 1.945 2.945 3.945 ], 'xticklabel', [0 1 2 3 4], 'xlim', [0 3.5], 'clim', [-3 3])
        set(gca,'XTick',[], 'YTick', [])
        set(gca,'xticklabel',[], 'FontSize', 12)
        colormap jet
    end

end
        

f2p = 'myFig.png';
exportgraphics(gcf, f2p, 'Resolution', 300)





%% plot NICELY

h = zeros(52, 45);

figure(); set(gcf, 'Position', [100 100 1500 1000]);
% Scheme|'BrBG'|'PRGn'|'PiYG'|'PuOr'|'RdBu'|'RdGy'|'RdYlBu'|'RdYlGn'|'Spectral'|
myCmap = colormap(brewermap([],'*Spectral'));
times = 0:.01:4.49;
%freqs = [.1:.1:29.9 30:.5:150]; 
freqs = 1:520;
contourf(times, freqs, myresizem(t, 10), 100, 'linecolor', 'none'); hold on; 
contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 10); %colorbar
plot([.45 .45], get(gca, 'ylim'), 'k:', 'LineWidth', 12);
plot([.95 .95], get(gca, 'ylim'), 'k:', 'LineWidth', 12);
%plot([1.25 1.25], get(gca, 'ylim'), 'k:', 'LineWidth', 12);
plot(get(gca, 'xlim'),[300 300], 'k:', 'LineWidth', 12);
set(gca, 'xtick', [-0.055 .945 1.945 2.945 3.945 ], 'xticklabel', [0 1 2 3 4], 'xlim', [0 3.5], 'clim', [-4 4]); colorbar
set(gca,'XTick',[], 'YTick', [])
set(gca,'xticklabel',[], 'FontSize', 12)
    
exportgraphics(gcf, 'allM.png', 'Resolution', 300);



%% Frequency resolved RNN analysis for all trials and 8th time point only 

clearvars -except act_CH act_FR 

%lois         = [1 9 17 25 33 41 49]; 
lois         = [8 16 24 32 40 48 56]; 
f2sav       = 'RNN_pfc_Maint_Av_54_Real_nT_loi_8th'
load_parameters_WMFRA;
loadNet_WM;

all_r_Times    = zeros(nSubj, length(lois), nFreqs, nTimes);

tic

for subji=1:nSubj
    disp (['Subj: ' num2str(subji)])
    load(sublist{subji});
    cfg_contrasts.oneListPow = cfg_contrasts.oneListPow(:, :, :, resTime); 
    cfg_contrasts = normalize_WM(cfg_contrasts, 1, 'blo', [3 7]);
    if strcmp(avMeth,'pow')
       cfg_contrasts = average_repetitions(cfg_contrasts);
    end
    
    for freqi = 1:length(freqs2test)
        f  = freqs2test(freqi);
        rdm_prev(freqi) = create_rdms(cfg_contrasts, f, it, 5, 1);
    end
    
    [allS ids act_FR2 act_CH2] = getIDsAndnRDMs(rdm_prev, act_CH, act_FR); 
   
    
    parfor layi = 1:length(lois)
        if subji < subj_ch_fr
            M =  squeeze(act_FR2(lois(layi),:,:)); 
        else
            M =  squeeze(act_CH2(lois(layi),:,:)); 
        end
        M(M ==1) = 1000;M(M ==0) = 2000;
        M = tril(M, -1);M(M==0) =[];
        M(M==1000) = 1;M(M==2000) = 0;

        allTEst = corr(allS, M', 'type', 's');

        all_r_Times(subji, layi,:, : ) = reshape(allTEst, nFreqs, nTimes );  
    end
 
 
end 

cd (currentF) 
toc 
save(f2sav, 'all_r_Times');

%% Plot in 3 different time periods vertical separately for first and last time points in each layer

%allTC= squeeze(mean(all_r_Times(:, [1 9 17 25 33 41 49], :, :), 4, 'omitnan')); % first layer
%allTC= squeeze(mean(all_r_Times(:, [8 16 24 32 40 48 56], :, :), 4, 'omitnan')); % lastlayer
allTC= squeeze(mean(all_r_Times, 4, 'omitnan'));
[h p ci ts] = ttest(allTC);
t = squeeze(ts.tstat); h = squeeze(h);

mAllTC = squeeze(mean(allTC, 'omitnan'))';



figure(); set(gcf, 'Position', [100 100 100 500])
subplot(211)
imagesc(flipud(myresizem(mAllTC, 10))); colorbar; hold on; 
%set(gca, 'xtick', [], 'xticklabels', [], 'ytick', [], 'yticklabels', [])
plot ([0 560],  [540-300 540-300],'k:', 'LineWidth', 2);
set(gca, 'xlim', [0 70])

subplot(212)
imagesc(flipud(myresizem(t', 10))); colorbar; hold on;
contour(flipud(myresizem(h', 10)), 1, 'Color', [0, 0, 0], 'LineWidth', 1);
plot ([0 560],  [540-300 540-300],'k:', 'LineWidth', 2);
set(gca, 'xlim', [0 70])




exportgraphics(gcf, 'vertical_plots.png', 'Resolution', 300); 








%% Plot in 3 different time periods vertical

nLays = 56; 

allTC= squeeze(mean(all_r_Times, 4, 'omitnan'));
[h p ci ts] = ttest(allTC);
t = squeeze(ts.tstat); h = squeeze(h);

mAllTC = squeeze(mean(allTC, 'omitnan'))';


figure()
subplot(211)
imagesc(flipud(myresizem(mAllTC, 10))); colorbar; hold on; 
plot ([80 80], get(gca, 'xlim'), 'LineWidth', 1,'Color', [0 0 0 ] );
plot ([160 160],get(gca, 'xlim'),  'LineWidth', 1,'Color', [0 0 0 ] );
plot ([240 240], get(gca, 'xlim'), 'LineWidth', 1,'Color', [0 0 0 ] );
plot ([320 320], get(gca, 'xlim'), 'LineWidth', 1,'Color', [0 0 0 ] );
plot ([400 400], get(gca, 'xlim'), 'LineWidth', 1,'Color', [0 0 0 ] );
plot ([480 480], get(gca, 'xlim'), 'LineWidth', 1,'Color', [0 0 0 ] );
set(gca, 'xtick', [], 'xticklabels', [], 'ytick', [], 'yticklabels', [])
plot ([0 560],  [540-300 540-300],'k:', 'LineWidth', 2);axis equal
set(gca, 'xlim', [0 560])

subplot(212)
imagesc(flipud(myresizem(t', 10))); colorbar; hold on;
contour(flipud(myresizem(h', 10)), 1, 'Color', [0, 0, 0], 'LineWidth', 1);
plot ([80 80], get(gca, 'xlim'), 'LineWidth', 1,'Color', [0 0 0 ] );
plot ([160 160],get(gca, 'xlim'),  'LineWidth', 1,'Color', [0 0 0 ] );
plot ([240 240], get(gca, 'xlim'), 'LineWidth', 1,'Color', [0 0 0 ] );
plot ([320 320], get(gca, 'xlim'), 'LineWidth', 1,'Color', [0 0 0 ] );
plot ([400 400], get(gca, 'xlim'), 'LineWidth', 1,'Color', [0 0 0 ] );
plot ([480 480], get(gca, 'xlim'), 'LineWidth', 1,'Color', [0 0 0 ] );
set(gca, 'xtick', [], 'xticklabels', [], 'ytick', [], 'yticklabels', [])
plot ([0 560],  [540-300 540-300],'k:', 'LineWidth', 2);axis equal




















 
