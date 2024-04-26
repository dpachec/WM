%%
%%  Load and plot results 
%Network_ROI_ER_layers_freqs_avRepet_avTFV_fRes(0-1)_fitMode(0:noTrials; 1:Trials)_timeRes_win_mf
clear , clc

%f2sav = 'Alex_vvs_E123_[1-8]_3-54_1_0_1_0_.1_5_1'
%f2sav = 'BLNETi_pfc_E123_[8-8-56]_3-54_1_0_1_0_.1_5_1'
f2sav = 'CORrt_pfc_E123_[2-2-8]_3-54_1_0_1_0_.1_5_1'

%f2sav = 'BLNETi_vvs_M123_[8-8-56]_3-54_1_0_1_0_.1_5_1_PCC'
%f2sav = 'BLNETi_vvs_M123_[8-8-56]_3-54_0_0_1_0_.1_5_1_MASK'

%f2sav = 'ITM_vvs_M123_[1]_3-54_0_0_1_0_.1_5_1'
%f2sav = 'BLNETi_vvs_M123_[8-8-56]_3-54_1_0_1_0_.1_5_1'
%f2sav = 'Alex_pfc_M123_[1-8]_3-54_0_0_1_0_.1_5_1_MASK'

%f2sav = 'Alex_vvs_M13_[1-8]_3-54_1_0_1_0_.1_5_1'

%f2sav = 'BLNETi_pfc_M123_[8-8-56]_3-54_1_0_1_0_.1_5_1_LATERAL';
%f2sav = 'CAT_pfc_M123_[1]_3-54_1_0_1_0_.1_5_1_LATERAL';
%f2sav = 'Alex_pfc_M123_[1-8]_3-54_1_0_1_0_.1_5_1_LATERAL';
%f2sav = 'BLNETi_pfc_M123_[8-8-56]_3-54_1_0_1_0_.1_5_1_LATERAL';
%f2sav = 'Alex_pfc_E123_[1-8]_3-54_1_0_1_0_.1_5_1_LATERAL';
%f2sav = 'BLNETi_pfc_M11_[8-8-56]_3-54_1_0_1_0_.1_5_1'

%f2sav = 'BFNETi_vvs_M123_[1-7]_3-54_1_0_1_0_.1_5_1';
%f2sav = 'BKNETi_vvs_M123_[1-7]_3-54_1_0_1_0_.1_5_1';
%f2sav = 'BDNETi_vvs_M123_[1-13]_3-54_1_0_1_0_.1_5_1';

%f2sav = 'BLNETi_pfc_M123_[8-8-56]_3-54_1_0_1_0_.1_5_1_LATERAL';
%f2sav =  'BLNETi_pfc_M123_[8-8-56]_3-54_1_0_1_0_.1_5_1';

%f2sav = 'CAT_vvs_E123_[1]_3-54_1_0_1_0_.1_5_1';


cLim = [-5 5]; 

cfg = getParams(f2sav);
if strcmp(cfg.brainROI, 'vvs')
    sub2exc = [18 22];
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
end

for layi = 1:size(nnFit{2}, 1) % nnFit{1} is empty in PFC
    ax1 = nexttile;
    clear nnH
    for subji = 1:length(nnFit)
         if ~isempty(nnFit{subji, 1})
           if strcmp(cfg.period(1), 'M')
             nnH(subji, : ,:) = atanh(nnFit{subji, 1}(layi,:,1:40));
           else
             nnH(subji, : ,:) = atanh(nnFit{subji, 1}(layi,:,1:15));
           end
         end
    end
    
    nnH(sub2exc, :, :) = []; 
    nnH = squeeze(nnH);
    %[h p ci ts] = ttest(nnH, 0, "Tail","right");
    [h p ci ts] = ttest(nnH);
    h = squeeze(h); t = squeeze(ts.tstat); 
    %h(:, 1:5) = 0; % only sum p-values in clusters after the baseline
    % if strcmp(cfg.period, 'M11') | strcmp(cfg.period, 'M12') | strcmp(cfg.period, 'M13')
    %     h(:, 1:13) = 0; % only sum p-values in clusters after the baseline
    % end
    
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
        set(gca, 'xlim', [1 40], 'clim', cLim, 'FontSize', 10);
        plot([5 5],get(gca,'ylim'), 'k:','lineWidth', 2);
    else
        set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
        set(gca, 'FontSize', 8, 'clim', cLim);
        plot([5 5],get(gca,'ylim'), 'k:','lineWidth', 2);
    end

    if strcmp(cfg.period, 'M11') | strcmp(cfg.period, 'M12') | strcmp(cfg.period, 'M13')
        set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
        set(gca, 'xlim', [1 40], 'clim', cLim, 'FontSize', 10);
        plot([5 5],get(gca,'ylim'), 'k:','lineWidth', 2);
        plot([5+8 5+8],get(gca,'ylim'), 'k:','lineWidth', 2);
        % set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
        % set(gca, 'xlim', [6 25], 'clim', [-5 5], 'FontSize', 10);
        
        
    end
    

end

%exportgraphics(gcf, [paths.results.DNNs 'myP.png'], 'Resolution', 300); 
exportgraphics(gcf, ['myP.png'], 'Resolution', 300); 



%% plot nicely

tiledlayout(8,8, 'TileSpacing', 'compact', 'Padding', 'compact');
if strcmp(cfg.period(1), 'M')
    set(gcf, 'Position', [100 100 1800 1300])
else
    set(gcf, 'Position', [100 100 670 1300])
end

isTrend = 0; 
for layi = 1:size(nnFit{2}, 1) % nnFit{1} is empty in PFC
    ax1 = nexttile;
    clear nnH
    for subji = 1:length(nnFit)
         if ~isempty(nnFit{subji, 1})
           if strcmp(cfg.period(1), 'M')
             nnH(subji, : ,:) = atanh(nnFit{subji, 1}(layi,:,1:40));
           else
             nnH(subji, : ,:) = atanh(nnFit{subji, 1}(layi,:,1:15));
           end
         end
    end
    
    nnH(sub2exc, :, :) = []; 
    nnH = squeeze(nnH);
    [h p ci ts] = ttest(nnH);
    h = squeeze(h); t = squeeze(ts.tstat); 
    %h(:, 1:5) = 0; % only sum p-values in clusters after the baseline
    if strcmp(cfg.period, 'M11') | strcmp(cfg.period, 'M12') | strcmp(cfg.period, 'M13')
        h(:, 1:13) = 0; % only sum p-values in clusters after the baseline
    end
    
    d2p = squeeze(mean(nnH, 'omitnan'));
    freqs = 1:520; 
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


    freqs = 1:520; 
    if strcmp(cfg.period(1), 'M')
        times = 1:size(t, 2)*10; 
        h = zeros(52, 40); 
        if strcmp(cfg.net2load, 'BLNETi') & layi == 7
            h(allClustInfo{7}.PixelIdxList{2}) = 1; 
            isTrend = 0; 
        end
        if strcmp(cfg.net2load, 'BLNETe') & layi == 7
            h(allClustInfo{7}.PixelIdxList{2}) = 1; 
            isTrend = 0; 
        end
        if strcmp(cfg.net2load, 'BKNETi') & layi == 4
            h(allClustInfo{4}.PixelIdxList{16}) = 1; 

        end
        if strcmp(cfg.net2load, 'BLNETi') & strcmp(cfg.meth, 'LATERAL') & layi == 7
            h(allClustInfo{7}.PixelIdxList{1}) = 1; 
        end


    else
        times = 1:150; 
        h = zeros(52, 15); 
        if strcmp(cfg.brainROI, 'vvs')
            if strcmp(cfg.net2load, 'Alex') & (layi == 1)
                h(allClustInfo{1}.PixelIdxList{1}) = 1; 
            end
            if strcmp(cfg.net2load, 'Alex') & (layi == 2)
                h(allClustInfo{2}.PixelIdxList{4}) = 1; 
            end
            if strcmp(cfg.net2load, 'Alex') & (layi == 3)
                h(allClustInfo{3}.PixelIdxList{5}) = 1; 
            end
            if strcmp(cfg.net2load, 'Alex') & (layi == 4)
                h(allClustInfo{4}.PixelIdxList{4}) = 1; 
                h(allClustInfo{4}.PixelIdxList{3}) = 1; 
            end
            if strcmp(cfg.net2load, 'Alex') & (layi == 5)
                h(allClustInfo{5}.PixelIdxList{3}) = 1; 
            end
            if strcmp(cfg.net2load, 'Alex') & (layi == 6)
                h(allClustInfo{layi}.PixelIdxList{4}) = 1; 
                h(allClustInfo{layi}.PixelIdxList{3}) = 1; 
            end
            if strcmp(cfg.net2load, 'Alex') & (layi ==7  | layi ==8)
                h(allClustInfo{layi}.PixelIdxList{3}) = 1; 
            end


            if strcmp(cfg.net2load, 'BLNETi') & (layi ==1  | layi ==2 )
                h(allClustInfo{layi}.PixelIdxList{2}) = 1; 
            end
            if strcmp(cfg.net2load, 'BLNETi') & (layi ==3  | layi ==4  | layi ==5  )
                h(allClustInfo{layi}.PixelIdxList{3}) = 1; 
            end
            if strcmp(cfg.net2load, 'BLNETi') & (layi ==6  | layi ==7)
                h(allClustInfo{layi}.PixelIdxList{4}) = 1; 
            end
            if strcmp(cfg.net2load, 'CORrt') & (layi ==1 )
                h(allClustInfo{layi}.PixelIdxList{1}) = 1; 
            end
            if strcmp(cfg.net2load, 'CORrt') & (layi ==2 )
                h(allClustInfo{layi}.PixelIdxList{3}) = 1; 
            end
            if strcmp(cfg.net2load, 'CORrt') & (layi ==3)
                h(allClustInfo{layi}.PixelIdxList{5}) = 1; 
            end
            if strcmp(cfg.net2load, 'CORrt') & (layi ==4)
                h(allClustInfo{layi}.PixelIdxList{2}) = 1; 
            end
   
            


        else


        end

        if strcmp(cfg.net2load, 'BLNETi') & strcmp(cfg.meth, 'LATERAL') & layi == 2
            h(allClustInfo{2}.PixelIdxList{2}) = 1; 
        end
        if strcmp(cfg.net2load, 'BLNETi') & strcmp(cfg.meth, 'LATERAL') & layi == 3
            h(allClustInfo{3}.PixelIdxList{2}) = 1; 
        end
        if strcmp(cfg.net2load, 'BLNETi') & strcmp(cfg.meth, 'LATERAL') & layi == 4
            h(allClustInfo{4}.PixelIdxList{2}) = 1; 
        end
        if strcmp(cfg.net2load, 'CAT') & layi == 1
            h(allClustInfo{1}.PixelIdxList{1}) = 1; 
            h(5, [1:3]) = 0; 
        end
    end
    myCmap = colormap(brewermap([],'*spectral'));
    colormap(myCmap)
    contourf(times, freqs,myresizem(t, 10), 100, 'linecolor', 'none'); hold on; %colorbar
    if ~isTrend
        contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 2);
    else
        contour(times, freqs, myresizem(h, 10), 1, ':', 'Color', [0, 0, 0], 'LineWidth', 2);
    end
    
    if strcmp(cfg.period(1), 'M') & ~strcmp(cfg.period, 'M11') & ~strcmp(cfg.period, 'M12') & ~strcmp(cfg.period, 'M13')
        set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
        set(gca, 'xlim', [1 400], 'clim', [-4 4], 'FontSize', 10);
        plot([45 45],get(gca,'ylim'), 'k:','lineWidth', 4);
    elseif strcmp(cfg.period(1), 'E')
        set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
        set(gca, 'FontSize', 8, 'clim', [-5 5]);
        plot([45 45],get(gca,'ylim'), 'k:','lineWidth', 4);
    end
    if strcmp(cfg.period, 'M11') | strcmp(cfg.period, 'M12') | strcmp(cfg.period, 'M13')
        set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', [], 'clim', [-5 5]); 
        plot([45 45],get(gca,'ylim'), 'k:','lineWidth', 3);
        plot([45+80 45+80],get(gca,'ylim'), 'k:','lineWidth', 3);
    end
    
    

end

%exportgraphics(gcf, [paths.results.DNNs 'myP.png'], 'Resolution', 300); 
exportgraphics(gcf, ['myP.png'], 'Resolution', 300); 
close all

%% plot nicely BLNETe
sub2exc = 1; 
tiledlayout(8,8, 'TileSpacing', 'compact', 'Padding', 'compact');
if strcmp(cfg.period(1), 'M')
    set(gcf, 'Position', [100 100 1800 1300])
else
    set(gcf, 'Position', [100 100 670 1300])
end

isTrend = 0; 
for layi = 7
    
    freqs = 1:520; 
    clustinfo = bwconncomp(h);


    freqs = 1:520; 
    
    times = 1:size(t, 2)*10; 
    h = zeros(52, 40); 
    h(allClustInfo{7}.PixelIdxList{2}) = 1; 
    isTrend = 0; 

    
    

    myCmap = colormap(brewermap([],'*spectral'));
    colormap(myCmap)
    contourf(times, freqs,myresizem(t, 10), 100, 'linecolor', 'none'); hold on; %colorbar
    if ~isTrend
        contour(times, freqs, myresizem(h, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 8);
    else
        contour(times, freqs, myresizem(h, 10), 1, ':', 'Color', [0, 0, 0], 'LineWidth', 8);
    end
    
    set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
    set(gca, 'xlim', [1 400], 'clim', [-5 5], 'FontSize', 20);
    plot([45 45],get(gca,'ylim'), 'k:','lineWidth', 15);
    

end

%exportgraphics(gcf, [paths.results.DNNs 'myP.png'], 'Resolution', 300); 
exportgraphics(gcf, ['myP.png'], 'Resolution', 300); 
close all



%% COMPUTE CLUSTERS in each permutation FREQUENCY RESOLVED - loads the same file as in the previous block with p1000
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

p (p==1.0010 | p == 1 | p == 1.0100) = nan; 
p_ranked = p; p_ranked(isnan(p_ranked)) = []; 
p_ranked = sort(p_ranked(:));
p

%% p last layer only
%p = p (end,:);
p_ranked = p; p_ranked(p_ranked == 0 | isnan(p_ranked)) = []; 
p_ranked = sort(p_ranked(:))




%%  plot SEPARATELY FOR CORRECT AND INCORRECT starting from the trial level fits
%Network_ROI_ER_layers_freqs_avRepet_avTFV_fRes(0-1)_fitMode(0:noTrials; 1:Trials)_timeRes_win_mf
clear , clc

%f2sav = 'Alex_pfc_M123_[1-8]_3-54_0_0_1_1_.1_5_1';
%f2sav = 'BLNETi_vvs_M123_[48]_3-54_0_0_1_1_.1_5_1';
%f2sav = 'BLNETi_vvs_M123_[8-8-56]_3-54_0_0_1_1_.1_5_1';
%f2sav = 'BLNETi_vvs_M123_[8-8-56]_3-54_0_0_1_1_.1_5_1';
f2sav = 'BLNETi_vvs_M123_[8-8-56]_3-54_0_0_1_1_.1_5_1';
cfg = getParams(f2sav);
paths = load_paths_WM(cfg.brainROI, cfg.net2load);
load([paths.results.DNNs f2sav '.mat']);


%%
cond2plot = 'CC'; %allT CC CI
minSubCrit = 5; 


if strcmp(cfg.brainROI, 'vvs')
    sub2exc = [18 22];
elseif strcmp(cfg.brainROI, 'pfc')
    sub2exc = [1];
end


tiledlayout(8,9, 'TileSpacing', 'tight', 'Padding', 'compact');
if strcmp(cfg.period(1), 'M')
    set(gcf, 'Position', [100 100 1200 1000])
else
    set(gcf, 'Position', [100 100 700 1200])
    %set(gcf, 'Position', [100 100 1800 1000])
end

for layi = 1:size(nnFit{2}, 1)
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
                ids = ids(:,20)==1; 
            elseif strcmp(cond2plot, 'IC')
                ids = ids(:,20)==0; 
            elseif strcmp(cond2plot, 'CI')
                ids = ids(:,19)==1; 
            elseif strcmp(cond2plot, 'II')
                ids = ids(:,19)==0; 
            elseif strcmp(cond2plot, 'allT')
                ids = logical(ones(1, size(ids, 1))); %just takes them all
            end

             nTR(subji,:) = sum(ids==1); 
             nnH(subji, : ,:) = squeeze(mean(avTR(:, ids,:,:), 2)); 

       end
    end
    
    sub2exc2 = find(nTR<minSubCrit); 
    sub2exc = union(sub2exc, sub2exc2);
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

    [tObs id] = max(allTObs');
    
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

%% plot nicely

cond2plot = 'CC'; %allT CC CI
minSubCrit = 5; 


if strcmp(cfg.brainROI, 'vvs')
    sub2exc = [18 22];
elseif strcmp(cfg.brainROI, 'pfc')
    sub2exc = [1];
end

tiledlayout(8,8, 'TileSpacing', 'compact', 'Padding', 'compact');
if strcmp(cfg.period(1), 'M')
    set(gcf, 'Position', [100 100 1800 1300])
else
    set(gcf, 'Position', [100 100 670 1300])
end

for layi = 1:size(nnFit{2}, 1)
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
                ids = ids(:,20)==1; 
            elseif strcmp(cond2plot, 'IC')
                ids = ids(:,20)==0; 
            elseif strcmp(cond2plot, 'CI')
                ids = ids(:,19)==1; 
            elseif strcmp(cond2plot, 'II')
                ids = ids(:,19)==0; 
            elseif strcmp(cond2plot, 'allT')
                ids = logical(ones(1, size(ids, 1))); %just takes them all
            end

             nTR(subji,:) = sum(ids==1); 
             nnH(subji, : ,:) = squeeze(mean(avTR(:, ids,:,:), 2)); 

       end
    end
    
    sub2exc2 = find(nTR<minSubCrit); 
    sub2exc = union(sub2exc, sub2exc2);
    nnH(sub2exc, :, :) = []; 
    nnH = squeeze(nnH);
    [h p ci ts] = ttest(nnH);
    h = squeeze(h); t = squeeze(ts.tstat); 
    %h(:, 1:5) = 0; % only sum p-values in clusters after the baseline
    
    d2p = squeeze(mean(nnH, 'omitnan'));
    
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

    [tObs id] = max(allTObs');
    
    freqs = 1:520; 
    if strcmp(cfg.period(1), 'M')
        times = 1:size(t, 2)*10; 
        if strcmp(cfg.net2load, 'BLNETi') & strcmp(cfg.brainROI, 'pfc') & strcmp(cond2plot, 'CC')
            h = zeros(52, 40);
            if layi == 7 
                h(clustinfo.PixelIdxList{8}) = 1;
            end
        end
        if strcmp(cond2plot, 'IC')
            h = zeros(52, 40);
        end
        if strcmp(cfg.net2load, 'BLNETi') & strcmp(cfg.brainROI, 'vvs') & strcmp(cond2plot, 'CC')
            h = zeros(52, 40);
            if layi == 6
                h(clustinfo.PixelIdxList{12}) = 1;
            end
        end
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

%exportgraphics(gcf, [paths.results.DNNs 'myP.png'], 'Resolution', 300); 
exportgraphics(gcf, ['myP.png'], 'Resolution', 300); 
close all



%% COMPUTE CLUSTERS in each permutation FOR THE TRIAL LEVEL FITS
% % % first load because nnFitPerm is big
clc
%Network_ROI_ER_layers_freqs_avRepet_avTFV_fRes(0-1)_fitMode(0:noTrials; 1:Trials)__timeRes__win__mf

% use the same name as the previously plotted file
f2sav = [f2sav '_200p.mat'];
%f2sav = [f2sav '_100p.mat'];
cfg = getParams(f2sav);
paths = load_paths_WM(cfg.brainROI, cfg.net2load);

cd ([paths.results.DNNs])
load(f2sav);

%% no need to change the condition to plot

cond2plot = 'CC'; %allT CC CI

if strcmp(cfg.brainROI, 'pfc')
    sub2exc = [1];
elseif strcmp(cfg.brainROI, 'vvs')
    sub2exc = [18 22];
end

f2t = strsplit(f2sav, '_');
nPerm = double(string((f2t{end}(1:end-5))));
nLays = size(nnFitPerm{1, 2}, 1);

for permi = 1:nPerm
    
    nnFitP = nnFitPerm(permi, :);
    for layi = 1:nLays
        
        clear nnHPerm
        for subji = 1:length(nnFitP)
           if ~isempty(nnFitP{subji})
               if strcmp(cfg.period(1), 'M')
                 avTR = atanh(nnFitP{subji}(layi,:,:,:));
                 ids = cellfun(@(x) strsplit(x, ' '), nnFit{subji, 2}, 'un', 0); %str2num does not work
                 ids = double(string(cat(1, ids{:})));%str2num does not work
               elseif strcmp(cfg.period(1), 'E') % for clarity
                 avTR = atanh(nnFitP{subji}(layi,:,:,1:15));
                 ids= str2num(cell2mat(nnFit{subji, 2}));
               end
                
                if strcmp(cond2plot, 'CC')
                    ids = ids(:,20)==1; 
                elseif strcmp(cond2plot, 'IC')
                    ids = ids(:,20)==0; 
                elseif strcmp(cond2plot, 'CI')
                    ids = ids(:,19)==1; 
                elseif strcmp(cond2plot, 'II')
                    ids = ids(:,19)==0; 
                elseif strcmp(cond2plot, 'allT')
                    ids = logical(ones(1, size(ids, 1))); %just takes them all
                end
    
                 nTR(subji,:) = sum(ids==1); 
                 nnHPerm(subji, : ,:) = squeeze(mean(avTR(:, ids,:,:), 2)); 
    
           end
        end
    
        sub2exc2 = find(nTR<minSubCrit); 
        sub2exc = union(sub2exc, sub2exc2);
        nnHPerm(sub2exc, :, :) = []; 
        nnHPerm = squeeze(nnHPerm);
        %[h p ci ts] = ttest(nnH, 0, "Tail","right");
        [h p ci ts] = ttest(nnHPerm);
        h = squeeze(h); t = squeeze(ts.tstat); 
        %h(:, 1:5) = 0; % only sum p-values in clusters after the baseline
        
        d2p = squeeze(mean(nnHPerm, 'omitnan'));
        
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

for layi = 1:size(tObs, 2)
    clear mcsR mcsP
    mcsR = tObs(layi); 
    mcsP = squeeze(max_clust_sum_perm(:,layi));
    
    %allAb = mcsP(abs(mcsP) > abs(mcsR));
    allAb = mcsP(mcsP > mcsR);
    p(layi, :) = 1 - ((nPerm-1) - (length (allAb)))  / nPerm;
    
end

p_ranked = p; p_ranked(isnan(p_ranked)) = []; 
p_ranked = sort(p_ranked(:));
p

%% p last layer only
p = p (end,:);
p_ranked = p; p_ranked(p_ranked == 0 | isnan(p_ranked)) = []; 
p_ranked = sort(p_ranked(:))







%% load files to plot all CATEGORY RESULTS  -- IN one PLOT ONLY 

clear, clc 

listF2sav = {       
    
                    % 'CAT_pfc_E123_[1]_3-8_1_0_0_0_.1_5_1.mat';
                    % 'CAT_pfc_E123_[1]_9-12_1_0_0_0_.1_5_1.mat';
                    % 'CAT_pfc_E123_[1]_13-29_1_0_0_0_.1_5_1.mat';
                    % 'CAT_pfc_E123_[1]_30-38_1_0_0_0_.1_5_1.mat';
                    % 'CAT_pfc_E123_[1]_39-54_1_0_0_0_.1_5_1.mat';


                    'CAT_pfc_M123_[1]_3-8_1_0_0_0_.1_5_1.mat';
                    'CAT_pfc_M123_[1]_9-12_1_0_0_0_.1_5_1.mat';
                    'CAT_pfc_M123_[1]_13-29_1_0_0_0_.1_5_1.mat';
                    'CAT_pfc_M123_[1]_30-38_1_0_0_0_.1_5_1.mat';
                    'CAT_pfc_M123_[1]_39-54_1_0_0_0_.1_5_1.mat';


%                     'CAT_vvs_E123_[1]_3-8_1_0_0_0_.1_5_1.mat';
%                     'CAT_vvs_E123_[1]_9-12_1_0_0_0_.1_5_1.mat';
%                     'CAT_vvs_E123_[1]_13-29_1_0_0_0_.1_5_1.mat';
%                     'CAT_vvs_E123_[1]_30-38_1_0_0_0_.1_5_1.mat';
%                     'CAT_vvs_E123_[1]_39-54_1_0_0_0_.1_5_1.mat';

                    % 'CAT_vvs_M123_[1]_3-8_1_0_0_0_.1_5_1.mat';
                    % 'CAT_vvs_M123_[1]_9-12_1_0_0_0_.1_5_1.mat';
                    % 'CAT_vvs_M123_[1]_13-29_1_0_0_0_.1_5_1.mat';
                    % 'CAT_vvs_M123_[1]_30-38_1_0_0_0_.1_5_1.mat';
                    % 'CAT_vvs_M123_[1]_39-54_1_0_0_0_.1_5_1.mat';

                    % 'BLNETi_pfc_M123_[8-8-56]_3-8_1_0_0_0_.1_5_1'
                    % 'BLNETi_pfc_M123_[8-8-56]_9-12_1_0_0_0_.1_5_1'
                    % 'BLNETi_pfc_M123_[8-8-56]_13-29_1_0_0_0_.1_5_1'
                    % 'BLNETi_pfc_M123_[8-8-56]_30-38_1_0_0_0_.1_5_1'
                    % 'BLNETi_pfc_M123_[8-8-56]_39-54_1_0_0_0_.1_5_1'


             };   

for listi = 1:length(listF2sav)
    disp(['File > ' num2str(listi) '      ' listF2sav{listi}]);
    clearvars -except listF2sav listi allFnnH
        
    f2sav       = listF2sav{listi}; 
    cfg = getParams(f2sav);
    paths = load_paths_WM(cfg.brainROI,  cfg.net2load);
    filelistSess = getFiles([paths.results.DNNs]);
    
    for sessi= 1:length(filelistSess) %this one starts at 1 and not at 3
        load([paths.results.DNNs listF2sav{listi}]);   
        
        nnH = cat(1, nnFit{:,1});
        
    end

    allFnnH(listi, :, :) = nnH; 


end
cfg = getParams(f2sav);
if strcmp(cfg.brainROI, 'vvs')
    %sub2exc = [18 22];
    sub2exc = [18];
elseif strcmp(cfg.brainROI, 'pfc')
    sub2exc = [];
end
%allFnnH(:, sub2exc, :) = []; 

%% plot


for freqi = 1:size(allFnnH, 1)
    if strcmp(cfg.period(1), 'M')
        nnH = allFnnH(freqi, :, 1:40);
    else
        nnH = allFnnH(freqi, :, 1:15);
    end
    nnH = squeeze(nnH);
    nnH(sub2exc, :, :) = []; 

    [h p ci ts] = ttest(nnH);
    h = squeeze(h); t = squeeze(ts.tstat);
    clustinfo = bwconncomp(h);
    
    if ~isempty(clustinfo.PixelIdxList)
        for pixi = 1:length(clustinfo.PixelIdxList)
             %if length(clustinfo.PixelIdxList{pixi}) > 1
                allTObs(freqi, pixi, :) = sum(t(clustinfo.PixelIdxList{pixi}));
             %end        
        end
    else
        allTObs(freqi, :, :) = 0;
    end
    
    hb = h; hb(h==0) = nan; hb(hb==1) = 0; 
    clear hbl
    hbL(freqi, :) = hb; 
    
    mART(freqi, :) = squeeze(mean(nnH)); 
    stdART = squeeze(std(nnH)); 
    seART = stdART/ sqrt(size(nnH, 1));
    

end


if strcmp(cfg.period(1), 'M')
    set(gcf, 'Position', [100 100 500 350])
    myCmap = colormap(brewermap(14,'RdPu'));
    myCmap = myCmap([6 8 10 12 14],:)
    times = 1:40;
    x = (-.015:-.003:-.028)';
    %hbL([1:5], :) = nan;     
    hbL = hbL+x;
    plot (times, hbL, 'Linewidth', 4); hold on; 
    plot(times, mART, 'Linewidth', 3); hold on; 
    set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', [], 'xlim', [1 37]); 
    set(gca, 'FontSize', 12, 'ylim', [-.0375 .0375]);
    plot([3 3],get(gca,'ylim'), 'k:','lineWidth',2);
    plot(get(gca,'xlim'), [0 0],'k:','lineWidth',2);

elseif strcmp(cfg.period(1), 'E')
    
    figure()
    set(gcf, 'Position', [100 100 150 300])
    myCmap = colormap(brewermap(14,'RdPu'));
    myCmap = myCmap([6 8 10 12 14],:)
    %colormap(jet(5));
    times = 1:15;
    plot(times, mART, 'Linewidth', 5); hold on; 
    x = (-.015:-.0053:-.039)';
    %hbL(4, 1:3) = nan; 

    hbL = hbL+x;
    plot (times, hbL, 'Linewidth', 5); hold on; 
    set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
    set(gca, 'FontSize', 12, 'ylim', [-.04 .15]);
    plot([3 3],get(gca,'ylim'), 'k:','lineWidth',3);
    plot(get(gca,'xlim'), [0 0],'k:','lineWidth',3);
    %legend

end


set(gca, 'ColorOrder', myCmap)


%exportgraphics(gcf, [paths.results.DNNs 'myP.png'], 'Resolution', 300); 
exportgraphics(gcf, ['myP.png'], 'Resolution', 300); 


%% PLOT BANDS ONE BY ONE - CATEGORY MODEL 

clear, clc 
   
f2sav = 'CAT_pfc_E123_[1]_3-8_1_0_0_0_.1_5_1.mat';


cfg = getParams(f2sav);
paths = load_paths_WM(cfg.brainROI,  cfg.net2load);

load([paths.results.DNNs f2sav]);   

cfg = getParams(f2sav);
if strcmp(cfg.brainROI, 'vvs')
    sub2exc = [18 22];
elseif strcmp(cfg.brainROI, 'pfc')
    sub2exc = [1];
end


for subji = 1:length(nnFit)
     if ~isempty(nnFit{subji, 1})
       if strcmp(cfg.period(1), 'M')
         nnH(subji, : ) = atanh(nnFit{subji, 1}(:,1:40));
       else
         nnH(subji, : ) = atanh(nnFit{subji, 1}(:,1:15));
       end
     end
end


nnH = squeeze(nnH);
nnH(sub2exc, :, :) = []; 

[h p ci ts] = ttest(nnH);
h = squeeze(h); t = squeeze(ts.tstat);
clustinfo = bwconncomp(h);

if ~isempty(clustinfo.PixelIdxList)
    for pixi = 1:length(clustinfo.PixelIdxList)
         %if length(clustinfo.PixelIdxList{pixi}) > 1
            allTObs(pixi, :) = sum(t(clustinfo.PixelIdxList{pixi}));
         %end        
    end
else
    allTObs = 0;
end

tObs = max(allTObs);
hb = h; hb(h==0) = nan; hb(hb==1) = 0; 
mART = squeeze(mean(nnH)); 
stdART = squeeze(std(nnH)); 
seART = stdART/ sqrt(size(nnH, 1));





if strcmp(cfg.period(1), 'M')
    set(gcf, 'Position', [100 100 500/1.5 300/1.7])
    myCmap = colormap(brewermap(14,'*RdPu'));
    times = 1:40;
    x = (-.015:-.003:-.028)';
    
    plot (times, hb, 'Linewidth', 4); hold on; 
    plot(times, mART, 'Linewidth', 3); hold on; 
    set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', [], 'xlim', [1 37]); 
    set(gca, 'FontSize', 12, 'ylim', [-.0375 .0375]);
    plot([5 5],get(gca,'ylim'), 'k:','lineWidth',2);
    plot(get(gca,'xlim'), [0 0],'k:','lineWidth',2);

elseif strcmp(cfg.period(1), 'E')
    
    figure()
    set(gcf, 'Position', [100 100 150 300])
    myCmap = colormap(brewermap(14,'RdPu'));
    myCmap = myCmap([6 8 10 12 14],:);
    %colormap(jet(5));
    times = 1:15;
    plot(times, mART, 'Linewidth', 5); hold on; 
    x = (-.015:-.0053:-.039)';
    
    plot (times, hb, 'Linewidth', 5); hold on; 
    set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', [], 'xlim', [1 11]); 
    set(gca, 'FontSize', 12, 'ylim', [-.04 .15]);
    plot([5 5],get(gca,'ylim'), 'k:','lineWidth',3);
    plot(get(gca,'xlim'), [0 0],'k:','lineWidth',3);
    %legend

end


set(gca, 'ColorOrder', myCmap)


%exportgraphics(gcf, [paths.results.DNNs 'myP.png'], 'Resolution', 300); 
exportgraphics(gcf, ['myP.png'], 'Resolution', 300); 



%% compute clusters in each permutation BANDS for all layers
%Network_ROI_ER_layers_freqs_avRepet_avTFV_fRes(0-1)_fitMode(0:noTrials; 1:Trials)__timeRes__win__mf
%f2sav =  'RNN_pfc_M123_[8-8-56]_13-29_1_1_0_0_.1_5_1_1000p.mat'; 
%f2sav = 'CORrt_vvs_E123_[1-8]_39-54_1_0_0_0_.1_5_1_1000p.mat';
%f2sav = 'CAT_pfc_E123_[1]_3-8_1_0_0_0_.1_5_1_1000p.mat';


f2sav = [f2sav '_1000p.mat'];

cd ([paths.results.DNNs])
load(f2sav);


if strcmp(cfg.brainROI, 'pfc')
    sub2exc = [1];
elseif strcmp(cfg.brainROI, 'vvs')
    sub2exc = [18 22];
end


f2t = strsplit(f2sav, '_');
nPerm = double(string((f2t{13}(1:end-5))));


for permi = 1:nPerm
    
    dataP = squeeze(nnFitPerm(permi, :, :));
    dataP(sub2exc, :,:) = []; 
    [h p ci ts] = ttest(dataP);
    h = squeeze(h); t = squeeze(ts.tstat);

    
    clear allSTs  
    clustinfo = bwconncomp(h);
    for pxi = 1:length(clustinfo.PixelIdxList)
       allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
    end
    
    %sort all clusters 
    if exist('allSTs')
        [max2u id] = max(allSTs);
        max_clust_sum_perm(permi, :) = allSTs(id); 
    else
        max_clust_sum_perm(permi, :) = 0; 
    end


end


%allAB = max_clust_sum_perm(abs(max_clust_sum_perm) > abs(tObs));
allAB = max_clust_sum_perm(max_clust_sum_perm > tObs);
p = 1 - ((nPerm-1) - (length (allAB)))  / nPerm



cd (paths.github)



%% plot final figure only last layer / time point ENCODING
times = 1:150;
freqs = 1:520; 
h = zeros(52, 15); 

%h(clustinfo.PixelIdxList{1}) = 1; % Alexnet VVS Layer 1-3, 5

%h(clustinfo.PixelIdxList{1}) = 1; % Alexnet VVS Layer 4 and 6
%h(clustinfo.PixelIdxList{2}) = 1; % Alexnet VVS Layer 4 and 6


%h(clustinfo.PixelIdxList{5}) = 1; % BLNETi PFC
%h(clustinfo.PixelIdxList{1}) = 1; % BLNETi VVS

%h(clustinfo.PixelIdxList{1}) = 1; % BLNETe cluster 1
%h(clustinfo.PixelIdxList{2}) = 1; % BLNETe cluster 2

%h(clustinfo.PixelIdxList{2}) = 1; %Cornet
%h(clustinfo.PixelIdxList{4}) = 1; %Cornet cluster 2

%h(clustinfo.PixelIdxList{6}) = 1; % ITEM MODEL VVS
h(clustinfo.PixelIdxList{6}) = 1; % ITEM MODEL PFC

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
h = zeros(52, 40); 

%h(clustinfo.PixelIdxList{2}) = 1; %BLNETi PFC

%h(clustinfo.PixelIdxList{3}) = 1; %CORNET PFC
%h(clustinfo.PixelIdxList{15}) = 1; %CORNET VVS

%h(clustinfo.PixelIdxList{5}) = 1; %category model

%h(clustinfo.PixelIdxList{2}) = 1; %pfc - Cornet
%h(clustinfo.PixelIdxList{23}) = 1; %pfc - Cornet

%h(clustinfo.PixelIdxList{8}) = 1; %pfc - ecoset
%h(clustinfo.PixelIdxList{10}) = 1; %vvs1 - ecoset
%h(clustinfo.PixelIdxList{23}) = 1; %vvs2 - ecoset

h(clustinfo.PixelIdxList{27}) = 1; %BLNETi PFC

figure; set(gcf, 'Position', [1000 918 560 420])
myCmap = colormap(brewermap([],'*Spectral'));
colormap(myCmap)
%contourf(times, freqs, t, 100, 'linecolor', 'none'); hold on; %colorbar
%contour(times, freqs, h, 1, 'Color', [0, 0, 0], 'LineWidth', 2);
contourf(times, freqs, myresizem(t, 10), 100, 'linecolor', 'none'); hold on; %colorbar
contour(times, freqs, myresizem(h, 10), 1, ':', 'Color', [0, 0, 0], 'LineWidth', 4);

set(gca, 'ytick', [], 'yticklabels', [], 'xtick', [], 'xticklabels', []); 
set(gca, 'xlim', [1 390], 'clim', [-4 4], 'FontSize', 10);
%set(gca, 'clim', [-4 4], 'FontSize', 10);
plot([45 45],get(gca,'ylim'), 'k:','lineWidth', 5); 
%colorbar

%exportgraphics(gcf, [paths.results.DNNs 'myP.png'], 'Resolution', 300); 
exportgraphics(gcf, ['myP.png'], 'Resolution', 300); 




%%  plot all layers FREQUENCY RESOLVED FANCY PLOT FOR MASK ANALYSIS
%Network_ROI_ER_layers_freqs_avRepet_avTFV_fRes(0-1)_fitMode(0:noTrials; 1:Trials)_timeRes_win_mf
clear , clc
f2sav = 'BLNETi_pfc_M123_[8-8-56]_3-54_0_0_1_0_.1_5_1_MASK'; 
%f2sav = 'BLNETi_pfc_E123_[8-8-56]_3-54_0_0_1_0_.1_5_1_MASK'; 
%f2sav = 'Alex_pfc_E123_[1-8]_3-54_0_0_1_0_.1_5_1_MASK'; 

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
        
        
        if strcmp(cfg.net2load, 'BLNETi') & strcmp(cfg.brainROI, 'vvs') & layi == 5
            h(clustinfo.PixelIdxList{23}) = 1;
        end
        if strcmp(cfg.net2load, 'BLNETi') & strcmp(cfg.brainROI, 'vvs') & layi == 6
            h(clustinfo.PixelIdxList{17}) = 1;
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
        if strcmp(cfg.net2load, 'Alex') & strcmp(cfg.brainROI, 'vvs') 
            if layi == 1
                h(clustinfo.PixelIdxList{1}) = 1;
            end
            if layi ==5
                h(clustinfo.PixelIdxList{4}) = 1;
            end
            if layi == 6
                h(clustinfo.PixelIdxList{5}) = 1;
            end
            if layi == 7 
                h(clustinfo.PixelIdxList{3}) = 1;
            end
        end
        if strcmp(cfg.net2load, 'Alex') & strcmp(cfg.brainROI, 'pfc') 
            if layi == 1
                h(clustinfo.PixelIdxList{4}) = 1;
            end
            
        end

        if strcmp(cfg.net2load, 'BLNETi') & strcmp(cfg.brainROI, 'vvs') 
            if layi == 6
                h(clustinfo.PixelIdxList{4}) = 1;
            end
            if layi == 7 
                h(clustinfo.PixelIdxList{5}) = 1;
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

%exportgraphics(gcf, [paths.results.DNNs 'myP.png'], 'Resolution', 300); 
exportgraphics(gcf, ['myP.png'], 'Resolution', 300); 
close all; 








%%










