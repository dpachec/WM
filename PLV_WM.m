%% PLV inter area (PFV-VVS)
%% 

clearvars 

currentF = pwd;
vvs_link = 'D:\_WM\analysis\out_contrasts\raw_traces\allTrials\vvs';
pfc_link = 'D:\_WM\analysis\out_contrasts\raw_traces\allTrials\pfc';


pfc_ids = [2 3  5  9 10 11 12 14 15 16];
vvs_ids = [7 9 13 18 19 20 21 23 27 28];

tic
for subji = 1:10 
    subji
    cd(vvs_link);sublist = dir('*contr.mat');sublist = {sublist.name};
    load (sublist{vvs_ids(subji)})
    c_vvs = cfg_contrasts;
    cd(pfc_link);sublist = dir('*contr.mat');sublist = {sublist.name};
    load (sublist{pfc_ids(subji)})
    c_pfc = cfg_contrasts;
    count = 0;
    for triali = 1:size(c_vvs.oneListTraces,3)
        id2u = strsplit(c_vvs.oneListIds_c{triali});
        if strcmp(id2u{1},'7') & strcmp(id2u{2},'4') 
            count = count+1;
            for chani = 1:size(c_vvs.chanNames,1)
                for chanj = 1:size(c_pfc.chanNames,1)
                   t_vvs = squeeze(c_vvs.oneListTraces(chani,2000:3000,triali));
                   t_pfc = squeeze(c_pfc.oneListTraces(chanj,2000:3000,triali));

                   EEG_vvs.data    = t_vvs;
                   EEG_vvs.trials  = 1; EEG_vvs.srate   = 1000; EEG_vvs.nbchan  = 1; EEG_vvs.pnts = size(t_vvs,2);EEG_vvs.event   = [];
                   EEG_vvs         = pop_eegfiltnew (EEG_vvs, 16,29);
                   data_vvs        = squeeze(EEG_vvs.data); 
                   dataHA_vvs      = angle(hilbert(data_vvs));

                   EEG_pfc.data    = t_pfc;
                   EEG_pfc.trials  = 1; EEG_pfc.srate   = 1000; EEG_pfc.nbchan  = 1; EEG_pfc.pnts = size(t_vvs,2);EEG_pfc.event   = [];
                   EEG_pfc         = pop_eegfiltnew (EEG_pfc, 16,29);
                   data_pfc        = squeeze(EEG_pfc.data); 
                   dataHA_pfc      = angle(hilbert(data_pfc));


                   phase_data(1,:) = dataHA_vvs;
                   phase_data(2,:) = dataHA_pfc;

                   PLVch{subji,:}(count, chani, chanj ) = abs(mean(exp(1i*(diff(phase_data,1)))));
                end
            end
        end
    end

    
end


cd (currentF)
save('PLVch', 'PLVch');



toc

%% compute mean
clearvars -except PLVch
for subji = 1:10
    plv2u = PLVch{subji};
    for triali = 1:size(plv2u, 1)
        plvTR(triali,:) = mean(plv2u(triali, :, :), 'all');
    end 
    plvSUBJ(subji, :) = mean(plvTR);
end 

%% 2Bar 
data.data = [plvSUBJ_MI plvSUBJ_SI]; 


figure(2); set(gcf,'Position', [0 0 500 650]); 
mean_S = mean(data.data, 1);
std_S = std(data.data, [], 1);
se_S = std_S / sqrt(13); 
hb = plot ([1 2], data.data); hold on;
set(hb, 'lineWidth', 3, 'Marker', '.', 'MarkerSize',35);hold on;
h = bar (mean_S);hold on;
%h1=errorbar(mean_S, se_S,'c'); hold on;
set(h,'FaceColor', 'none', 'lineWidth', 3);
%set(h1, 'Color','k','linestyle','none', 'lineWidth', 2);
set(gca,'XTick',[1 2],'XTickLabel',{'', ''}, ...
    'FontSize', 30, 'linew',2, 'xlim', [0 3], 'ylim', [.18 .5] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 3);

[h p ci t] = ttest (data.data(:,1), data.data(:,2));
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth', 3);






%% TIME FREQUENCY reSOLVED

clearvars 

currentF = pwd;
vvs_link = 'D:\_WM\analysis\out_contrasts\raw_traces\allTrials\vvs';
pfc_link = 'D:\_WM\analysis\out_contrasts\raw_traces\allTrials\pfc';


pfc_ids = [2 3  5  9 10 11 12 14 15 16];
vvs_ids = [7 9 13 18 19 20 21 23 27 28];

tic
for subji = 1:10 
    subji
    cd(vvs_link);sublist = dir('*contr.mat');sublist = {sublist.name};
    load (sublist{vvs_ids(subji)})
    c_vvs = cfg_contrasts;
    cd(pfc_link);sublist = dir('*contr.mat');sublist = {sublist.name};
    load (sublist{pfc_ids(subji)})
    c_pfc = cfg_contrasts;
    count = 0;
    for triali = 1:size(c_vvs.oneListTraces,3)
        id2u = strsplit(c_vvs.oneListIds_c{triali});
        if strcmp(id2u{1},'7') & strcmp(id2u{2},'4') 
            count = count+1;
            for chani = 1:size(c_vvs.chanNames,1)
                for chanj = 1:size(c_pfc.chanNames,1)
                   t_vvs = squeeze(c_vvs.oneListTraces(chani,2000:3000,triali));
                   t_pfc = squeeze(c_pfc.oneListTraces(chanj,2000:3000,triali));

                   EEG_vvs.data    = t_vvs;
                   EEG_vvs.trials  = 1; EEG_vvs.srate   = 1000; EEG_vvs.nbchan  = 1; EEG_vvs.pnts = size(t_vvs,2);EEG_vvs.event   = [];
                   EEG_vvs         = pop_eegfiltnew (EEG_vvs, 16,29);
                   data_vvs        = squeeze(EEG_vvs.data); 
                   dataHA_vvs      = angle(hilbert(data_vvs));

                   EEG_pfc.data    = t_pfc;
                   EEG_pfc.trials  = 1; EEG_pfc.srate   = 1000; EEG_pfc.nbchan  = 1; EEG_pfc.pnts = size(t_vvs,2);EEG_pfc.event   = [];
                   EEG_pfc         = pop_eegfiltnew (EEG_pfc, 16,29);
                   data_pfc        = squeeze(EEG_pfc.data); 
                   dataHA_pfc      = angle(hilbert(data_pfc));


                   phase_data(1,:) = dataHA_vvs;
                   phase_data(2,:) = dataHA_pfc;

                   PLVch{subji,:}(count, chani, chanj ) = abs(mean(exp(1i*(diff(phase_data,1)))));
                end
            end
        end
    end

    
end


cd (currentF)
save('PLVch', 'PLVch');



toc



















%% 


clear
load pfc_elec
load RNN_pfc_M&C_noAv_54_Real_T_1-49_C_aFV
load PLVch_1-8Hz_2500-5500

pfc_ids = [2 3  5  9 10 11 12 14 15 16];
pfc = all_r_Times_Trials(pfc_ids)';
tmp = cell2mat(cellfun(@size, chanNames_all, 'un', 0));
s2e_pfc = ~(tmp(pfc_ids,1) > 1); 
pfc(s2e_pfc) = [];
PLVch(s2e_pfc) = [];
x = cellfun(@size, pfc, 'un', 0)

clear allR
for subji = 1:length(pfc)
   
    plvT = squeeze(mean(mean(PLVch{subji},3),2));
    plvT(plvT==0) = [];
    tmp = squeeze(pfc{subji}(7,:,30:52,13:25)); %last layer in beta and time period of interest
    tmp = squeeze(mean(mean(tmp,3),2));
    ch2(subji,:) = mean(tmp);
    %tmp = reshape (tmp, size(tmp, 1), []);
    
    
    allR(subji,:) = corr(plvT, tmp, 'type','s');
         
    
end

[h p ci ts] = ttest(allR);
disp(['t >>  ' num2str(ts.tstat)]);



%% 





















%% build 2 compare each frequency
clearvars 
currentF = pwd;
vvs_link = 'D:\_WM\analysis\out_contrasts\raw_traces\vvs\allTrials';
pfc_link = 'D:\_WM\analysis\out_contrasts\raw_traces\pfc\allTrials';

pfc_ids = [2 3  5  9 10 11 12 14 15 16];
vvs_ids = [7 9 13 18 19 20 21 23 27 28];

freqs2use = [3:29 30:5:150];

tic
for subji = 1:10 
    subji
    cd(vvs_link);sublist = dir('*contr.mat');sublist = {sublist.name};
    load (sublist{vvs_ids(subji)})
    c_vvs = cfg_contrasts;
    cd(pfc_link);sublist = dir('*contr.mat');sublist = {sublist.name};
    load (sublist{pfc_ids(subji)})
    c_pfc = cfg_contrasts;
    count = 1;
    for triali = 1:size(c_vvs.oneListTraces,3)
        id2u = strsplit(c_vvs.oneListIds_c{triali});
        if strcmp(id2u{1},'7') & ~strcmp(id2u{2},'4') 
            for freqi = 1:length(freqs2use)
                f = freqs2use(freqi);
                for chani = 1:size(c_vvs.chanNames,1)
                    for chanj = 1:size(c_pfc.chanNames,1)
                       t_vvs = squeeze(c_vvs.oneListTraces(chani,2500:5500,triali));
                       t_pfc = squeeze(c_pfc.oneListTraces(chanj,2500:5500,triali));

                       EEG_vvs.data    = t_vvs;
                       EEG_vvs.trials  = 1; EEG_vvs.srate   = 1000; EEG_vvs.nbchan  = 1; EEG_vvs.pnts = size(t_vvs,2);EEG_vvs.event   = [];
                       EEG_vvs         = pop_eegfiltnew (EEG_vvs, f,f);
                       data_vvs        = squeeze(EEG_vvs.data); 
                       phase_data(1,:)      = angle(hilbert(data_vvs));

                       EEG_pfc.data    = t_pfc;
                       EEG_pfc.trials  = 1; EEG_pfc.srate   = 1000; EEG_pfc.nbchan  = 1; EEG_pfc.pnts = size(t_vvs,2);EEG_pfc.event   = [];
                       EEG_pfc         = pop_eegfiltnew (EEG_pfc, f,f);
                       data_pfc        = squeeze(EEG_pfc.data); 
                       phase_data(2,:) = angle(hilbert(data_pfc));

                       PLVch{subji}(freqi, count, chani, chanj ) = abs(mean(exp(1i*(diff(phase_data,1)))));
                    end
                end
            end
            count = count+1;
        end
    end

    
end


cd (currentF)
save('PLVch_allF_2500-5500', 'PLVch');



toc



%% plot for each freq
load pfc_elec
load RNN_pfc_M&C_noAv_54_Real_T_1-49_C_aFV
load PLVch_allF_2500-5500

pfc_ids = [2 3  5  9 10 11 12 14 15 16];
pfc = all_r_Times_Trials(pfc_ids)';
tmp = cell2mat(cellfun(@size, chanNames_all, 'un', 0));
s2e_pfc = ~(tmp(pfc_ids,1) > 1); 
pfc(s2e_pfc) = [];
PLVch(s2e_pfc) = [];


clear allR
for subji = 1:length(pfc)
   
    for freqi = 1:52
        conT = squeeze(mean(mean(PLVch{subji}(freqi,:,:,:),4),3))';
        tmp = squeeze(pfc{(subji)}(1,:,11:27,7:12)); %last layer in beta and time period of interest
        tmp = squeeze(mean(mean(tmp,3),2));
        %ch2(subji,:) = mean(tmp);
        %tmp = reshape (tmp, size(tmp, 1), []);
        allR(subji,freqi) = corr(conT, tmp, 'type','s');
    end 
    
end

[h p ci ts] = ttest(allR)
 
plot(h)


%% plot for each freq 2
load pfc_elec
load RNN_pfc_M&C_noAv_54_Real_T_1-49_C_aFV
load PLVch_allF_2500-5500
%load PLVch_allF

pfc_ids = [2 3  5  9 10 11 12 14 15 16];
pfc = all_r_Times_Trials(pfc_ids)';
tmp = cell2mat(cellfun(@size, chanNames_all, 'un', 0));
s2e_pfc = ~(tmp(pfc_ids,1) > 2); 
s2e_pfc(4) = 1;
pfc(s2e_pfc) = [];
PLVch(s2e_pfc) = [];



clear allR
for subji = 1:length(pfc)
   for freqi = 1:52
    for freqi2 = 1:52
        conT = squeeze(mean(mean(PLVch{subji}(freqi,:,:,:),4),3))';
        %tmp = squeeze(pfc{(subji)}(1,:,freqi2:freqi2,7:12)); %last layer in beta and time period of interest
        
        %nornmalized to baseline
        tmp = squeeze(pfc{(subji)}(1,:,freqi2:freqi2,7:12)); %last layer in beta and time period of interest
        %tmpb = mean(squeeze(pfc{(subji)}(1,:,freqi2:freqi2,1:5)),2); %last layer in beta and time period of interest
        %tmp = tmp-tmpb;
        
        tmp = squeeze(mean(mean(tmp,3),2));
        %ch2(subji,:) = mean(tmp);
        %tmp = reshape (tmp, size(tmp, 1), []);
        allR(subji,freqi, freqi2) = corr(conT, tmp, 'type','s');
    end 
   end
end

[h p ci ts] = ttest(allR)
 

figure()
contourf(myresizem(squeeze(ts.tstat), 10), 40, 'linecolor', 'none'); axis equal, hold on; colorbar
contour(myresizem(squeeze(h), 10), 1, 'Color', [0, 0, 0], 'LineWidth', 2);



%% 
d2p = squeeze(PLVch{1}(22,:,:));
figure();
imagesc(d2p); colorbar



%% load trialIDs

d2p = squeeze(all_r_Times(16,56,:,:));
imagesc(d2p)



%% all_r-times rho values show expected increases


d2p = squeeze(mean(all_r_Times(pfc_ids,56,:,:),'omitnan'));
figure()
imagesc(d2p); colorbar


d2p2 = all_r_Times_Trials(pfc_ids)









































%%