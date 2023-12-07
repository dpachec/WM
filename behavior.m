%% check behavioral data
% first load events

clear

[all_events allTrlInfo ] = loadLogsWM;
allTrlInfo = [allTrlInfo{:}];%% has to be done twice because of the nested 
allTrlInfo = [allTrlInfo{:}];%% has to be done twice because of the nested 
allTrlInfo = allTrlInfo(~cellfun('isempty',allTrlInfo));
allTrlInfo = allTrlInfo';


%% check number of single and multi item trials for every session

for sessi = 1:length(allTrlInfo)
    
    trlinfo = allTrlInfo{sessi};
    nSITR(sessi, :) = length(trlinfo(trlinfo(:, 10) ~= 4));
    nMITR(sessi, :) = length(trlinfo(trlinfo(:, 10) == 4));


end

nSITR
%% correct item-cat and by pos

clearvars -except allTrlInfo
tic 

for vp=1:length(allTrlInfo)
    
    trialInfo = allTrlInfo{vp};
    
    cues = trialInfo(:,10);
    
    indall = find(cues == 4);
    ind123 = find(cues ~= 4);
    
    %all
    corr_all(1)=sum(trialInfo(indall,40))./numel(indall); %item
    corr_all(2)=sum(trialInfo(indall,41))./numel(indall); %category
    gData.corr_all(vp, :) = corr_all;
    
    %123
    corr_123(1)=sum(trialInfo(ind123,40))./numel(ind123); %item
    corr_123(2)=sum(trialInfo(ind123,41))./numel(ind123); %category
    gData.corr_123(vp, :) = corr_123;
    
    %all trials in each position 
    corr_all_pos(1,:)=nansum(trialInfo(indall,34:36))./numel(indall);
    corr_all_pos(2,:)=nansum(trialInfo(indall,37:39))./numel(indall);
    gData.corr_all_pos{vp} = corr_all_pos;
    %123 trials in each position 
    corr_123_pos(1,:)=nansum(trialInfo(ind123,34:36))./(numel(ind123)./3);
    corr_123_pos(2,:)=nansum(trialInfo(ind123,37:39))./(numel(ind123)./3);
    gData.corr_123_pos{vp} = corr_123_pos;
    
end

toc



%% average subjects not sessions

corr_all = gData.corr_all;
corr_123 = gData.corr_123;
corr_all_pos = cell2mat (gData.corr_all_pos); 
corr_all_pos_it  = corr_all_pos(1,:); corr_all_pos_it1 = reshape(corr_all_pos_it', 3, [])';
corr_all_pos_cat = corr_all_pos(2,:); corr_all_pos_cat1 = reshape(corr_all_pos_cat', 3, [])';
corr_all_pos1(:, 1:3) = corr_all_pos_it1; corr_all_pos1(:, 4:6) = corr_all_pos_cat1;
corr_123_pos = cell2mat (gData.corr_123_pos); 
corr_123_pos_it  = corr_123_pos(1,:); corr_123_pos_it1 = reshape(corr_123_pos_it', 3, [])';
corr_123_pos_cat = corr_123_pos(2,:); corr_123_pos_cat1 = reshape(corr_123_pos_cat', 3, [])';
corr_123_pos1(:, 1:3) = corr_123_pos_it1; corr_123_pos1(:, 4:6) = corr_123_pos_cat1;


region = 'all';
sub2exc = [];
corr_all = average_xGM (corr_all, region); 
corr_123 = average_xGM (corr_123, region); 
corr_all_pos = average_xGM (corr_all_pos1, region); 
corr_123_pos = average_xGM (corr_123_pos1, region); 

m_corr_all = mean(corr_all) ; std_corr_all = std(corr_all); se_corr_all = std_corr_all / sqrt ( size(corr_all, 1));
m_corr_123 = mean(corr_123) ; std_corr_123 = std(corr_123); se_corr_123 = std_corr_123 / sqrt ( size(corr_123, 1));
m_corr_all_pos = mean(corr_all_pos) ; std_corr_all_pos = std(corr_all_pos); se_corr_all_pos = std_corr_all_pos / sqrt ( size(corr_all_pos, 1));
m_corr_123_pos = mean(corr_123_pos) ; std_corr_123_pos = std(corr_123_pos); se_corr_123_pos = std_corr_123_pos / sqrt ( size(corr_123_pos, 1));

allPerf = [corr_123 corr_all]

%% 2Bar 
clear data
data.data = [allPerf(:, [1 3])] 

diff = data.data(:,1)- data.data(:,2);

figure(1); set(gcf,'Position', [0 0 550 650]); 
mean_S = mean(data.data, 1);
std_S = std(data.data, [], 1);
se_S = std_S / sqrt(size(data.data,1)); 
%hb = plot ([1 2], data.data, 'Color', [0 0 .1]); hold on;
hb = plot ([1 2], data.data); hold on;
%set(hb, 'lineWidth', 3, 'Marker', '.', 'MarkerSize',35, 'MarkerEdgeColor', [0 0 0]);hold on;
set(hb, 'lineWidth', 3, 'Marker', '.', 'MarkerSize',35);hold on;
h = bar (mean_S);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 3);
%set(hb,'linestyle','none', 'lineWidth', 3);
set(gca,'XTick',[1 2],'XTickLabel',{'', ''}, ...
    'FontSize', 30, 'linew',2, 'xlim', [0.25 2.75], 'ylim', [0 1.3] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 3);hold on;
%set(gca, 'LineWidth', 3);

%[p h stats] = signrank(data.data(:,1), data.data(:,2));
%disp (['W = ' num2str(stats.signedrank) '  ' ' p = ' num2str(p)]);

[h p ci t] = ttest(diff); % for the 123 trials
%[h p ci t] = ttest(diff, -0.0833); % for the 123 trials
%[h p ci t] = ttest(data.data(:,1), 1/6); % for the 123 trials
%[h p ci t] = ttest(data.data(:,2), 1/4); % for the 123 trials
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);


exportgraphics(gcf, '_21.png','Resolution', '150');

%% Simulation for chance levels in the 123 trials

clearvars -except allTrlInfo
trialinfo = allTrlInfo{4};
nPerm = 1000;

clear sim_res trRes clear allTRCNI

for permi = 1:nPerm
    countT = 0; 
    for triali = 1:size(trialinfo, 1)
        idCorr            = trialinfo(triali, 10);
        if idCorr ~= 4
            allItms           = [trialinfo(triali, 4:6) trialinfo(triali, 12:14)];
            allCats           = floor(allItms/100);
            corrItm           = trialinfo(triali, 3+idCorr);
            corrCat           = floor(corrItm/100);
            rand_it           = allItms(randsample(6, 1));
            rand_cat          = floor(rand_it/100);
            trResItm(triali,:)   = double(corrItm == rand_it);
            trResCat(triali,:)   = double(corrCat == rand_cat);
            countT = countT+1; 
        else
            trResItm(triali,:)    = nan;
            trResCat(triali,:)    = nan;
            trCatnotIt(triali, :) = nan;
        end
    end
    allTR_IT(permi,:)  = sum(trResItm, 'omitnan')/countT; 
    allTR_CAT(permi,:)  = sum(trResCat, 'omitnan')/countT; 
end

mean(allTR_IT)
mean(allTR_CAT)


%% Simulation for chance levels for ALL SSESSIONS

clearvars -except allTrlInfo
nPerm = 1000;

for sessi = 1:length(allTrlInfo)

    trialinfo = allTrlInfo{sessi};
    
    for permi = 1:nPerm
        countT = 0; 
        for triali = 1:size(trialinfo, 1)
            idCorr            = trialinfo(triali, 10);
            if idCorr ~= 4
                allItms           = [trialinfo(triali, 4:6) trialinfo(triali, 12:14)];
                allCats           = floor(allItms/100);
                corrItm           = trialinfo(triali, 3+idCorr);
                corrCat           = floor(corrItm/100);
                rand_it           = allItms(randsample(6, 1));
                rand_cat          = floor(rand_it/100);
                trResItm(triali,:)   = double(corrItm == rand_it);
                trResCat(triali,:)   = double(corrCat == rand_cat);
                countT = countT+1; 
            else
                trResItm(triali,:)    = nan;
                trResCat(triali,:)    = nan;
                trCatnotIt(triali, :) = nan;
            end
        end
        allTR_IT(permi,:)  = sum(trResItm, 'omitnan')/countT; 
        allTR_CAT(permi,:)  = sum(trResCat, 'omitnan')/countT; 
    end
    
    allSIT(sessi, :) = mean(allTR_IT);
    allSCAT(sessi, :) = mean(allTR_CAT);
end


%% 
mean (allSIT)
std(allSIT)
mean (allSCAT)
std(allSCAT)



%% Simulation for chance levels: all trials

clearvars -except allTrlInfo
trialinfo = allTrlInfo{1};
nPerm = 1000;

clear sim_res trResItm trResCat
for permi = 1:nPerm
    countT = 0; 
    for triali = 1:size(trialinfo, 1)
        idCorr                = trialinfo(triali, 10);
        if idCorr == 4
            allItms             = [trialinfo(triali, 4:6) trialinfo(triali, 12:14)];
            allCats             = floor(allItms/100);
            corr_seq_it         = allItms(1:3);
            corr_seq_cat        = floor(corr_seq_it/100);
            rand_seq_it         = allItms(randsample(6, 3));
            rand_seq_cat        = floor(rand_seq_it/100);
            trResItm(triali,:)  = sum(double(corr_seq_it == rand_seq_it));
            trResCat(triali,:)  = sum(double(corr_seq_cat == rand_seq_cat));
            countT = countT+1; 
        else
            trResItm(triali,:)  = nan;
            trResCat(triali,:)  = nan;
        end
    end
    allTR_IT(permi,:)  = sum(trResItm == 3)/countT; 
    allTR_CAT(permi,:)  = sum(trResCat == 3)/countT; 
end


mean(allTR_IT)
mean(allTR_CAT)



%% Simulation for chance levels: ALL SESSIONS

clearvars -except allTrlInfo
nPerm = 1000;

for sessi = 1:length(allTrlInfo)
    trialinfo = allTrlInfo{sessi};
    for permi = 1:nPerm
        countT = 0; 
        for triali = 1:size(trialinfo, 1)
            idCorr                = trialinfo(triali, 10);
            if idCorr == 4
                allItms             = [trialinfo(triali, 4:6) trialinfo(triali, 12:14)];
                allCats             = floor(allItms/100);
                corr_seq_it         = allItms(1:3);
                corr_seq_cat        = floor(corr_seq_it/100);
                rand_seq_it         = allItms(randsample(6, 3));
                rand_seq_cat        = floor(rand_seq_it/100);
                trResItm(triali,:)  = sum(double(corr_seq_it == rand_seq_it));
                trResCat(triali,:)  = sum(double(corr_seq_cat == rand_seq_cat));
                countT = countT+1; 
            else
                trResItm(triali,:)  = nan;
                trResCat(triali,:)  = nan;
            end
        end
        allTR_IT(permi,:)  = sum(trResItm == 3)/countT; 
        allTR_CAT(permi,:)  = sum(trResCat == 3)/countT; 
    end
    allSIT(sessi, :) = mean(allTR_IT);
    allSCAT(sessi, :) = mean(allTR_CAT);
end


%% 
mean (allSIT)
std(allSIT)
mean (allSCAT)
std(allSCAT)


%%
clear anCorr 
for permi = 1:nPerm
    nItCorr = sum(sim_res(permi, :,:, 1),3);
    an3Corr_it(permi,:) = sum(nItCorr == 3) /30 ;
    an2Corr_it(permi,:) = sum(nItCorr == 2) /30 ;
    an1Corr_it(permi,:) = sum(nItCorr == 1) /30 ;
    an0Corr_it(permi,:) = sum(nItCorr == 0) /30 ;
    
    nCatCorr = sum(sim_res(permi, :,:, 2),3);
    an3Corr_cat(permi,:) = sum(nCatCorr == 3) /30 ;
    an2Corr_cat(permi,:) = sum(nCatCorr == 2) /30 ;
    an1Corr_cat(permi,:) = sum(nCatCorr == 1) /30 ;
    an0Corr_cat(permi,:) = sum(nCatCorr == 0) /30 ;
    
end

mCorrItm = mean(an2Corr_it)
mCorrCat = mean(an2Corr_cat)


%%
mean(an3Corr_it)
(1/6) * (1/5) * (1/4)







%%


