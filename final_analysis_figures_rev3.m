%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reviwer 3 comment 1 - SLIDE 22 - PFC ENC-MAINT WITHIN AND BETWEEN CORRELATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check behavioral data LME

























%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reviwer 3 comment 2 - SLIDE 20 - HIGHER ORDER REPRESENTATIONS IN PFC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check behavioral data LME
clear

[all_events allTrlInfo ] = loadLogsWM;
allTrlInfo = [allTrlInfo{:}];%% has to be done twice because of the nested 
allTrlInfo = [allTrlInfo{:}];%% has to be done twice because of the nested 
allTrlInfo = allTrlInfo(~cellfun('isempty',allTrlInfo));
allTrlInfo = allTrlInfo';




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
    %den2use = sum(~isnan(trialInfo(ind123,37:39))); 
    %corr_123_pos(1,:)=nansum(trialInfo(ind123,34:36))./den2use;
    %corr_123_pos(2,:)=nansum(trialInfo(ind123,37:39))./den2use;
    gData.corr_123_pos{vp} = corr_123_pos;
    
end

toc


%% average subjects not sessions

corr_all = gData.corr_all;
corr_123 = gData.corr_123;
corr_all_pos = cell2mat (gData.corr_all_pos); 
corr_all_pos_it  = corr_all_pos(1,:); corr_all_pos_it1 = reshape(corr_all_pos_it', 3, [])';
corr_all_pos_cat = corr_all_pos(2,:); corr_all_pos_cat1 = reshape(corr_all_pos_cat', 3, [])';
corr_all_pos1(:, 1:3) = corr_all_pos_it1; 
corr_all_pos1(:, 4:6) = corr_all_pos_cat1;
corr_123_pos = cell2mat (gData.corr_123_pos); 
corr_123_pos_it  = corr_123_pos(1,:); corr_123_pos_it1 = reshape(corr_123_pos_it', 3, [])';
corr_123_pos_cat = corr_123_pos(2,:); corr_123_pos_cat1 = reshape(corr_123_pos_cat', 3, [])';
corr_123_pos1(:, 1:3) = corr_123_pos_it1; 
corr_123_pos1(:, 4:6) = corr_123_pos_cat1;


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

allPerf = [corr_123 corr_all];
allPerfPos = [corr_123_pos corr_all_pos];



%% 6 Bar 
clear data
data.data = [corr_123_pos(:, 1:3) corr_all_pos(:, 1:3)]; 


figure(1); set(gcf,'Position', [0 0 550 650]); 
mean_S = mean(data.data);


h = bar (mean_S,'FaceColor','flat', 'linew', 2);hold on;
h.CData = [1 1 1; 1 1 1; 1 1 1; .5 .5 .5; .5 .5 .5; .5 .5 .5]; 
hb = plot ([1:6], data.data, 'k'); hold on; 
set(hb, 'lineWidth', 1, 'Marker', 'o', 'MarkerSize',4);hold on;
set(hb,'linestyle','none', 'linew', 2);
set(gca,'XTick',[1 6],'XTickLabel',{'', ''}, 'FontSize', 22, 'xlim', [0.25 6.75], 'ylim', [.3 1.3] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 1);hold on;
set(gca, 'LineWidth', 1);
box on

exportgraphics(gcf, 'myP.png','Resolution', '150');



%% LME
clc
allD = data.data(:); 
positions = [[ones(1, 32) ones(1, 32)*2 ones(1, 32)*3]' ; [ones(1, 32) ones(1, 32)*2 ones(1, 32)*3]'];
singORMult = [ones(1, 96) ones(1, 96)*2]';
subID = [1:32 1:32 1:32 1:32 1:32 1:32]';
d4LME = [allD positions singORMult subID];

tbl = table(d4LME(:,1), d4LME(:,2), d4LME(:,3), d4LME(:,4), ...
    'VariableNames',{'performance','position','trial_type', 'subID'});
lme = fitlme(tbl,'performance ~ position + trial_type + (1|subID)'); % random intercept model

lme

%% 
clc

[h p ci ts] = ttest(corr_123_pos(:, 1), corr_123_pos(:, 2));
disp (['t(' num2str(ts.df) ')= ' num2str(ts.tstat, 3) ',' ' p = ' num2str(p, 3)]);

[h p ci ts] = ttest(corr_123_pos(:, 1), corr_123_pos(:, 3));
disp (['t(' num2str(ts.df) ')= ' num2str(ts.tstat, 3) ',' ' p = ' num2str(p, 3)]);

[h p ci ts] = ttest(corr_123_pos(:, 2), corr_123_pos(:, 3));
disp (['t(' num2str(ts.df) ')= ' num2str(ts.tstat, 3) ',' ' p = ' num2str(p, 3)]);

[h p ci ts] = ttest(corr_all_pos(:, 1), corr_all_pos(:, 2));
disp (['t(' num2str(ts.df) ')= ' num2str(ts.tstat, 3) ',' ' p = ' num2str(p, 3)]);

[h p ci ts] = ttest(corr_all_pos(:, 1), corr_all_pos(:, 3));
disp (['t(' num2str(ts.df) ')= ' num2str(ts.tstat, 3) ',' ' p = ' num2str(p, 3)]);

[h p ci ts] = ttest(corr_all_pos(:, 2), corr_all_pos(:, 3));
disp (['t(' num2str(ts.df) ')= ' num2str(ts.tstat, 3) ',' ' p = ' num2str(p, 3)]);


[h p ci ts] = ttest(corr_123_pos(:, 1), mean([corr_123_pos(:, 2) corr_123_pos(:, 3)], 2) );
disp (['t(' num2str(ts.df) ')= ' num2str(ts.tstat, 3) ',' ' p = ' num2str(p, 3)]);

[h p ci ts] = ttest(corr_all_pos(:, 1), mean([corr_all_pos(:, 2) corr_all_pos(:, 3)], 2) );
disp (['t(' num2str(ts.df) ')= ' num2str(ts.tstat, 3) ',' ' p = ' num2str(p, 3)]);

singleItems = mean([corr_123_pos(:, 1) corr_123_pos(:, 2) corr_123_pos(:, 3)], 2) ; 
multiItems = mean([corr_all_pos(:, 1) corr_all_pos(:, 2) corr_all_pos(:, 3)], 2) ; 
[h p ci ts] = ttest(singleItems, multiItems);
disp (['t(' num2str(ts.df) ')= ' num2str(ts.tstat, 3) ',' ' p = ' num2str(p, 3)]);



%% ANOVA VERSION

d4ANOVA = [[1:32]' data.data(:, 1) data.data(:, 2) data.data(:, 3)]; 
% organize the data in a table
T = array2table(d4ANOVA(:,2:end));
T.Properties.VariableNames = {'Pos1' 'Pos2' 'Pos3'};
% create the within-subjects design
withinDesign = table([1 2 3]','VariableNames',{'Positions'});
withinDesign.Model = categorical(withinDesign.Positions);
% create the repeated measures model and do the anova
rm = fitrm(T,'Pos1-Pos3 ~ 1','WithinDesign',withinDesign);
AT = ranova(rm,'WithinModel','Model'); % remove comma to see ranova's table
%tbl = multcompare(rm, 'Model', 'ComparisonType', 'tukey-kramer'); % see: help RepeatedMeasuresModel/multcompare;  'tukey-kramer' (default), 'dunn-sidak', 'bonferroni','scheffe'
tbl = multcompare(rm, 'Model', 'ComparisonType', 'bonferroni'); % see: help RepeatedMeasuresModel/multcompare;  'tukey-kramer' (default), 'dunn-sidak', 'bonferroni','scheffe'


% output a conventional anova table
disp(anovaTable(AT, 'Measure (units)'));



%%





















%%

































%%