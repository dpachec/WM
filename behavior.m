%% check behavioral data
% first load events

clear

[allEventInfo all_events allTrlInfo subsessblo] = loadLogsWM;
allEventInfo = [allEventInfo{:}];
allEventInfo = [allEventInfo{:}];%% has to be done twice because of the nested 
allEventInfo = allEventInfo(~cellfun('isempty',allEventInfo));
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



%%
figure(); set(gcf, 'Position', [ 100 100 400 800])
subplot (3, 1, 1)
bar(m_corr_all); hold on; 
errorbar(1:2, m_corr_all, se_corr_all, 'k', 'LineWidth', 2);
title('Correct responses')
set(gca,'Xtick',1:4,'XtickLabel',{'all','123'}, 'FontSize', 12)
ylabel('Proportion correct')
set(gca, 'ylim', [0.1 .9])
subplot (3, 1, 2)
bar(m_corr_123); 
subplot (3, 1, 3)
bar(m_corr_all_pos); 


%% 
x1 = corr_123(:,1); x2 = ones(28, 1); 
y1 = corr_123(:,2); y2 = ones(28, 1)+1; 
slopes =(x2-x1)./(y2-y1)


slopeT = (1-(1/6)) / (2-(1/4))


[h p ci t] = ttest(slopes, slopeT)


%% 
x1 = corr_all(:,1); x2 = ones(28, 1); 
y1 = corr_all(:,2); y2 = ones(28, 1)+1; 
slopes =(x2-x1)./(y2-y1)


slopeT = (1-(1/6)) / (2-(1/4))


[h p ci t] = ttest(slopes, slopeT)
%% 
scatter([x2 y2], y1, y2)
hb = plot ([1 2], data.data); hold on;

%% 
diff12 = corr_123(:,2) - corr_123(:,1);



%%
mean(corr_all)
std(corr_all)


%% 2Bar 
clear data
data.data = [corr_all] 

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


%% 6Bar 
clear data
data.data = corr_123_pos; 

figure(1); set(gcf,'Position', [0 0 550 650]); 
mean_S = mean(data.data, 1);
std_S = std(data.data, [], 1);
se_S = std_S / sqrt(size(data.data,1)); 
%hb = plot ([1 2], data.data, 'Color', [0 0 .1]); hold on;
hb = plot ([1 2 3 4 5 6], data.data); hold on;
set(hb, 'lineWidth', 3, 'Marker', '.', 'MarkerSize',35);hold on;
h = bar (mean_S);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 3);
set(hb,'linestyle','none', 'lineWidth', 3);
set(gca,'XTick',[1 2 3 4 5 6],'XTickLabel',{'', ''}, ...
    'FontSize', 30, 'linew',2, 'xlim', [0.25 6.75], 'ylim', [0 1.3] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 3);hold on;
%set(gca, 'LineWidth', 3);

%[p h stats] = signrank(data.data(:,1), data.data(:,2));
%disp (['W = ' num2str(stats.signedrank) '  ' ' p = ' num2str(p)]);
[h p ci t] = ttest(data.data(:,1), data.data(:,2));
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);


export_fig(2, '_2.png','-transparent', '-r150');


%% 
[p,tbl,stats] = anova1(corr_123_pos(:,4:6));


%% 
clear S
Y = [corr_123(:) ; corr_all(:)];
S = [1:32' 1:32' 1:32' 1:32']';
F1 = [ones(32,1) ; zeros(32,1) ; ones(32,1) ; zeros(32,1)];
F2 = [ones(32,1) ; ones(32,1) ; zeros(32,1) ; zeros(32,1)];
FACTNAMES = {'Item-Category' '123-All'};

stats = rm_anova2(Y,S,F1,F2,FACTNAMES);



%% anova from matlab
clear S
Y = [corr_123 ; corr_all];
[~,~,stats] = anova2(Y,32);

c = multcompare(stats)
%% DECISION MAKERS : simulation for the 123 trials

clearvars -except allTrlInfo
trialinfo = allTrlInfo{9};
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
    allTR_IT(permi,:)  = sum(trResItm(~isnan(trResItm)))/countT; 
    allTR_CAT(permi,:)  = sum(trResCat(~isnan(trResCat)))/countT; 
end

mean(allTR_IT)
mean(allTR_CAT)



%% DECISION MAKERS : simulation for the all trials

clearvars -except allTrlInfo
trialinfo = allTrlInfo{20};
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














































%% backup

clear gData

for vp=1:length(allTrlInfo)
    
    trialInfo = allTrlInfo{vp};
    
    cues = trialInfo(:,10);
    
    indall = find(cues == 4);
    ind123 = find(cues ~= 4);
    ind1=find(cues ==1);
    ind2=find(cues ==2);
    ind3=find(cues ==3);

    
    % RTs
    RTs_all=zeros(size(trialInfo,1),3);
    RTs_all = trialInfo(:,42:44);
    rt123(ind123,:) = RTs_all(ind123,:);
    rtall(indall,:) = RTs_all(indall,:);


    % get behav data in nice way
    f=figure(1); 
   
 % plot bargraph completely correct ('token') trials vs type correct trials
    disp(['> indall = ' num2str(numel(indall)) '> ind123 = ' num2str(numel(ind123))]); 
    
    % all
    disp(['All item> ' num2str(sum(trialInfo(indall,40))) '/' num2str(numel(indall))]);
    disp(['All category > ' num2str(sum(trialInfo(indall,41))) '/' num2str(numel(indall))]);
    corr_all(1)=sum(trialInfo(indall,40))./numel(indall);
    corr_all(2)=sum(trialInfo(indall,41))./numel(indall);
    gData.corr_all_a(vp, :) = corr_all;
    
    %123
    disp(['123 item> ' num2str(sum(trialInfo(ind123,40))) '/' num2str(numel(ind123))]);
    disp(['123 category > ' num2str(sum(trialInfo(ind123,41))) '/' num2str(numel(ind123))]);
    corr_123(1)=sum(trialInfo(ind123,40))./numel(ind123);
    corr_123(2)=sum(trialInfo(ind123,41))./numel(ind123);
    gData.corr_123_a(vp, :) = corr_123;

   subplot(3,2,1)
    bar([corr_all;corr_123]);
    title('correct responses')
    set(gca,'Xtick',1:4,'XtickLabel',{'all','123'})
    legend('correct token','correct type')
    ylabel('relative correct reponses')
    set(gca, 'ylim', [0 1]);
    
    clear corr_all corr_123
    % plot correct type/token response for each pos*cue
    corr_all(1,:)=nansum(trialInfo(indall,34:36))./numel(indall);
    corr_all(2,:)=nansum(trialInfo(indall,37:39))./numel(indall);
    gData.corr_all_aPos{vp} = corr_all;
    %123
    corr_123(1,:)=nansum(trialInfo(ind123,34:36))./(numel(ind123)./3);
    corr_123(2,:)=nansum(trialInfo(ind123,37:39))./(numel(ind123)./3);
    gData.corr_123_aPos{vp} = corr_123;
    
    subplot(3,2,2)
    bar([corr_all;corr_123]);
    title('correct responses: across position')
    set(gca,'Xtick',1:4,'XtickLabel',{'all token','all type','123 token','123type'})
    legend('pos1','pos2','pos3')
    ylabel('relative correct reponses')
    clear corr_all corr_123

 % plot performance for each type
    type_enco=round(trialInfo(:,4:6)./100);
    type123=type_enco(ind123,:);
    typeall=type_enco(indall,:);

    token123_corr=type123((trialInfo(ind123,34:36))==1);
    token123_notcorr=type123((trialInfo(ind123,34:36))==0);
    tokenall_corr=typeall((trialInfo(indall,34:36))==1);
    tokenall_notcorr=typeall((trialInfo(indall,34:36))==0);
   
    
    for t=1:6
       token123_relcorr(t)= sum(token123_corr==t)./(sum(token123_corr==t)+sum(token123_notcorr==t));
       tokenall_relcorr(t)= sum(tokenall_corr==t)./(sum(tokenall_corr==t)+sum(tokenall_notcorr==t));
    end
    gData.token123_relcorr_a(vp, :) = token123_relcorr;
    gData.tokenall_relcorr_a(vp,:) = tokenall_relcorr;
    
    type123_corr=type123((trialInfo(ind123,37:39))==1);
    type123_notcorr=type123((trialInfo(ind123,37:39))==0);
    typeall_corr=typeall((trialInfo(indall,37:39))==1);
    typeall_notcorr=typeall((trialInfo(indall,37:39))==0);

    for t=1:6
       type123_relcorr(t)= sum(type123_corr==t)./(sum(type123_corr==t)+sum(type123_notcorr==t));
       typeall_relcorr(t)= sum(typeall_corr==t)./(sum(typeall_corr==t)+sum(typeall_notcorr==t));
    end
    gData.type123_relcorr_a(vp,:) = type123_relcorr;
    gData.typeall_relcorr_a(vp,:) = typeall_relcorr;
    
    subplot(3,2,3:4)
    bar([token123_relcorr;tokenall_relcorr;type123_relcorr;typeall_relcorr]);
    title('correct responses: across types')
    set(gca,'Xtick',1:4,'XtickLabel',{'123 token','all token','123 type ','all type'})
    legend('tree','robot','hand', 'house','planet','merkel')
 
    % plot performance across blocks
%     
%     for b=1:max(trialInfo(:,47))
%         tmp123=find(cues(1:length(rand_vec))<4 & (trialInfo(:,47)==b)');
%         tmpall=find(cues(1:length(rand_vec))==4 & (trialInfo(:,47)==b)');
%         
%         corr_all(b)=sum(trialInfo(tmpall,42))./numel(tmpall);
%         corr_123(b)=sum(trialInfo(tmp123,42))./numel(tmp123);
%         corr_typeall(b)=sum(trialInfo(tmpall,43))./numel(tmpall);
%         corr_type123(b)=sum(trialInfo(tmp123,43))./numel(tmp123);
%         
%         gData.corr_all_blocks(vp, :) = [corr_all zeros(1, 6-length(corr_all))];
%         gData.corr_123_blocks(vp, :) = [corr_123 zeros(1, 6-length(corr_123))];
%         gData.corr_typeall_blocks(vp, :) = [corr_typeall zeros(1, 6-length(corr_typeall))];
%         gData.corr_type123_blocks(vp, :) = [corr_type123 zeros(1, 6-length(corr_type123))];
%         
%     end
%     subplot(4,2,5:6)
%      bar([corr_123; corr_all;corr_type123; corr_typeall]);
%     title('correct responses: across blocks')
%     set(gca,'Xtick',1:4,'XtickLabel',{'123 token','all token','123 type ','all type'})
%     legend('block1','block2','block3','block4');

    
    % histograms RTs
    subplot(3,2,5)
    histogram(RTs_all(ind1,1),1000:1000:40000)
    hold on
    histogram(RTs_all(ind2,2),1000:1000:40000)
    hold on
    histogram(RTs_all(ind3,3),1000:1000:40000)
    legend('position1','position2','position3')
    xlabel('RT')
    ylabel('number of responses')
    title('RTs 123')
    
    subplot(3,2,6)
    histogram(RTs_all(indall,1),1000:1000:40000)
    hold on
    histogram(RTs_all(indall,2),1000:1000:40000)
    hold on
    histogram(RTs_all(indall,3),1000:1000:40000)
    legend('position1','position2','position3')
    xlabel('RT')
    ylabel('number of responses')
    title('RTs all')

%saveas(f,strcat(selvps{vp},'_behavsummary'))
%filename = strcat(selvps{vp},'_behavsummary');
%export_fig(1, filename, '-r300', '-transparent');

%close all

end




%% average subjects not sessions

corr_all_a = gData.corr_all_a;
corr_123_a = gData.corr_123_a;
gData.corr_all_aPosR = cell2mat (gData.corr_all_aPos); gData.corr_all_aPosR = reshape(gData.corr_all_aPosR, [], 6);
corr_all_aPosR = gData.corr_all_aPosR;
gData.corr_123_aPosR = cell2mat (gData.corr_123_aPos); gData.corr_123_aPosR = reshape(gData.corr_123_aPosR, [], 6);
corr_123_aPosR = gData.corr_123_aPosR;
type123_relcorr_a = gData.type123_relcorr_a;
typeall_relcorr_a = gData.typeall_relcorr_a;
token123_relcorr_a = gData.token123_relcorr_a;
tokenall_relcorr_a = gData.tokenall_relcorr_a;


corr_allA(1,:) = corr_all_a(1,:);
corr_allA(2,:) = mean(corr_all_a(2:3,:));
corr_allA(3,:) = mean(corr_all_a(4:5,:));
corr_allA(4,:) = mean(corr_all_a(6:7,:));
corr_allA(5,:) = mean(corr_all_a(8:9,:));
corr_allA(6,:) = mean(corr_all_a(10:14,:));
corr_allA(7,:) = mean(corr_all_a(15:18,:));
corr_allA(8,:) = mean(corr_all_a(19:20,:));
corr_allA(9,:) = mean(corr_all_a(21:22,:));
corr_allA(10,:) = mean(corr_all_a(23:24,:));
corr_allA(11,:) = mean(corr_all_a(25:26,:));
corr_allA(12,:) = mean(corr_all_a(27:28,:));
corr_allA(13,:) = mean(corr_all_a(29:30,:));
corr_allA(14,:) = mean(corr_all_a(31:32,:));
corr_allA(15,:) = mean(corr_all_a(33:37,:));
corr_allA(16,:) = mean(corr_all_a(38:39,:));
corr_allA(17,:) = mean(corr_all_a(40:41,:));
corr_allA(18,:) = mean(corr_all_a(42:43,:));
corr_allA(19,:) = mean(corr_all_a(44:45,:));
corr_allA(20,:) = mean(corr_all_a(46:47,:));
corr_allA(21,:) = mean(corr_all_a(48:49,:));
corr_allA(22,:) = mean(corr_all_a(50:51,:));
corr_allA(23,:) = mean(corr_all_a(52:53,:));
corr_allA(24,:) = mean(corr_all_a(54:55,:));
corr_allA(25,:) = mean(corr_all_a(56:57,:));
corr_allA(26,:) = mean(corr_all_a(58:61,:));
corr_allA(27,:) = mean(corr_all_a(62:64,:));
corr_allA(28,:) = mean(corr_all_a(65:68,:));
corr_allA(29,:) = mean(corr_all_a(69:72,:));
corr_allA(30,:) = mean(corr_all_a(73:76,:));


corr_123A(1,:) = corr_123_a(1,:);
corr_123A(2,:) = mean(corr_123_a(2:3,:));
corr_123A(3,:) = mean(corr_123_a(4:5,:));
corr_123A(4,:) = mean(corr_123_a(6:7,:));
corr_123A(5,:) = mean(corr_123_a(8:9,:));
corr_123A(6,:) = mean(corr_123_a(10:14,:));
corr_123A(7,:) = mean(corr_123_a(15:18,:));
corr_123A(8,:) = mean(corr_123_a(19:20,:));
corr_123A(9,:) = mean(corr_123_a(21:22,:));
corr_123A(10,:) = mean(corr_123_a(23:24,:));
corr_123A(11,:) = mean(corr_123_a(25:26,:));
corr_123A(12,:) = mean(corr_123_a(27:28,:));
corr_123A(13,:) = mean(corr_123_a(29:30,:));
corr_123A(14,:) = mean(corr_123_a(31:32,:));
corr_123A(15,:) = mean(corr_123_a(33:37,:));
corr_123A(16,:) = mean(corr_123_a(38:39,:));
corr_123A(17,:) = mean(corr_123_a(40:41,:));
corr_123A(18,:) = mean(corr_123_a(42:43,:));
corr_123A(19,:) = mean(corr_123_a(44:45,:));
corr_123A(20,:) = mean(corr_123_a(46:47,:));
corr_123A(21,:) = mean(corr_123_a(48:49,:));
corr_123A(22,:) = mean(corr_123_a(50:51,:));
corr_123A(23,:) = mean(corr_123_a(52:53,:));
corr_123A(24,:) = mean(corr_123_a(54:55,:));
corr_123A(25,:) = mean(corr_123_a(56:57,:));
corr_123A(26,:) = mean(corr_123_a(58:61,:));
corr_123A(27,:) = mean(corr_123_a(62:64,:));
corr_123A(28,:) = mean(corr_123_a(65:68,:));
corr_123A(29,:) = mean(corr_123_a(69:72,:));
corr_123A(30,:) = mean(corr_123_a(73:76,:));




corr_all_aPos(1,:) = corr_all_aPosR(1,:);
corr_all_aPos(2,:) = mean(corr_all_aPosR(2:3,:));
corr_all_aPos(3,:) = mean(corr_all_aPosR(4:5,:));
corr_all_aPos(4,:) = mean(corr_all_aPosR(6:7,:));
corr_all_aPos(5,:) = mean(corr_all_aPosR(8:9,:));
corr_all_aPos(6,:) = mean(corr_all_aPosR(10:14,:));
corr_all_aPos(7,:) = mean(corr_all_aPosR(15:18,:));
corr_all_aPos(8,:) = mean(corr_all_aPosR(19:20,:));
corr_all_aPos(9,:) = mean(corr_all_aPosR(21:22,:));
corr_all_aPos(10,:) = mean(corr_all_aPosR(23:24,:));
corr_all_aPos(11,:) = mean(corr_all_aPosR(25:26,:));
corr_all_aPos(12,:) = mean(corr_all_aPosR(27:28,:));
corr_all_aPos(13,:) = mean(corr_all_aPosR(29:30,:));
corr_all_aPos(14,:) = mean(corr_all_aPosR(31:32,:));
corr_all_aPos(15,:) = mean(corr_all_aPosR(33:37,:));
corr_all_aPos(16,:) = mean(corr_all_aPosR(38:39,:));
corr_all_aPos(17,:) = mean(corr_all_aPosR(40:41,:));
corr_all_aPos(18,:) = mean(corr_all_aPosR(42:43,:));
corr_all_aPos(19,:) = mean(corr_all_aPosR(44:45,:));
corr_all_aPos(20,:) = mean(corr_all_aPosR(46:47,:));
corr_all_aPos(21,:) = mean(corr_all_aPosR(48:49,:));
corr_all_aPos(22,:) = mean(corr_all_aPosR(50:51,:));
corr_all_aPos(23,:) = mean(corr_all_aPosR(52:53,:));
corr_all_aPos(24,:) = mean(corr_all_aPosR(54:55,:));
corr_all_aPos(25,:) = mean(corr_all_aPosR(56:57,:));
corr_all_aPos(26,:) = mean(corr_all_aPosR(58:61,:));
corr_all_aPos(27,:) = mean(corr_all_aPosR(62:64,:));
corr_all_aPos(28,:) = mean(corr_all_aPosR(65:68,:));
corr_all_aPos(29,:) = mean(corr_all_aPosR(69:72,:));
corr_all_aPos(30,:) = mean(corr_all_aPosR(73:76,:));


corr_123_aPos(1,:) = corr_123_aPosR(1,:);
corr_123_aPos(2,:) = mean(corr_123_aPosR(2:3,:));
corr_123_aPos(3,:) = mean(corr_123_aPosR(4:5,:));
corr_123_aPos(4,:) = mean(corr_123_aPosR(6:7,:));
corr_123_aPos(5,:) = mean(corr_123_aPosR(8:9,:));
corr_123_aPos(6,:) = mean(corr_123_aPosR(10:14,:));
corr_123_aPos(7,:) = mean(corr_123_aPosR(15:18,:));
corr_123_aPos(8,:) = mean(corr_123_aPosR(19:20,:));
corr_123_aPos(9,:) = mean(corr_123_aPosR(21:22,:));
corr_123_aPos(10,:) = mean(corr_123_aPosR(23:24,:));
corr_123_aPos(11,:) = mean(corr_123_aPosR(25:26,:));
corr_123_aPos(12,:) = mean(corr_123_aPosR(27:28,:));
corr_123_aPos(13,:) = mean(corr_123_aPosR(29:30,:));
corr_123_aPos(14,:) = mean(corr_123_aPosR(31:32,:));
corr_123_aPos(15,:) = mean(corr_123_aPosR(33:37,:));
corr_123_aPos(16,:) = mean(corr_123_aPosR(38:39,:));
corr_123_aPos(17,:) = mean(corr_123_aPosR(40:41,:));
corr_123_aPos(18,:) = mean(corr_123_aPosR(42:43,:));
corr_123_aPos(19,:) = mean(corr_123_aPosR(44:45,:));
corr_123_aPos(20,:) = mean(corr_123_aPosR(46:47,:));
corr_123_aPos(21,:) = mean(corr_123_aPosR(48:49,:));
corr_123_aPos(22,:) = mean(corr_123_aPosR(50:51,:));
corr_123_aPos(23,:) = mean(corr_123_aPosR(52:53,:));
corr_123_aPos(24,:) = mean(corr_123_aPosR(54:55,:));
corr_123_aPos(25,:) = mean(corr_123_aPosR(56:57,:));
corr_123_aPos(26,:) = mean(corr_123_aPosR(58:61,:));
corr_123_aPos(27,:) = mean(corr_123_aPosR(62:64,:));
corr_123_aPos(28,:) = mean(corr_123_aPosR(65:68,:));
corr_123_aPos(29,:) = mean(corr_123_aPosR(69:72,:));
corr_123_aPos(30,:) = mean(corr_123_aPosR(73:76,:));


type123_relcorrA(1,:) = type123_relcorr_a(1,:);
type123_relcorrA(2,:) = mean(type123_relcorr_a(2:3,:));
type123_relcorrA(3,:) = mean(type123_relcorr_a(4:5,:));
type123_relcorrA(4,:) = mean(type123_relcorr_a(6:7,:));
type123_relcorrA(5,:) = mean(type123_relcorr_a(8:9,:));
type123_relcorrA(6,:) = mean(type123_relcorr_a(10:14,:));
type123_relcorrA(7,:) = mean(type123_relcorr_a(15:18,:));
type123_relcorrA(8,:) = mean(type123_relcorr_a(19:20,:));
type123_relcorrA(9,:) = mean(type123_relcorr_a(21:22,:));
type123_relcorrA(10,:) = mean(type123_relcorr_a(23:24,:));
type123_relcorrA(11,:) = mean(type123_relcorr_a(25:26,:));
type123_relcorrA(12,:) = mean(type123_relcorr_a(27:28,:));
type123_relcorrA(13,:) = mean(type123_relcorr_a(29:30,:));
type123_relcorrA(14,:) = mean(type123_relcorr_a(31:32,:));
type123_relcorrA(15,:) = mean(type123_relcorr_a(33:37,:));
type123_relcorrA(16,:) = mean(type123_relcorr_a(38:39,:));
type123_relcorrA(17,:) = mean(type123_relcorr_a(40:41,:));
type123_relcorrA(18,:) = mean(type123_relcorr_a(42:43,:));
type123_relcorrA(19,:) = mean(type123_relcorr_a(44:45,:));
type123_relcorrA(20,:) = mean(type123_relcorr_a(46:47,:));
type123_relcorrA(21,:) = mean(type123_relcorr_a(48:49,:));
type123_relcorrA(22,:) = mean(type123_relcorr_a(50:51,:));
type123_relcorrA(23,:) = mean(type123_relcorr_a(52:53,:));
type123_relcorrA(24,:) = mean(type123_relcorr_a(54:55,:));
type123_relcorrA(25,:) = mean(type123_relcorr_a(56:57,:));
type123_relcorrA(26,:) = mean(type123_relcorr_a(58:61,:));
type123_relcorrA(27,:) = mean(type123_relcorr_a(62:64,:));
type123_relcorrA(28,:) = mean(type123_relcorr_a(65:68,:));
type123_relcorrA(29,:) = mean(type123_relcorr_a(69:72,:));
type123_relcorrA(30,:) = mean(type123_relcorr_a(73:76,:));


typeall_relcorrA(1,:) = typeall_relcorr_a(1,:);
typeall_relcorrA(2,:) = mean(typeall_relcorr_a(2:3,:));
typeall_relcorrA(3,:) = mean(typeall_relcorr_a(4:5,:));
typeall_relcorrA(4,:) = mean(typeall_relcorr_a(6:7,:));
typeall_relcorrA(5,:) = mean(typeall_relcorr_a(8:9,:));
typeall_relcorrA(6,:) = mean(typeall_relcorr_a(10:14,:));
typeall_relcorrA(7,:) = mean(typeall_relcorr_a(15:18,:));
typeall_relcorrA(8,:) = mean(typeall_relcorr_a(19:20,:));
typeall_relcorrA(9,:) = mean(typeall_relcorr_a(21:22,:));
typeall_relcorrA(10,:) = mean(typeall_relcorr_a(23:24,:));
typeall_relcorrA(11,:) = mean(typeall_relcorr_a(25:26,:));
typeall_relcorrA(12,:) = mean(typeall_relcorr_a(27:28,:));
typeall_relcorrA(13,:) = mean(typeall_relcorr_a(29:30,:));
typeall_relcorrA(14,:) = mean(typeall_relcorr_a(31:32,:));
typeall_relcorrA(15,:) = mean(typeall_relcorr_a(33:37,:));
typeall_relcorrA(16,:) = mean(typeall_relcorr_a(38:39,:));
typeall_relcorrA(17,:) = mean(typeall_relcorr_a(40:41,:));
typeall_relcorrA(18,:) = mean(typeall_relcorr_a(42:43,:));
typeall_relcorrA(19,:) = mean(typeall_relcorr_a(44:45,:));
typeall_relcorrA(20,:) = mean(typeall_relcorr_a(46:47,:));
typeall_relcorrA(21,:) = mean(typeall_relcorr_a(48:49,:));
typeall_relcorrA(22,:) = mean(typeall_relcorr_a(50:51,:));
typeall_relcorrA(23,:) = mean(typeall_relcorr_a(52:53,:));
typeall_relcorrA(24,:) = mean(typeall_relcorr_a(54:55,:));
typeall_relcorrA(25,:) = mean(typeall_relcorr_a(56:57,:));
typeall_relcorrA(26,:) = mean(typeall_relcorr_a(58:61,:));
typeall_relcorrA(27,:) = mean(typeall_relcorr_a(62:64,:));
typeall_relcorrA(28,:) = mean(typeall_relcorr_a(65:68,:));
typeall_relcorrA(29,:) = mean(typeall_relcorr_a(69:72,:));
typeall_relcorrA(30,:) = mean(typeall_relcorr_a(73:76,:));


token123_relcorrA(1,:) = token123_relcorr_a(1,:);
token123_relcorrA(2,:) = mean(token123_relcorr_a(2:3,:));
token123_relcorrA(3,:) = mean(token123_relcorr_a(4:5,:));
token123_relcorrA(4,:) = mean(token123_relcorr_a(6:7,:));
token123_relcorrA(5,:) = mean(token123_relcorr_a(8:9,:));
token123_relcorrA(6,:) = mean(token123_relcorr_a(10:14,:));
token123_relcorrA(7,:) = mean(token123_relcorr_a(15:18,:));
token123_relcorrA(8,:) = mean(token123_relcorr_a(19:20,:));
token123_relcorrA(9,:) = mean(token123_relcorr_a(21:22,:));
token123_relcorrA(10,:) = mean(token123_relcorr_a(23:24,:));
token123_relcorrA(11,:) = mean(token123_relcorr_a(25:26,:));
token123_relcorrA(12,:) = mean(token123_relcorr_a(27:28,:));
token123_relcorrA(13,:) = mean(token123_relcorr_a(29:30,:));
token123_relcorrA(14,:) = mean(token123_relcorr_a(31:32,:));
token123_relcorrA(15,:) = mean(token123_relcorr_a(33:37,:));
token123_relcorrA(16,:) = mean(token123_relcorr_a(38:39,:));
token123_relcorrA(17,:) = mean(token123_relcorr_a(40:41,:));
token123_relcorrA(18,:) = mean(token123_relcorr_a(42:43,:));
token123_relcorrA(19,:) = mean(token123_relcorr_a(44:45,:));
token123_relcorrA(20,:) = mean(token123_relcorr_a(46:47,:));
token123_relcorrA(21,:) = mean(token123_relcorr_a(48:49,:));
token123_relcorrA(22,:) = mean(token123_relcorr_a(50:51,:));
token123_relcorrA(23,:) = mean(token123_relcorr_a(52:53,:));
token123_relcorrA(24,:) = mean(token123_relcorr_a(54:55,:));
token123_relcorrA(25,:) = mean(token123_relcorr_a(56:57,:));
token123_relcorrA(26,:) = mean(token123_relcorr_a(58:61,:));
token123_relcorrA(27,:) = mean(token123_relcorr_a(62:64,:));
token123_relcorrA(28,:) = mean(token123_relcorr_a(65:68,:));
token123_relcorrA(29,:) = mean(token123_relcorr_a(69:72,:));
token123_relcorrA(30,:) = mean(token123_relcorr_a(73:76,:));


tokenall_relcorrA(1,:) = tokenall_relcorr_a(1,:);
tokenall_relcorrA(2,:) = mean(tokenall_relcorr_a(2:3,:));
tokenall_relcorrA(3,:) = mean(tokenall_relcorr_a(4:5,:));
tokenall_relcorrA(4,:) = mean(tokenall_relcorr_a(6:7,:));
tokenall_relcorrA(5,:) = mean(tokenall_relcorr_a(8:9,:));
tokenall_relcorrA(6,:) = mean(tokenall_relcorr_a(10:14,:));
tokenall_relcorrA(7,:) = mean(tokenall_relcorr_a(15:18,:));
tokenall_relcorrA(8,:) = mean(tokenall_relcorr_a(19:20,:));
tokenall_relcorrA(9,:) = mean(tokenall_relcorr_a(21:22,:));
tokenall_relcorrA(10,:) = mean(tokenall_relcorr_a(23:24,:));
tokenall_relcorrA(11,:) = mean(tokenall_relcorr_a(25:26,:));
tokenall_relcorrA(12,:) = mean(tokenall_relcorr_a(27:28,:));
tokenall_relcorrA(13,:) = mean(tokenall_relcorr_a(29:30,:));
tokenall_relcorrA(14,:) = mean(tokenall_relcorr_a(31:32,:));
tokenall_relcorrA(15,:) = mean(tokenall_relcorr_a(33:37,:));
tokenall_relcorrA(16,:) = mean(tokenall_relcorr_a(38:39,:));
tokenall_relcorrA(17,:) = mean(tokenall_relcorr_a(40:41,:));
tokenall_relcorrA(18,:) = mean(tokenall_relcorr_a(42:43,:));
tokenall_relcorrA(19,:) = mean(tokenall_relcorr_a(44:45,:));
tokenall_relcorrA(20,:) = mean(tokenall_relcorr_a(46:47,:));
tokenall_relcorrA(21,:) = mean(tokenall_relcorr_a(48:49,:));
tokenall_relcorrA(22,:) = mean(tokenall_relcorr_a(50:51,:));
tokenall_relcorrA(23,:) = mean(tokenall_relcorr_a(52:53,:));
tokenall_relcorrA(24,:) = mean(tokenall_relcorr_a(54:55,:));
tokenall_relcorrA(25,:) = mean(tokenall_relcorr_a(56:57,:));
tokenall_relcorrA(26,:) = mean(tokenall_relcorr_a(58:61,:));
tokenall_relcorrA(27,:) = mean(tokenall_relcorr_a(62:64,:));
tokenall_relcorrA(28,:) = mean(tokenall_relcorr_a(65:68,:));
tokenall_relcorrA(29,:) = mean(tokenall_relcorr_a(69:72,:));
tokenall_relcorrA(30,:) = mean(tokenall_relcorr_a(73:76,:));


m_corr_all_a = mean(corr_allA);std_corr_all_a = std(corr_allA); se_corr_all_a = std_corr_all_a / sqrt(length(corr_allA)); 
m_corr_123_a = mean(corr_123A);std_corr_123_a = std(corr_123A); se_corr_123_a = std_corr_123_a / sqrt(length(corr_123A)); 
m_corr_all_aPos = mean(corr_all_aPos);std_corr_all_aPos = std(corr_all_aPos); se_corr_all_aPos = std_corr_all_aPos / sqrt(length(corr_all_aPos)); 
m_corr_123_aPos = mean(corr_123_aPos);std_corr_123_aPos = std(corr_123_aPos); se_corr_123_aPos = std_corr_123_aPos / sqrt(length(corr_123_aPos)); 
m_type123_relcorrA = mean(type123_relcorrA); std_type123_relcorrA = std(type123_relcorrA);se_type123_relcorrA = std_type123_relcorrA  / sqrt(length(type123_relcorrA)); 
m_typeall_relcorrA = mean(typeall_relcorrA); std_typeall_relcorrA = std(typeall_relcorrA);se_typeall_relcorrA = std_typeall_relcorrA  / sqrt(length(typeall_relcorrA)); 
m_token123_relcorrA = mean(token123_relcorrA); std_token123_relcorrA = std(token123_relcorrA);se_token123_relcorrA = std_token123_relcorrA  / sqrt(length(token123_relcorrA)); 
m_tokenall_relcorrA = mean(tokenall_relcorrA); std_tokenall_relcorrA = std(tokenall_relcorrA);se_tokenall_relcorrA = std_tokenall_relcorrA  / sqrt(length(tokenall_relcorrA)); 


%% stats 

[h p ci t ] = ttest(corr_allA(:,1), corr_allA(:,2));

%% one way anova for the items
%[p,tbl,stats] = anova1(corr_123_aPos);
[p,tbl,stats] = anova1(corr_123_aPos(4:6));

multcompare(stats)
%% non-parametric
p = kruskalwallis(corr_all_aPos)

%% plot 

ftsze = 14;
f = figure(1); set(f, 'Position', [0 0 1200 800])
subplot(3,2,1)
x = [m_corr_all_a;m_corr_123_a];
err = [se_corr_all_a; se_corr_all_a];
bar(x);hold on;
myErrorBar(x, err);
title('Correct responses')
set(gca,'Xtick',1:4,'XtickLabel',{'all','123'}, 'FontSize', ftsze)
legend('correct token','correct type', 'Location', 'southeast')
ylabel('Proportion correct')

subplot(3,2,2)
y = [m_corr_all_aPos(1:3);m_corr_all_aPos(4:6);m_corr_123_aPos(1:3);m_corr_123_aPos(4:6)];
err = [se_corr_all_aPos(1:3); se_corr_all_aPos(4:6); se_corr_all_aPos(1:3); se_corr_all_aPos(4:6)];
bar(y);hold on;
myErrorBar(y, err);
title('Performance across position')
set(gca,'Xtick',1:4,'XtickLabel',{'all token','all type','123 token','123type'}, 'FontSize', ftsze)
legend('pos1','pos2','pos3', 'Location', 'southeast')
ylabel('Proportion correct')

subplot(3,2,3:4)
y = [m_token123_relcorrA; m_tokenall_relcorrA; m_type123_relcorrA; m_typeall_relcorrA];
err = [se_token123_relcorrA;se_tokenall_relcorrA;se_type123_relcorrA;se_typeall_relcorrA];
bar(y); hold on
myErrorBar(y, err);
title('correct responses: across types')
set(gca,'Xtick',1:4,'XtickLabel',{'123 token','all token','123 type ','all type'})
legend('tree','robot','hand', 'house','planet','merkel', 'Location', 'southeast')



%% plot group data

m_corr_all_a = mean(gData.corr_all_a); std_corr_all_a = std(gData.corr_all_a); se_corr_all_a = std_corr_all_a / sqrt(19); 
m_corr_123_a = mean(gData.corr_123_a); std_corr_123_a = std(gData.corr_123_a); se_corr_123_a = std_corr_123_a / sqrt(19); 

gData.corr_all_aPosR = cell2mat (gData.corr_all_aPos); gData.corr_all_aPosR = reshape(gData.corr_all_aPosR, [], 6);
m_corr_all_aPos = mean (gData.corr_all_aPosR); std_corr_all_aPos = std(gData.corr_all_aPosR); se_corr_all_aPos = std_corr_all_aPos / sqrt(19); 
gData.corr_123_aPosR = cell2mat (gData.corr_123_aPos); gData.corr_123_aPosR = reshape(gData.corr_123_aPosR, [], 6);
m_corr_123_aPos = mean (gData.corr_123_aPosR); std_corr_123_aPos = std(gData.corr_123_aPosR); se_corr_123_aPos = std_corr_123_aPos / sqrt(19); 

% plot performance for each type
m_token123_relcorr_a = mean(gData.token123_relcorr_a); std_token123_relcorr_a = std(gData.token123_relcorr_a); se_token123_relcorr_a = std_token123_relcorr_a / sqrt(19); 
m_tokenall_relcorr_a = mean(gData.tokenall_relcorr_a); std_tokenall_relcorr_a = std(gData.tokenall_relcorr_a); se_tokenall_relcorr_a = std_tokenall_relcorr_a / sqrt(19); 

m_type123_relcorr_a = mean(gData.type123_relcorr_a); std_type123_relcorr_a = std(gData.type123_relcorr_a); se_type123_relcorr_a = std_type123_relcorr_a / sqrt(19); 
m_typeall_relcorr_a = mean(gData.typeall_relcorr_a); std_typeall_relcorr_a = std(gData.typeall_relcorr_a); se_typeall_relcorr_a = std_typeall_relcorr_a / sqrt(19); 

% performance by blocks
% gData.corr_all_blocks(gData.corr_all_blocks == 0) = nan; m_corr_all_blocks = mean(gData.corr_all_blocks, 'omitnan'); 
% gData.corr_123_blocks(gData.corr_123_blocks == 0) = nan; m_corr_123_blocks = mean(gData.corr_123_blocks, 'omitnan');
% nsu = sum(~isnan(gData.corr_all_blocks));nsu1 = sum(~isnan(gData.corr_123_blocks));
% std_corr_all_blocks = std(gData.corr_all_blocks, 'omitnan'); se_corr_all_blocks = std_corr_all_blocks ./ sqrt(nsu); 
% std_corr_123_blocks = std(gData.corr_123_blocks, 'omitnan'); se_corr_123_blocks = std_corr_123_blocks ./ sqrt(nsu1); 

ftsze = 14;
f = figure(1); set(f, 'Position', [0 0 800 800])
subplot(3,2,1)
y = [m_corr_all_a;m_corr_123_a];
err = [se_corr_all_a; se_corr_all_a];
bar(y);hold on;
myErrorBar(y, err);
title('Correct responses')
set(gca,'Xtick',1:4,'XtickLabel',{'all','123'}, 'FontSize', ftsze)
legend('correct token','correct type', 'Location', 'southeast')
ylabel('Proportion correct')

subplot(3,2,2)
y = [m_corr_all_aPos(1:3);m_corr_all_aPos(4:6);m_corr_123_aPos(1:3);m_corr_123_aPos(4:6)];
err = [std_corr_all_aPos(1:3); std_corr_all_aPos(4:6); std_corr_all_aPos(1:3); std_corr_all_aPos(4:6)];
bar(y);hold on;
myErrorBar(y, err);
title('Performance across position')
set(gca,'Xtick',1:4,'XtickLabel',{'all token','all type','123 token','123type'}, 'FontSize', ftsze)
legend('pos1','pos2','pos3', 'Location', 'southeast')
ylabel('Proportion correct')
clear corr_all corr_123

subplot(3,2,3:4)
y = [m_corr_all_aPos;m_corr_123_aPos];
err = [std_corr_all_aPos; std_corr_all_aPos];
bar(y);hold on;
myErrorBar(y, err);
title('correct responses: across types')
set(gca,'Xtick',1:4,'XtickLabel',{'123 token','all token','123 type ','all type'}, 'FontSize', ftsze)
legend('tree','robot','hand', 'house','planet','merkel', 'Location', 'southeast')


subplot(3,2,3:4)
bar([token123_relcorr;tokenall_relcorr;type123_relcorr;typeall_relcorr]);
title('correct responses: across types')
set(gca,'Xtick',1:4,'XtickLabel',{'123 token','all token','123 type ','all type'})
legend('tree','robot','hand', 'house','planet','merkel')


% subplot(3,2,5:6)
% y = [m_corr_all_blocks;m_corr_123_blocks];
% err = [std_corr_all_blocks; std_corr_123_blocks];
% bar(y);hold on;
% myErrorBar(y, err);
% title('correct responses: across blocks')
% set(gca,'Xtick',1:4,'XtickLabel',{'123 token','all token','123 type ','all type'}, 'FontSize', ftsze)
% legend('block1','block2','block3','block4');




%% decision makers 
% % % for the 123 trials chance is easy to determine
trialinfo = allTrlInfo{3};


clear nLures chanceall_cat chance123_cat
for triali = 1:size(trialinfo, 1)
    
    allItms = [trialinfo(triali, 4:6) trialinfo(triali, 12:14)];
    idCorr = trialinfo(triali, 10);
    allCats = floor(allItms/100);
    nLures = length(unique(allCats));
       
    if idCorr ~= 4
       corrItm = trialinfo(triali, 3+idCorr);
       %catUs = floor(corrItm/100)
       chance123_cat = 1/nLures; 
       chanceall_cat = 1/6; 
    end
        
end


%% 







%%


