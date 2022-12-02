%% load cfg_contrasts
%% 

clearvars
region              = 'pfc';
paths = load_paths_WM(region);
filelistSess = getFiles(paths.out_contrasts);


frequncies2test = [{3:54} {3:8} {9:12} {13:29} {30:38} {39:54} ]';
fnames = {'3-150Hz' '3-8Hz' '9-12Hz' '13-29Hz' '30-75Hz' '75-150Hz' }'; fnames = fnames';


%frequncies2test = [{15:28}]';
%fnames = {'15-28Hz'}'; fnames = fnames';


win_width           = 5; 
mf                  = 1; 
meanInTime          = 0; 
avMeth              = 'pow';  % average across image repetitions or not
meanInFreq          = 0; 
takeElec            = 0; 
takeFreq            = 0;
TG                  = 1; %temporal generalization
contr2save          = {'ALL_EE'};
%contr2save          = {'SISC_EE' 'DISC_EE' 'DIDC_EE' 'SISC_EM2' 'DISC_EM2' 'DIDC_EM2' 'DISC_M2M2' 'DIDC_M2M2'}; %{};
%contr2save          = {'DISC_M2123V1' 'DIDC_M2123V1' 'DISC_M2123V2' 'DIDC_M2123V2' 'DISC_M2123CNCV1' ...
%                          'DIDC_M2123CNCV1' 'DISC_M2123CNCV2' 'DIDC_M2123CNCV2' 'DISC_M2123NC' 'DIDC_M2123NC' ...
%                          'DISC_M2A' 'DIDC_M2A' 'DISC_M2A123' 'DIDC_M2A123'}; 
bline               = [3 7];
acrossTrials        = 1;
batch_bin           = 1000;
%n2s                 = 1000000;
n2s                 = 10000;
loadSurr            = 0; 
zScType             = 'sess'; %'blo''sess' % 'allTrials' = all trials from all sessions and blocks

 
%diary([paths.results_path 'rsa_log.txt']); diary on; disp(string(datetime));
 
 
 
for sessi= 1:length(filelistSess) %this one starts at 1 and not at 3
    disp(['File > ' num2str(sessi)]);
    load([paths.out_contrasts filelistSess{sessi}]);   
    
    disp ([ 'fnames = ' fnames{:} newline ...
            'win_width = ' num2str(win_width) newline ...
            'mf = ' num2str(mf) newline ...
            'meanInTime = ' num2str(meanInTime) newline ...
            'meanInFreq = ' num2str(meanInFreq) newline ...
            'TG = ' num2str(TG) newline ...
            'bline = ' num2str(bline) newline ...
            'acrossTrials = ' num2str(acrossTrials) newline ...
            'batch_bin = ' num2str(batch_bin);
            ]);
        
    
    cfg_contrasts = normalize_WM(cfg_contrasts, acrossTrials, zScType, bline);
 
    cfg_contrasts.contr2save = contr2save';
    cfg_contrasts.n2s = n2s;
    cfg_contrasts.loadSurr = loadSurr;
    cfg_contrasts.batch_bin = batch_bin;
    

    if strcmp(avMeth,'pow')
        cfg_contrasts = average_repetitions(cfg_contrasts);
    end
    
    [out_contrasts] = create_contrasts_WM (cfg_contrasts);
        
        
    if ~exist('idxCH')
       idxCH = []; 
    end
    if ~exist('idxF')
       idxF = []; 
    end
  
    for freqi = 1:length(frequncies2test)
        %create folder
        fname = fnames{freqi};
        mkdir ([paths.results.bands fname]);
        cd ([paths.results.bands fname]);
        f           = frequncies2test{freqi}; 
 
        rsa_WM (out_contrasts, win_width, mf, f, meanInTime, meanInFreq, takeElec, takeFreq, idxCH, idxF, sessi, TG, 0)
 
        cd ..
    end
 
end

cd .. 

 
%%process Folders

clearvars -except region
paths = load_paths_WM(region); 
currentDir = pwd; 
mkdir(paths.results.bands)
cd (paths.results.bands)
fold = dir(); dirs = find(vertcat(fold.isdir));
fold = fold(dirs);
 
for foldi = 3:length(fold) % start at 3 cause 1 and 2 are . and ...
    
    clearvars -except contrasts fold foldi cmaps2use perms2use t2use cmapi region dupSym2use frequncies2test currentDir
    
    direct = fold(foldi);
    cd (direct.name)
 
    processFoldersWM; 

    cd .. 
end

cd (currentDir);
disp ('done')







 
 
%% Plot all trials
%%load cfg_contrasts
clear
tic
sublist = dir('*contr.mat');
sublist = {sublist.name};
disp (['measurements -> ' num2str(length(sublist))]);
 
for subji=1:length(sublist)
    disp(['File > ' num2str(subji)]);
    load(sublist{subji});
    s_all{subji} = cfg_contrasts; 
    %size(cfg_contrasts.oneListPow)
end
s_all = s_all';
toc
 
%% plot by trials
bline               = [3 7];
acrossTrials       = 1;
zScType = 'allTrials'; 
 
for subji = 1:length(sublist)
   cfg_contrasts = s_all{subji};    
   cfg_contrasts = normalize_WM(cfg_contrasts, acrossTrials, zScType, bline)
   oneListPow = cfg_contrasts.oneListPow;
   for triali = 1:size(oneListPow, 1)
       chanN = size(oneListPow, 2);
       pag2plot = ceil(chanN/100);
       cstri = 0;
       for pagi = 1:pag2plot
           figure(pagi);set(gcf,'Position', [0 0 2000 1500]); 
           csubi = 1;
           startChan = (cstri*100) +1
           endChan = (cstri*100) +100
           if endChan > chanN
               endChan = chanN
           end
           
           for chani = startChan:endChan
                subplot (10, 10, csubi)
                imagesc(squeeze(oneListPow(triali,chani,:,:))); 
                title(num2str(triali))
                set(gca,'YDir','normal')
                csubi = csubi + 1; 
           end
           
           cstri = cstri + 1;
           filename = [num2str(subji) '_tr_' num2str(triali) '_' num2str(pagi)];
           export_fig(pagi, filename,'-r200');    
           close all;  
       end
    end
   
end
 
 
 
 

%% plot by channels
bline               = [3 7];
acrossTrials       = 1;
zScType = 'allTrials'; 
 
for subji = 18:18 %1:length(sublist)
   cfg_contrasts = s_all{subji};    
   cfg_contrasts = normalize_WM(cfg_contrasts, acrossTrials, zScType, bline)
   oneListPow = cfg_contrasts.oneListPow;
   for chani = 1:size(oneListPow, 2)
       trialN = size(oneListPow, 1);
       pag2plot = ceil(trialN/100);
       cstri = 0;
       for pagi = 1:pag2plot
           figure(pagi);set(gcf,'Position', [0 0 2000 1500]); 
           csubi = 1;
           startTri = (cstri*100) +1
           endTri = (cstri*100) +100
           if endTri > trialN
               endTri = trialN
           end
           
           for triali = startTri:endTri
                subplot (10, 10, csubi)
                imagesc(squeeze(oneListPow(triali,chani,:,:))); 
                title(num2str(triali))
                set(gca,'YDir','normal')
                csubi = csubi + 1; 
           end
           
           cstri = cstri + 1;
           filename = [num2str(subji) '_' num2str(chani) '_' num2str(pagi)];
           export_fig(pagi, filename,'-r200');    
           close all;  
       end
    end
   
end
 
 
%% plot by channels ONE
bline               = [3 7];
takeAllTrials       = 1;
 
cfg_contrasts = normalize_WM(cfg_contrasts, takeAllTrials, 'blo', bline);
oneListPow = cfg_contrasts.oneListPow;
for chani = 1:size(oneListPow, 2)
   trialN = size(oneListPow, 1);
   pag2plot = ceil(trialN/100);
   cstri = 0;
   for pagi = 1:pag2plot
       figure(pagi);set(gcf,'Position', [0 0 2000 1500]); 
       csubi = 1;
       startTri = (cstri*100) +1
       endTri = (cstri*100) +100
       if endTri > trialN
           endTri = trialN
       end

       for triali = startTri:endTri
            subplot (10, 10, csubi)
            imagesc(squeeze(oneListPow(triali,chani,:,:))); 
            title(num2str(triali))
            set(gca,'YDir','normal')
            csubi = csubi + 1; 
       end

       cstri = cstri + 1;
       filename = [num2str(subji) '_' num2str(chani) '_' num2str(pagi)];
       export_fig(pagi, filename,'-r200');    
       close all;  
   end
    
 
 
end



%%
 
imagesc(squeeze(oneListPow(10,20,:,:))); colorbar;
 
 
 
 
 
 
 
 
%% count electrodes
%%load cfg_contrasts
clearvars
tic
sublist = dir('*contr.mat');
sublist = {sublist.name};
disp (['measurements -> ' num2str(length(sublist))]);
 
for subji=1:length(sublist)
    disp(['File > ' num2str(subji)]);
    load(sublist{subji});
    %s_all{subji} = size(cfg_contrasts.oneListPow, 2); 
    s_all{subji} = cfg_contrasts.chanNames; 
    %size(cfg_contrasts.oneListPow)
end
s_all = s_all';

elLength = cell2mat(cellfun(@length, s_all, 'un', 0))
sum(elLength)

toc

%% count hippocapmal electrodes

clear s_ids
for subji = 1:length(sublist)

    chanNames = allChans{subji};
    c = string(chanNames(:,5));
    ids = strfind(c, '38');
    id = ~cellfun('isempty',ids);
    if ~isempty(find(id))
        s_ids(subji,:) = 1; 
    else
        s_ids(subji,:) = 0;
    end
    
end

sum(s_ids)
 
 
%%


 
 
 
 
%%
