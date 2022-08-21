% Load All Subjects  

clear
diary('_console_log.txt'); diary on; disp(string(datetime));
tic 
 
%global parameters % no ; for log
subjfir = 1;
subjend = 83;


%%power parameters
eLim            =       [-2 7]; %[-3 7];            %   max 5 
%eLim            =       [-6 6]; %[-3 7];            %   max 5 
xlimE           =       [-500 4000];       %   for cleaning at enc
xlimR           =       [-500 4000];      %   for cleaning at ret
cleaning        =       1;                 %      
takeMean        =       0;                 %   0 = median+nIQR ; 1 = mean+sd ; for cleaning
rereference     =       1;                 %
montage         =       'bipo';            %   'aver', 'bipo'
finalCut        =       [-.5 5];            %   in secs
timeRes         =       0.01; %0.01;               %    0.1 = 100ms; 0.01 = 10ms; 0.05 = 50ms or 'all' 
takeAllTrials   =       0; 
region          =       'hipp';
rawD            =       1; %power or raw amplitude time series
  
disp ([ 'eLim = ' num2str(eLim) newline ...
        'takeMean = ' num2str(takeMean) newline ...
        'rereference = ' num2str(rereference)  newline 'montage = ' montage newline ...
        'finalCut = ' num2str(finalCut) newline ...
        'timeRes = ' num2str(timeRes) newline ...
        'region = ' region newline ...
        ]);
 
%1) first load all events
% allEventInfo contains info for each trial in a table and all_events is in EEG.event format
[allEventInfo all_events] = loadLogsWM;
    
 
for subji = [1 8:9 25:26 29:39 51:54 59:60 63:64 65:83] %subjfir:subjend % subjs with hipp [1 8:9 25:26 29:39 51:54 59:60 63:64 65:83]
    disp (['Subj: ' num2str(subji)]);
    
    events = all_events{subji};
    EEG = load_data_WM (subji, events);
    
    if rereference
        EEG = re_reference_WM(EEG, montage); % chans2exc, EEG, average, bipolar    
    end
    
    EEG = select_electrodes_WM(EEG, region); %selects electrodes and assign labels to the subset
    
    if isempty(EEG.data)
       disp(['Subject ' num2str(subji) ' has no matching electrodes']);
       continue %skips the rest and begins the next iteration. Avoids subjects with no channels
    end
    
       
    if cleaning
        EEG.markers = zeros(size(EEG.data));
        for mi = 1:size(EEG.art_rejec, 1)
            time_per(1) = EEG.art_rejec(mi, 1);
            time_per(2) = EEG.art_rejec(mi, 2);
            EEG.markers(:, time_per(1):time_per(2)) = NaN;
        end 
    else
        EEG.markers = EEG.data;
    end
    
    [oneListTraces oneListIds oneListMarkers] = epoch_WM (EEG, eLim);
    
    %remove noisy trials 
    if cleaning
        [tr2exc_auto ] = remTriwithNans_WM (oneListTraces, oneListIds, oneListMarkers, xlimE, eLim);
                                            
               
        %first remove trials detected automatically
        if exist ('tr2exc_auto')
            disp(['oneListTraces before removing: ' num2str(size(oneListTraces))]);
            disp (['removing trials : ' num2str(tr2exc_auto) ]);
            oneListTraces_c  = oneListTraces(:,:,~tr2exc_auto); %I do this after decomposition to avoid NaNs
            oneListIds_c     = oneListIds(~tr2exc_auto);
            oneListMarkers_c = oneListMarkers(:,:,~tr2exc_auto);
 
            disp ([num2str(length(find(tr2exc_auto))) ' trials were excluded automatically']);
        end
        
        %exclude manually detected trilals
        t2excM = EEG.tr2exc_manu;
        [C ia ib] = intersect(t2excM, oneListIds_c);
        oneListTraces_c(:,:,ib) = []; 
        oneListMarkers_c(:,:,ib) = []; 
        oneListIds_c(ib) = [];
 
        disp ([num2str(length(find(ib))) ' trials were excluded manually']);
        disp(['oneListTraces after removing: ' num2str(size(oneListTraces_c))]);
    else
        oneListTraces_c  = oneListTraces; 
        oneListIds_c     = oneListIds;
        oneListMarkers_c = oneListMarkers; 
    end
    
    disp (['size onelistTraces_c > ' num2str(size(oneListTraces_c)) ]);
    
    if rawD
        cfg_contrasts.oneListIds_c        =       oneListIds_c; 
        cfg_contrasts.oneListTraces       =       oneListTraces_c; 
        cfg_contrasts.chanNames           =       struct2cell(EEG.chanlocs)'; %[{EEG.chanlocs.labels}']; 
        cfg_contrasts.subj                =       EEG.subj; 
        
    else
        [oneListPow] = extract_power_WM (oneListTraces_c, oneListIds_c, timeRes);
        cfg_contrasts.oneListIds_c        =       oneListIds_c; 
        cfg_contrasts.oneListPow          =       oneListPow; 
        cfg_contrasts.chanNames           =       struct2cell(EEG.chanlocs)'; %[{EEG.chanlocs.labels}']; 
        cfg_contrasts.subj                =       EEG.subj; 
    
    end
    
    filename = [EEG.subj '_out_contr'];
    varinfo=whos('cfg_contrasts');saveopt=''; if (varinfo.bytes >= 2^31) saveopt ='-v7.3'; else saveopt ='';end
    save(filename, 'cfg_contrasts', saveopt);
    
    
end


% % % create folder to store the subject data
fname = 'D:\_WM\analysis\out_contrasts\raw_traces\pfc\allTrials';
mkdir(fname);
sublist = dir('*contr.mat');
sublist = {sublist.name};
disp (['measurements -> ' num2str(length(sublist))]);

sessblo = cellfun(@(x) strsplit(string(x),'_'), sublist, 'UniformOutput', false);
sessblo = cellfun(@(x) x(1:3), sessblo, 'UniformOutput', false);
sessblo = cell2mat(cellfun(@(x) double(string(regexp(x, '\d+', 'match'))), sessblo, 'UniformOutput', false)');

[C, iaF, icF] = unique(sessblo(:,1),'first');   [C, iaL, icL] = unique(sessblo(:,1),'last');
for i = 1:length(C)
   dataF{i, 1} = C(i);dataF{i, 2} = iaF(i):iaL(i); %cell array that stores item id sorted 
end                                                   %in column 1 and instances in column 2 


for subji = 1:length(dataF) 
    
    d2m = sublist(dataF{subji, 2});
    s2u{subji,:} = d2m{1}(1:3);
    sb = sessblo(dataF{subji, 2},2:3);
    allC{subji,:} = merge_sessions_WM(d2m, sb);   %load files 1 by one
    
end

cd (fname)

for subji = 1:length(allC)

disp(['Subj: ' num2str(subji)]);
%filename = ['s' num2str(subji,'%02.f') '_out_contr']
filename = [ s2u{subji} '_out_contr'];
cfg_contrasts = allC{subji};
save(filename, 'cfg_contrasts', '-v7.3')

end

toc
 
cd .. 

