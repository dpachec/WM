%%
clear, clc
%addpath('/Users/danielpacheco/Documents/LFP_ANALYSIS/_WM/data')
%all_pat_folder='E:\WM_paradigm\iEEG\WM_Freiburg\data\iEEG\';
if ~ismac
    all_pat_folder='D:\_WM\data\iEEG\';
else 
    all_pat_folder='/Users/danielpacheco/Documents/LFP_ANALYSIS/_WM/data/iEEG/';
end
cd(all_pat_folder);
%all_pat=ls('sub*');
all_pat = dir('sub*');all_pat = {all_pat.name};

%% read in logfiles and create trialinfo for each session

tic
for sub=19:19 %size(all_pat,2)
sel_sub=char(all_pat(sub));
fname1 = string(strcat(all_pat_folder,sel_sub));
cd(fname1);
all_sess=dir('Sess*');all_sess = {all_sess.name};
    for ses=1:size(all_sess,2)
    	sel_sess=all_sess(ses);
        fname = string(strcat(all_pat_folder,sel_sub,'/',sel_sess,'/Beh'));
        cd(fname);
        all_logs = dir('*log');all_logs = {all_logs.name};
        for log=1:size(all_logs, 2)        
            sel_log=char(all_logs(log)); 
            sel_log(isspace(sel_log))=[]; % delete blanks for better indei
            [trlinfo,datelog]=log2trlinfo(sel_log);
            fname2 = strcat(sel_sub(1:5),'_',sel_log(1:end-4), '_trialinfo');
            save(fname2,'trlinfo','datelog');
        end
    end
end
cd(all_pat_folder);
toc
disp('loaded sessions')

%% read in data
% check for multiple sess
% loop sub, sess, logfiles
tic

clear allEventInfo allD
for sub=18:18%size(all_pat,2)
sel_sub=all_pat(sub);
fname = char(strcat(all_pat_folder,sel_sub));
cd(fname);
all_sess=ls('Sess*');
all_sess=dir('Sess*');all_sess = {all_sess.name};
    for ses=2:2 %size(all_sess, 2)
        sel_sess=all_sess(ses);

        % load data
        folder_str = string(strcat(all_pat_folder,sel_sub,'/',sel_sess,'/Macro/'));
        data     = data_readin_lukas(folder_str);

        sr_log=10000;
        sr_eeg=data.fsample;
        %sr log filehigher:upsample eeg trigger
        trig_channel_eeg=data.trial{1,1}(strcmp('Trigger',data.label),:);
        %     %checkphotodiode values
        %    all_trig_value=unique(trig_channel_eeg);
        %    %check for trigger def
        %    for t=1:numel(all_trig_value)
        %    trig_channel_eeg_tmp=trig_channel_eeg==(all_trig_value(t));
        %    % plotfor sanity check
        %    figure;plot(trig_channel_eeg_tmp)
        %   trig_val_ind(t)= input('trigger def photodiode ok?')
        %    close all
        %    end
        % trig_channel_eeg=trig_channel_eeg==all_trig_value(find(str2double(trig_val_ind)));
        % for now consider all >0 trigger
        trig_channel_eeg=trig_channel_eeg>0;
        trig_channel_eeg=[0, trig_channel_eeg];
        trig_channel_eeg=diff(trig_channel_eeg);% get all offsets
        trig_channel_eeg=double(trig_channel_eeg==1);
    
        trig_sp_eeg=find(trig_channel_eeg);
        diff_trig_eeg=diff(trig_sp_eeg).*(sr_log/sr_eeg);
        %use crosscorrelation to find match(problem:clocks drift
        %apart,matchingdifference might be more stable)
        trig_channel_eeg_rs=resample(trig_channel_eeg,sr_log,sr_eeg);
    
        fname = string(strcat(all_pat_folder,sel_sub,'/',sel_sess,'/Beh'));
        cd(fname);
        trl_files=dir('*trialinfo.mat'); trl_files = {trl_files.name};
    
    
        % load trialinfo files
        for trl=1:size(trl_files,2)
            sel_log=char(trl_files(trl));
            load(sel_log);
            % construct trigger channels from trlinfo
            % 1: block_type 2=WM
            % 2: WM block
            % 3: trialnumber
            % 4:9 enco item & times
            % 10:11cue & times
            % 12:15: test item token + time ( at test always 6 items are presented, the three enco items plus three added items (12:14)
            % 16:21 which response codes for which item (randvec: position of x in randvec codes for item: ie randvec(4)=2 means response 4 means enco item 2, also item shown at position 4 on test screen)
            % 22-27 recode randvec (column 22: correct response for item 1, 25: correct response for item 2 etc...)
            % 28-30: responses
            % 31-33: lure pos (there are always some some similar lures in the response display, this codes on which position similar lures to the respective enco item 1-3 are shown, should make sense when you check the token vs type correct)
            % 34-36. all correct (token correct)
            % 37:39: correct on type level
            % 40: completely correct all token
            % 41: completely correct all types
            % 42-44: RT for response relative to test onset
            if ~isempty(strfind(sel_log,'wmblock'));
                trlinfo = getPresentedItems(trlinfo);
                log_type='wmblock';
                trig_log=[trlinfo(:,7);trlinfo(:,8);trlinfo(:,9);trlinfo(:,11)];
                trig_log=sort(trig_log);
                
                clear c1 c2 c3 c4 c5 c6 c7 c8 c91011 event_info
                c1 = [trlinfo(:,7);trlinfo(:,8);trlinfo(:,9);trlinfo(:,11); ...
                                   trlinfo(:,15); ... 
                                   trlinfo(:,15) + trlinfo(:,42); ...
                                   trlinfo(:,15) + trlinfo(:,43); ...
                                   trlinfo(:,15) + trlinfo(:,44)];
               %c2 = trial type
               c2 = 1:8;
               c2 = repmat(c2,1,60)';   
               c3 = trlinfo(:, 10); 
               c3 = repelem(c3,8,1);

               
               c4 = [trlinfo(:,4);trlinfo(:,5);trlinfo(:,6);NaN(60,1); ...
                                   NaN(60,1); ... 
                                   trlinfo(:,45); ...
                                   trlinfo(:,46); ...
                                   trlinfo(:,47)];
               c5 = [NaN(60,1);NaN(60,1);NaN(60,1);NaN(60,1); ...
                                   NaN(60,1); ... 
                                   trlinfo(:,28); ...
                                   trlinfo(:,29); ...
                                   trlinfo(:,30)];

               c6 = [ trlinfo(:,22);trlinfo(:,23);trlinfo(:,24);NaN(60,1); ...
                                   NaN(60,1); ... 
                                   trlinfo(:,22); ...
                                   trlinfo(:,23); ...
                                   trlinfo(:,24)];
                c7 = [NaN(60,1);NaN(60,1);NaN(60,1);NaN(60,1); ...
                                   NaN(60,1); ... 
                                   trlinfo(:,34); ...
                                   trlinfo(:,35); ...
                                   trlinfo(:,36)];
                c8 = [NaN(60,1);NaN(60,1);NaN(60,1);NaN(60,1); ...
                                   NaN(60,1); ... 
                                   trlinfo(:,37); ...
                                   trlinfo(:,38); ...
                                   trlinfo(:,39)];
                               
                c9 = trlinfo(:, 40);
                c9 = repelem(c9,8,1);
                c10 = trlinfo(:, 41); 
                c10 = repelem(c10,8,1);
                
                c11 = trlinfo(:, 1); 
                c11 = repelem(c11,8,1);
                c12 = trlinfo(:, 2); 
                c12 = repelem(c12,8,1);
                c13 = trlinfo(:, 3); 
                c13 = repelem(c13,8,1);
                 
                c141516 = trlinfo(:, 4:6);
                c141516 = repelem(c141516,8,1);
                c171819 = trlinfo(:, 12:14);
                c171819 = repelem(c171819,8,1); 
                
               c20 = [NaN(60,1);NaN(60,1);NaN(60,1);NaN(60,1); ...
                   NaN(60,1); ... 
                   trlinfo(:,42); ...
                   trlinfo(:,43); ...
                   trlinfo(:,44)];
                
                
                [c1 idx]=sort(c1);
                c4 = c4(idx);
                c5 = c5(idx);
                c6 = c6(idx);
                c7 = c7(idx);
                c8 = c8(idx);
                c20 = c20(idx);
                
                event_info = [c1 c2 c3 c4 c5 c6 c7 c8 c9 c10 c11 c12 c13 c141516 c171819 c20];
                
                [uni idx] = unique(c1);
                event_info = event_info(idx, :);
                
                %%event by type
                %3 = cue; 
                %4= id(1-3 = of pres, 6-8 of respnse)
                %5 = 6-8position response at test
                %6 = position of item at test (repeated in 1-3 just in case)
                %7 = correct / incorrect item
                %8 = correct / incorrect type
                %9 = all correct type
                %10 = all correct token
                %11 = block type (2 = WM)
                %12 = block number
                %13 = trial number
                %14-19 = all presented items (14-16 = presented, 16-19 new)
                %20 = latency of response with respect to test onset
                
                

    
         allTrigLog{trl} = trig_log;
         allEventInfo{sub,:}{ses,:}{trl} = event_info;
      
        diff_trig_log=diff(trig_log);


  
        trig_onoffset=trig_log(1);
        %construct log trigger channel 
        trig_channel_log=zeros(1,max(trig_log)+trig_log(1));
        trig_channel_log(trig_log)=1;
        [x,lags]= xcorr(trig_channel_eeg_rs,trig_channel_log);
        [maxr,ind]=max(x);
   
        max_lag=lags(ind);  
    
        %evaluate matchvisually
%         figure();
%         plot(trig_channel_eeg_rs(max_lag:max_lag+(numel(trig_channel_log))));
%         hold on
%         plot(trig_channel_log);
%         [check_visual]= 1; %input('match makes sense?,0 if not!')
%         check_visual=str2double(check_visual);
%     
%         %based max_lag match triggers
%         start_logindata=round((max_lag+trig_log(1))/(sr_log/sr_eeg));% downsample max lag (onset log_trig in data)
%         end_logindata=start_logindata+round((trig_log(end)-trig_log(1))/(sr_log/sr_eeg));
%         cut_trig_eeg=trig_channel_eeg(start_logindata:end_logindata);
% 
%         trig_sp_cut=find(cut_trig_eeg);
%         check_num=numel(trig_sp_cut)==numel(trig_log);
% 
%         % check fit by comparing trigger differences
%         eeg_diff= diff(trig_sp_cut);
%         log_diff=round(diff(trig_log)/(sr_log/sr_eeg));
%         %log_diff(end+1) = 0;
% 
%         fit_tolerance=20;% 10ms jitter tolerance
        %check_diff=any((eeg_diff'-log_diff)>=fit_tolerance);

%         if check_num==0|check_diff==1|check_visual==0
%         % write file: trigger sp with trigger definition
%         error(strcat('no trigger match',sel_log)) 
%         end

        %%%  for sub4extra fix needed (lose photodiode)
        
                        
            elseif ~isempty(strfind(sel_log,'practice'));    
%                 log_type='practice';
%                 trig_log=trlinfo(:,4);
%                 event_info = [];
            elseif ~isempty(strfind(sel_log,'decoding'));
%                 log_type='decoding';
%                 trig_log=trlinfo(:,8);
%                 event_info =[];
            else
            end
        end
        
        
    end
    allD{sub,:}{ses,:} = data;
end

disp('all_data_imported');



toc


%% create EEG
clearvars -except data allEventInfo allD

clear EEG
%data = allD{subji}{ses};
EEG.data = data.trial{1};
EEG.srate = data.fsample;
EEG.chanlocs = struct('labels', data.label);
EEG.trials = 1;
EEG.setname = 'Freiburg';
EEG.icawinv = [];
EEG.icaweights = [];
EEG.icasphere = [];
EEG.icaact = [];
EEG.nbchan = size(EEG.data,1);
EEG.pnts = size(EEG.data,2);
EEG.xmax = max(EEG.data(:));
EEG.xmin = min(EEG.data(:));


%% plot to select experiment data
%EEG = ALLEEG{2};
eventChannel = 'Trigger';TTL = strmatch(eventChannel, {EEG.chanlocs.labels}, 'exact');chanids = [TTL];
eegplot(EEG.data(chanids,:), 'srate', EEG.srate, 'eloc_file',EEG.chanlocs(chanids), ...
    'winlength', 50, 'spacing', 100000000);


%% 1) FIRST REMOVE EMPTY CHANNELS

EEG.data(57:58,:) = [];
EEG.chanlocs(57:58,:) = [];
EEG.nbchan = size(EEG.data,1);

%% 2) and non experimental data

EEG.data = EEG.data (:,2557*EEG.srate:3584*EEG.srate);
EEG.pnts = length(EEG.data);

disp('data extracted');





%% Extract triggers from TTL channel
eventChannel = 'Trigger';TTL = strmatch(eventChannel, {EEG.chanlocs.labels}, 'exact');
strEvent = 'X > 900';

EEG.event = [];
EEG = pop_chanevent (EEG, TTL, 'edge', 'leading', 'oper', strEvent, 'delchan', 'off', ...
    'delevent', 'on');

chanids = [TTL];
eegplot(EEG.data(chanids,:), 'srate', EEG.srate, 'eloc_file',EEG.chanlocs(chanids), ...
    'winlength', 100, 'spacing', 10000000, 'events', EEG.event);


%% build events

%for wm_blocks 
trig_log = allEventInfo{18}{2}{5}; %allTrigLog{4};
diff_12 =  8076; %s1 = 24172 % s4_sess1 = 20362; % s4_sess2 = 40400 % s5 = 73482 
%s3b1 = 13040 
%sub03_sess1_b2 =29020
%sub05_sess2_b1 6664
%sub06_sess1_b1 16834
%sub06_sess1_b2 19785
%sub07_sess1_b1 21036 %sub07_sess1_b2= 7854 %sub07_sess1_b3= 11929
%sub07_sess2_b1= 7835 %sub07_sess2_b2= 5933 
%sub11_sess1_b1 13846
%sub11_sess1_b2 49861
%sub11_sess2_b1 7616.8
%sub11_sess2_b2  6168.8
%sub12_sess1_b1 61497
%sub12_sess1_b2 25824
%sub13_sess1_b1 6957.3
%sub13_sess1_b2 4406
%sub14_sess1_b1 38012
%sub14_sess1_b2 36989
%sub18_sess1_b1 17572
%sub18_sess1_b2 -32909
%sub19_sess1_b1 22498
%sub19_sess1_b2 2101.2
%sub20_sess1_b1 1568.6
%sub20_sess1_b2 59345
%sub21_sess1_b1 2490.4
%sub21_sess1_b2 464110
%sub22_sess1_b1 10988;
%sub22_sess1_b2 26472;
%sub22_sess1_b3 -11351;
%sub22_sess2_b1 1003;
%sub22_sess2_b2 6621.6
%sub23_sess1_b1 25509;
%sub23_sess1_b2 13282;
%sub24_sess1_b1 51033
%sub24_sess1_b2 2.3242e+04
%sub24_sess1_b3 1.2168e+04
%sub25_sess1_b1 1.7303e+04
%sub25_sess1_b2 1.0182e+04;
%sub25_sess2_b1 6.4681e+04
%sub25_sess2_b2 8076

EEG.event = struct('latency', [], 'type', '', 'urevent', []);
for i = 1:length(trig_log)
    EEG.event(length(EEG.event) + 1) = struct('latency', trig_log(i)  /5 - diff_12, 'type', ...
        num2str(trig_log(i,2:end)),... 
        'urevent', 0);
end
EEG.event(1) = [];

%EEG.event(length(EEG.event) + 1) = struct('latency',[[[[[[[[[31092]]]]]]]]], 'type', 'S', 'urevent', 0); 


%remove repeated latencies for single-item trials
idL = [EEG.event.latency];
[uni idx] = unique(idL)
EEG.event = EEG.event(idx);

EEG.event = nestedSortStruct(EEG.event, 'latency', 'type');


chanids = [TTL];
eegplot(EEG.data(chanids,:), 'srate', EEG.srate, 'eloc_file',EEG.chanlocs(chanids), ...
    'winlength', 50, 'spacing', 100000000, 'events', EEG.event);


%% 3) Downsample (first downsample then filter (much faster))
%https://es.mathworks.com/help/signal/ug/changing-signal-sample-rate.html
tic
[P,Q] = rat(1000/EEG.srate);

clear data
for chani = 1:size(EEG.data,1)
    chani
    data(chani,:) = resample(double(EEG.data(chani,:)), P, Q);
end

EEG.data = data;
EEG.srate = 1000;
EEG.pnts = size(EEG.data,2);
EEG.xmax = size(EEG.data,2)/EEG.srate;
%EEG.times = 0:2:size(EEG.data,2)*2;

toc

%% rebuild events

x = ([EEG.event.latency] / 2)'
x = num2cell(x')
EEG.event = struct('latency',x, 'type', {EEG.event.type});

%% plot 
close all
eventChannel = 'Trigger';TTL = strmatch(eventChannel, {EEG.chanlocs.labels}, 'exact');
disp ('--> plotting'); chanids=[TTL];
eegplot(EEG.data(chanids,:,:), 'srate', EEG.srate, 'eloc_file',EEG.chanlocs(chanids),...
    'winlength', 25, 'spacing', 10000000, 'events', EEG.event);


%% fft
%Perform FFT ckeck 

H = 'TAR1'; H = strmatch(H, {EEG.chanlocs.labels}, 'exact');
if (isempty(H)) disp ('Channel does not exist!'); return; end  

signal = EEG.data(H,1:length(EEG.data));
%signal = EEG.data(H,1000:100000);
N       = length(signal);    % length of sequence
nyquist = EEG.srate/2;           % Nyquist frequency 
%-- highest frequency that can be measured in the data
frequencies = linspace(0,nyquist,floor(N/2)+1);

% Compute fourier transform and scale
fourierCoefsF = fft(signal) / N;

figure; plot(frequencies,abs(fourierCoefsF(1:length(frequencies)))*2,'ro')
set(gca,'xlim',[2 170], 'ylim', [0 1])






%% 4) filter & notch
tic; disp (['starting at ' datestr(now, 'HH:MM:SS')]);
EEG = pop_eegfiltnew (EEG, 1, 150);
EEG = pop_eegfiltnew (EEG, 49, 51, [],  1); %notch filter
EEG = pop_eegfiltnew (EEG, 99, 101, [],  1); %notch filter
EEG = pop_eegfiltnew (EEG, 149, 151, [],  1); %notch filter
toc;
disp('filtering done')

%% 5) SAVE
%cd '/Users/danielpachecoSpecs/Google Drive/HospitalDelMar/dataAnalysis/data/experiments/recognition/recognitionLFP';
cd 'D:\_WM\data\iEEG\set - monopolar';
%%
save('s25_sess2_blo2_EEG', 'EEG', '-v7.3');




%%

x = EEG.event(2).type
strsplit(x)


%%
[oneline, cellarray]=cuixuFindStructure([20 6 10; 30 9 12]);


