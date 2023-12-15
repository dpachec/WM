%% load edf
clearvars;
tic
eeglab;
cd 'C:\Users\Neuropsychology\Desktop\wm09_m_20210119_hangzhou\wm09_m_20210119_hangzhou\iEEG\'

EEG = pop_biosig('JY20210119.edf', 'importevent', 'off'); 
toc

%% plot (without filtering values can be off in the eegplot)
eventChannel = 'POL DC13';TTL = strmatch(eventChannel, {EEG.chanlocs.labels}, 'exact');
disp ('--> plotting'); chanids=[TTL];
eegplot(EEG.data(chanids,:,:), 'srate', EEG.srate, 'eloc_file',EEG.chanlocs(chanids),...
    'winlength', 100, 'spacing', 15000000);



%% 1) FIRST REMOVE EMPTY CHANNELS
chans2rem = [38:49 51:53 158:277];
EEG.data(chans2rem, :) = [];
EEG.chanlocs(chans2rem) = [];
EEG.nbchan = size(EEG.data,1);

%% 2) and non experimental data

EEG.data = EEG.data (:,2620*EEG.srate:3750*EEG.srate);
%s1= 2534:3570%s4sess1=1647:2637, s4_sess2 =1513:2578, s5 = 979:2006
%s3=1832:4292
EEG.pnts = length(EEG.data);

disp('data extracted');





%% Extract triggers from TTL channel
eventChannel = 'POL DC13';TTL = strmatch(eventChannel, {EEG.chanlocs.labels}, 'exact');
strEvent = 'X < -162000';

EEG.event = [];
EEG = pop_chanevent (EEG, TTL, 'edge', 'leading', 'oper', strEvent, 'delchan', 'off', ...
    'delevent', 'on');

chanids = [TTL];
eegplot(EEG.data(chanids,:), 'srate', EEG.srate, 'eloc_file',EEG.chanlocs(chanids), ...
    'winlength', 25, 'spacing', 150000, 'events', EEG.event);

%% get events from log files

[allEventInfo alleeg_event ] = loadLogsWM;
allEventInfo = allEventInfo';



%% build events

%for wm_blocks 
trig_log = allEventInfo{25}{1}{6}; %allTrigLog{4};
diff_12 =  13969;
%s31_sess1_b1 =27278  
%s31_sess1_b2 =17857  
%s32_sess1_b1 = 12628  
%s32_sess1_b2 = 64796
%s33_sess1_b1 = 15068
%s33_sess1_b2 = 7418.8
%s34_sess1_b1 = 29531
%s34_sess1_b2 = 25955
%s35_sess1_b1 = 25851
%s35_sess1_b2 = 15358
%s36_sess1_b1 = 32825
%s36_sess1_b2 = 46942
%s37_sess1_b1 = -3900.4
%s37_sess1_b2 = -11319
%s38_sess1_b1 = -14026
%s38_sess1_b2 = -4045.6; 
%s39_sess1_b1 = -9210.6;
%s39_sess1_b2 = 13969;

EEG.event = struct('latency', [], 'type', '', 'urevent', []);
for i = 1:length(trig_log)
    EEG.event(length(EEG.event) + 1) = struct('latency', trig_log(i)  /5 - diff_12, 'type', ...
        num2str(trig_log(i,2:end)),... 
        'urevent', 0);
end
EEG.event(1) = [];

%EEG.event(length(EEG.event) + 1) = struct('latency',[[[46066]]], 'type', 'S', 'urevent', 0); 

%remove repeated latencies for single-item trials
idL = [EEG.event.latency];
[uni idx] = unique(idL)
EEG.event = EEG.event(idx);

EEG.event = nestedSortStruct(EEG.event, 'latency', 'type');


chanids = [TTL];
eegplot(EEG.data(chanids,:), 'srate', EEG.srate, 'eloc_file',EEG.chanlocs(chanids), ...
    'winlength', 50, 'spacing', 10000000, 'events', EEG.event);


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
eventChannel = 'POL DC13';TTL = strmatch(eventChannel, {EEG.chanlocs.labels}, 'exact');
disp ('--> plotting'); chanids=[TTL];
eegplot(EEG.data(chanids,:,:), 'srate', EEG.srate, 'eloc_file',EEG.chanlocs(chanids),...
    'winlength', 25, 'spacing', 10000000, 'events', EEG.event);


%% fft
%Perform FFT ckeck 

H = 'POL B1'; H = strmatch(H, {EEG.chanlocs.labels}, 'exact');
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

%% 5) SAVE
cd 'D:\_WM\data\iEEG\set - monopolar';


%%
tic
disp('saving')
save('s25_sess1_blo2_EEG', 'EEG');%, '-v7.3'


toc

%% 



%%

x = EEG.event(2).type
strsplit(x)


%%
[oneline, cellarray]=cuixuFindStructure([20 6 10; 30 9 12]);


%%






%%






%%