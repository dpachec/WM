%% load edf
clearvars;
tic
%eeglab;
cd F:\Marie\WM_all_data\raw_data\WM_datasets\jings_data\iEEG\sub1_zhengweiwei

%EEG = pop_biosig('sub1_zhengweiwei', 'importevent', 'off'); 
EEG = pop_loadbv('.', 'wm_block4.vhdr');
toc

%% plot (without filtering values can be off in the eegplot)
chanids = [1:3];
eegplot(EEG.data(chanids,:,:), 'srate', EEG.srate, 'winlength', 100, 'spacing', 1500, 'events', EEG.event);



%% 3) Downsample to 2000Hz
%https://es.mathworks.com/help/signal/ug/changing-signal-sample-rate.html
tic
[P,Q] = rat(2000/EEG.srate);

clear data
for chani = 1:size(EEG.data,1)
    chani
    data(chani,:) = resample(double(EEG.data(chani,:)), P, Q);
end

EEG.data = data;
EEG.srate = 2000;
EEG.pnts = size(EEG.data,2);
EEG.xmax = size(EEG.data,2)/EEG.srate;
%EEG.times = 0:2:size(EEG.data,2)*2;



%%rebuild events
x = ([EEG.event.latency] *0.8)'
x = num2cell(x')
EEG.event = struct('latency',x, 'type', '----');


events_jin = EEG.event;


%%get events from log files

[allEventInfo alleeg_event ] = loadLogsWM;
allEventInfo = allEventInfo';

toc


%%

%%build events
EEG.event = events_jin;
%for wm_blocks 
trig_log = allEventInfo{28}{1}{10}; %allTrigLog{4};
diff_12 = -34050; 
%s51_sess1_b1 = -25189; 
%s51_sess1_b2 = -23465; 
%s51_sess1_b3 = -37847; 
%s51_sess1_b4 = -47219; 
%s52_sess1_b1 = -23451; 
%s52_sess1_b2 = 306270; 
%s52_sess1_b3 = -13927; 
%s53_sess1_b1 = -15666; 
%s53_sess1_b2 = 69429; 
%s53_sess1_b3 = -59014; 
%s53_sess1_b4 = -63686; 
%s54_sess1_b1 = -67823; 
%s54_sess1_b2 = -12403; 
%s54_sess1_b3 = -26549; 
%s54_sess1_b4 = -70942; 
%s55_sess1_b1 = -24069; 
%s55_sess1_b2 = -43360; 
%s55_sess1_b3 = -91856; 
%s55_sess1_b3 = -34050; 

EEG.event = struct('latency', [], 'type', []);
for i = 1:length(trig_log)
%     EEG.event(length(EEG.event) + 1) = struct('latency', trig_log(i)  /2 - diff_12, ...
%         'duration', 1, 'channel', 0, 'bvtime', [], 'bvmknum', [],'visible', [], 'type', num2str(trig_log(i,2:end)), 'code', 'Stimulus',...
%         'urevent', 0);
    EEG.event(length(EEG.event) + 1) = struct('latency', trig_log(i)  /5 - diff_12, 'type', num2str(trig_log(i,2:end)));


end

EEG.event(1) = [];

%EEG.event(length(EEG.event) + 1) = struct('latency',[[26886]], 'type', 'S', 'urevent', 0); 

%remove repeated latencies for single-item trials
idL = [EEG.event.latency];
[uni idx] = unique(idL)
EEG.event = EEG.event(idx);

EEG.event = nestedSortStruct(EEG.event, 'latency', 'type');


chanids = [1];
eegplot(EEG.data(chanids,:), 'srate', EEG.srate, 'eloc_file',EEG.chanlocs(chanids), ...
    'winlength', 50, 'spacing', 1000, 'events', EEG.event);


%% 3) Downsample to 1000
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

%%rebuild events
x = ([EEG.event.latency] / 2)'
x = num2cell(x')
EEG.event = struct('latency',x, 'type', {EEG.event.type});

%% plot 
close all
eventChannel = '1';TTL = strmatch(eventChannel, {EEG.chanlocs.labels}, 'exact');
disp ('--> plotting'); chanids=[TTL];
eegplot(EEG.data(chanids,:,:), 'srate', EEG.srate, 'eloc_file',EEG.chanlocs(chanids),...
    'winlength', 50, 'spacing', 1000, 'events', EEG.event);


%% fft
%Perform FFT ckeck 

H = '10'; H = strmatch(H, {EEG.chanlocs.labels}, 'exact');
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
save('s55_sess1_blo4_EEG', 'EEG');%, '-v7.3'


toc

%% 






%%






%%