%% SUB 1 
%% load edf
clearvars;
%eeglab;

%sub01
cd /Users/danielpacheco/Documents/iEEG_data_analysis/extinction/data/raw_data/c_sub01/ieeg
EEG = pop_biosig('DBX-Learning test.edf', 'importevent', 'off'); 

%cd 'D:\GCLRs\Stanford_Navigation_data\DATA\raw_data\S21_172\E21-1363_sess03-0124_object\behavior'



%% plot 
chanids = [49:49];
eegplot(EEG.data(chanids,:,:), 'srate', EEG.srate, 'winlength', 10, 'spacing', 5000000, 'events', EEG.event);

%%
cfg        = [];
cfg.metric = 'zvalue';  % use by default zvalue method
cfg.method = 'summary'; % use by default summary method
data       = ft_rejectvisual(cfg,data);

%% 
