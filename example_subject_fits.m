%% first load power data of one subject 

clear , clc
load example_power_data.mat   


%% perform fit of Alexnet representations during the maintenance period

cfg.period      = 'M123'; 
cfg.avRep       = 1; 
cfg.freqs       = 3:54; 
cfg.win_width   = 5; 
cfg.mf          = 1; 
cfg.fR          = 1; 
cfg.avTFV       = 0; 
cfg.brainROI    = 'pfc';
cfg.net2load    = 'Alex';
cfg.lays2load   = 1:8; 
cfg.fitMode     = 0; % trials (i.e., separately for each RSM row, or full RSM)

paths.stim = '/Users/danielpacheco/Documents/iEEG_data_analysis/WM/stimuli/'; % indicate here path to image stim

cfg_contrasts               = getIdsWM(cfg.period, cfg_contrasts);
cfg_contrasts               = average_repetitions(cfg, cfg_contrasts);
neuralRDMs                  = createNeuralRDMs(cfg, cfg_contrasts);
networkRDMs                 = createNetworkRDMs(cfg, cfg_contrasts, 14, paths);
nnFit                       = fitModel_WM(neuralRDMs, networkRDMs, cfg.fitMode); 



%% plot layer 8 time-frequency fit 
figure()
data4Plot = squeeze(nnFit(1, :, :)); 

times = 1:46; 
freqs = 1:52; 
contourf(times, freqs, data4Plot, 100, 'linecolor', 'none');  %colorbar




%%
























