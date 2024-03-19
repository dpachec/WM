

function [out_real] = plot_reinst_map_wm_no_plot (cfg)

runPerm = strsplit(cfg.res, '_'); 
if strcmp(runPerm{2}, 'perm')
    cfg.runperm = 1;
else
    cfg.runperm = 0;
end
 
%exclude subjects
if cfg.subj2exc > 0
    cfg.all_cond1(cfg.subj2exc) = []; 
    cfg.all_cond2(cfg.subj2exc) = []; 
end

[cfg_plot]   =      set_reinst_plot_wm (cfg);
 
% apply limits
for si = 1:size(cfg.all_cond1, 1);   
    cfg.all_cond1{si} = cfg.all_cond1{si}(:,cfg_plot.mlimE,cfg_plot.mlimR);
    cfg.all_cond2{si} = cfg.all_cond2{si}(:,cfg_plot.mlimE,cfg_plot.mlimR);
end
 

%calculate real differences
[out_real]         =   real_diff_reinst_wm(cfg);

cfg_plot.lwd1        =       2;%baseline
cfg_plot.lwd2        =       2; %cluster
cfg_plot.plotClust   =       1; %cluster
%cfg_plot.sigMH_thres =       out_real.sigMH_real;



end









