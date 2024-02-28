

function [out_real] = plot_reinst_map_wm (cfg)

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


if strcmp(cfg.cond1(end-2:end), 'EM2') | strcmp(cfg.cond1(end-2:end), 'EM1')
    dupSym      =      0; 
else
    dupSym      =      1; 
end




%load myCmap;
myCmap = colormap(brewermap([],'YlOrRd'));

meanReal_cond1      =   out_real.meanReal_cond1;
meanReal_cond2      =   out_real.meanReal_cond2;
sigMT_real          =   out_real.sigMT_real;
sigMH_real          =   out_real.sigMH_real;
plot1clust          =   cfg.plot1clust; 
clust2plot          =   cfg.clust2plot;
lwd1                =   cfg.lwd1;%baseline
lwd2                =   cfg.lwd2; %cluster

myPlotR= figure(1); set (gcf, 'Position', [100 100 1500 500]);

meanReal_cond1 = squeeze(mean(meanReal_cond1, 'omitnan' ));
meanReal_cond2 = squeeze(mean(meanReal_cond2, 'omitnan' ));

if dupSym
    meanReal_cond1 = triu(meanReal_cond1.',1) + tril(meanReal_cond1);
    meanReal_cond2 = triu(meanReal_cond2.',1) + tril(meanReal_cond2);
    sigMT_real = triu(sigMT_real.',1) + tril(sigMT_real);
    sigMH_thres = triu(sigMH_real.',1) + tril(sigMH_real);
    sigMH_real = triu(sigMH_real.',1) + tril(sigMH_real);
end


% % % uncomment to plot all 3 maps
% % subplot (131);
% % imagesc (flipud(myresizem(meanReal_cond1, 10)));axis equal; hold on; colorbar;
% % if cfg_plot.plotCueOnset
% %     plot(get(gca,'xlim'), [cfg_plot.limFE*10 cfg_plot.limFE*10],'k:', 'linewidth', lwd1); 
% %     plot([cfg_plot.limFR*10 cfg_plot.limFR*10], get(gca,'ylim'),'k:', 'linewidth', lwd1);  
% % end
% % t1 = title (cfg_plot.lbls3plot{1});
% % 
% % subplot (132);
% % imagesc (flipud(myresizem(meanReal_cond2, 10))); axis equal; hold on; colorbar;
% % if cfg_plot.plotCueOnset
% %     plot(get(gca,'xlim'), [cfg_plot.limFE*10 cfg_plot.limFE*10],'k:', 'linewidth', lwd1); 
% %     plot([cfg_plot.limFR*10 cfg_plot.limFR*10], get(gca,'ylim'),'k:', 'linewidth', lwd1); 
% % end
% % t2 = title (cfg_plot.lbls3plot{2});
% % colormap jet;


ax3 = subplot (133);

if plot1clust
    sigMH_real = zeros(size(sigMH_real)); 
    for pixi = 1:length(clust2plot)
        sigMH_real(out_real.clustInfoReal.PixelIdxList{clust2plot(pixi)}) = 1;
    end
    if dupSym
        sigMH_real = triu(sigMH_real.',1) + tril(sigMH_real);
    end
end
imagesc (flipud(myresizem(sigMT_real, 10))); axis equal;hold on; 
if cfg_plot.plotClust
    contour(flipud(myresizem(sigMH_real, 10)), 1, 'lineWidth', lwd2, 'linecolor', 'k'); %colorbar; 
end
if cfg_plot.plotCueOnset
    plot(get(gca,'xlim'), [cfg_plot.limFE*10 cfg_plot.limFE*10],'k:', 'linewidth', lwd1); 
    plot([cfg_plot.limFR*10 cfg_plot.limFR*10], get(gca,'ylim'),'k:', 'linewidth', lwd1); 
end
%t3 = title (cfg_plot.lbls3plot{3}); 
colormap (ax3, myCmap);

axesHandles = findall(0, 'type', 'axes');
labels_to_plotE = num2cell(cfg_plot.labels_to_plotE); 
labels_to_plotR = num2cell(cfg_plot.labels_to_plotR); 
labels_to_plotE(cfg_plot.l2excE) = {' '};
labels_to_plotR(cfg_plot.l2excR) = {' '};

d2uE = cfg_plot.binsE; d2uR = cfg_plot.binsR;
limE = [0 d2uE];
limR = [0 d2uR];

%set(axesHandles, 'clim', cfg_plot.clim, 'ytick', cfg_plot.placeTY*10, 'yticklabel', fliplr(labels_to_plotE),...
%'xtick', (cfg_plot.placeTX-1)*10, 'xticklabel', labels_to_plotR, 'FontSize', 12, 'xlim', limE*10, 'ylim', limR*10); 
set(axesHandles, 'clim', cfg_plot.clim, 'ytick', cfg_plot.placeTY*10, 'yticklabel', {},...
'xtick', (cfg_plot.placeTX-1)*10, 'xticklabel', {}, 'FontSize', 12, 'xlim', limE*10, 'ylim', limR*10); 

set(ax3, 'clim', cfg.climT);



for i=1:length(axesHandles)
    xlabel(axesHandles(i),{''})
    ylabel(axesHandles(i),{''})
end

if cfg.saveimg
%export_fig(2, cfg_plot.imageName,'-transparent', '-r300');
exportgraphics(gcf, cfg_plot.imageName, 'Resolution',150);
%close all;
end


end









