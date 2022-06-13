
function [myPlotR] = plot_reinst_map_wm (reinst_plot_cfg)

%load myCmap;
myCmap = colormap(brewermap([],'YlOrRd'));

cfg                 =   reinst_plot_cfg;
meanReal_cond1      =   cfg.out_real.meanReal_cond1;
meanReal_cond2      =   cfg.out_real.meanReal_cond2;
sigMT_real          =   cfg.out_real.sigMT_real;
plot1clust          =   cfg.plot1clust; 
clust2plot          =   cfg.clust2plot;
lwd1                =   cfg.lwd1;%baseline
lwd2                =   cfg.lwd2; %cluster

myPlotR= figure(1); set (gcf, 'Position', [100 100 1500 500]);
subplot (131);
meanReal_cond1 = squeeze(mean(meanReal_cond1, 'omitnan' ));
meanReal_cond2 = squeeze(mean(meanReal_cond2, 'omitnan' ));

if cfg.dupSym
    meanReal_cond1 = triu(meanReal_cond1.',1) + tril(meanReal_cond1);
    meanReal_cond2 = triu(meanReal_cond2.',1) + tril(meanReal_cond2);
    sigMT_real = triu(sigMT_real.',1) + tril(sigMT_real);
    cfg.sigMH_thres = triu(cfg.sigMH_thres.',1) + tril(cfg.sigMH_thres);
end


if cfg.square
    imagesc (flipud(myresizem(meanReal_cond1, 10)));axis equal; 
        hold on; colorbar;
    if cfg.plotCueOnset
        plot(get(gca,'xlim'), [cfg.limFE*10 cfg.limFE*10],'k:', 'linewidth', lwd1); 
        plot([cfg.limFR*10 cfg.limFR*10], get(gca,'ylim'),'k:', 'linewidth', lwd1);  
    end
    if cfg.plotD
        plot([length(meanReal_cond1) 1] + [0.5 -1.5],get(gca,'xlim') + [0 1], 'k:'); % diagonal
    end
else
    contourf (meanReal_cond1, 500, 'linecolor', 'none'); axis equal; 
        hold on; colorbar;
    if cfg.plotCueOnset
        plot(get(gca,'xlim'), [cfg.limFE cfg.limFE],'k:', 'linewidth', lwd1); 
        plot([cfg.limFR cfg.limFR], get(gca,'ylim'),'k:', 'linewidth', lwd1);  
    end
    if cfg.plotD
        plot(get(gca,'xlim'), [1 length(meanReal_cond1)],'k:'); % diagonal
    end
end
t1 = title (cfg.lbls3plot{1});

subplot (132);
if cfg.square
    imagesc (flipud(myresizem(meanReal_cond2, 10))); axis equal; 
    hold on; colorbar;
    if cfg.plotCueOnset
        plot(get(gca,'xlim'), [cfg.limFE*10 cfg.limFE*10],'k:', 'linewidth', lwd1); 
        plot([cfg.limFR*10 cfg.limFR*10], get(gca,'ylim'),'k:', 'linewidth', lwd1); 
    end
    if cfg.plotD
        plot([length(meanReal_cond1) 1] + [0.5 -1.5],get(gca,'xlim') + [0 1], 'k:'); % diagonal
    end
else
contourf (meanReal_cond2, 500, 'linecolor', 'none'); axis equal; 
    hold on; colorbar;
if cfg.plotCueOnset
    plot(get(gca,'xlim'), [cfg.limFE cfg.limFE],'k:', 'linewidth', lwd1); 
    plot([cfg.limFR cfg.limFR], get(gca,'ylim'),'k:', 'linewidth', lwd1);  
end
if cfg.plotD
    plot(get(gca,'xlim'), [1 length(meanReal_cond1)],'k:'); % diagonal
end
end
t2 = title (cfg.lbls3plot{2});
colormap jet;


ax3 = subplot (133);
if cfg.square
    if plot1clust
        cfg.sigMH_thres(cfg.out_real.clustInfoReal.PixelIdxList{clust2plot}) = -100;
    end
    imagesc (flipud(myresizem(sigMT_real, 10))); axis equal;hold on; % colorbar; 
    
    % % % take contour
    contour(flipud(myresizem(cfg.sigMH_thres, 10)), 1, 'lineWidth', lwd2, 'linecolor', 'k');colorbar; %axis square; 
    
    
    % % % %draw a square in each pixel
%     h = find(flipud(cfg.sigMH_thres));
%     for ci = 1:length(h)
%         [row, col] = ind2sub(size(cfg.sigMH_thres), h(ci));
%         coordX =  [(col-.5)  (col +.5) (col+.5) (col-.5)    ] ;
%         coordY =  [(row-.5)  (row -.5) (row+.5) (row+.5)    ] ;%
%         patch(coordX, coordY, 'red', 'FaceColor','none', 'LineWidth', 1); axis equal; colorbar
%     end
    
    if cfg.plotCueOnset
        plot(get(gca,'xlim'), [cfg.limFE*10 cfg.limFE*10],'k:', 'linewidth', lwd1); 
        plot([cfg.limFR*10 cfg.limFR*10], get(gca,'ylim'),'k:', 'linewidth', lwd1); 
    end
else
    if strcmp (cfg.diff, 'tmap')
        if plot1clust
            cfg.sigMH_thres(cfg.out_real.clustInfoReal.PixelIdxList{clust2plot}) = -100;
        end
        contourf(sigMT_real, 40, 'linecolor', 'none') ; axis equal; 
        hold on; colorbar;
        contour(cfg.sigMH_thres, 1, 'lineWidth', lwd1, 'linecolor', 'k:');colorbar; %axis square;        
    end
    if strcmp (cfg.diff, 'diff')
        diff = meanReal_cond1 - meanReal_cond2;
        contourf(diff, 40, 'linecolor', 'none') ;axis equal; hold on; colorbar;   
        contour(cfg.sigMH_thres, 1, 'lineWidth', lwd1, 'linecolor', 'k:');colorbar; %axis square;
        
    end
    if cfg.plotCueOnset
        plot(get(gca,'xlim'), [cfg.limFE cfg.limFE],'k:', 'linewidth', lwd1);
        plot([cfg.limFR cfg.limFR], get(gca,'ylim'),'k:', 'linewidth', lwd1); 
    end
    if cfg.plotD
        plot(get(gca,'xlim'), [1 length(meanReal_cond1)],'k:'); % diagonal
    end
end

%t3 = title (cfg.lbls3plot{3}); 
colormap (ax3, myCmap);

axesHandles = findall(0, 'type', 'axes');
labels_to_plotE = num2cell(cfg.labels_to_plotE); 
labels_to_plotR = num2cell(cfg.labels_to_plotR); 
labels_to_plotE(cfg.l2excE) = {' '};
labels_to_plotR(cfg.l2excR) = {' '};

d2uE = cfg.binsE; d2uR = cfg.binsR;
limE = [0 d2uE];
limR = [0 d2uR];
if cfg.square
    set(axesHandles, 'clim', cfg.clim, 'ytick', cfg.placeTY*10, 'yticklabel', fliplr(labels_to_plotE),...
    'xtick', cfg.placeTX*10, 'xticklabel', labels_to_plotR, 'FontSize', 12, 'xlim', limE*10, 'ylim', limR*10); 
    %set(gca,'XTick',[], 'YTick', [], 'xticklabel',[])
else
    set(axesHandles, 'clim', cfg.clim, 'ytick', cfg.placeTY, 'yticklabel', labels_to_plotE,...
    'xtick', cfg.placeTX, 'xticklabel', labels_to_plotR, 'FontSize', 12, 'xlim', limE, 'ylim', limR); 

    %set(t1,'Position',get(t1,'Position')+[0 1 0]);
    %set(t2,'Position',get(t2,'Position')+[0 1 0]);
    %set(t3,'Position',get(t3,'Position')+[0 1 0]);
    
end

if strcmp (cfg.diff, 'tmap') % if plotting differences let matlab set this value
    set(ax3, 'clim', cfg.climT);
end


for i=1:length(axesHandles)
    %xlabel(axesHandles(i),{'Retrieval time (s)'})
    %ylabel(axesHandles(i),{'Encoding time (s)'})
    xlabel(axesHandles(i),{''})
    ylabel(axesHandles(i),{''})
end



%set (gcf, 'Position', [300 300 820 500]);

if cfg.saveimg
export_fig(2, cfg.imageName,'-transparent', '-r300');
%exportgraphics(gcf, cfg.imageName, 'Resolution',150);
close all;
end


    
% figure(); hist (max_clust_sum); hold on; 
% plot([max_clust_sum_real max_clust_sum_real],get(gca,'ylim'),'g')
% plot([max_clust_sum_perm max_clust_sum_perm],get(gca,'ylim'),'r')

end









