function [cfg_plot] = set_reinst_plot_wm (cfg)
    
cfg_plot             =       [];
cfg_plot.remClust    =       cfg.remClust;
cfg_plot.clim        =       cfg.clim;
cfg_plot.climT       =       cfg.climT;
cfg_plot.saveimg     =       cfg.saveimg;
cfg_plot.plot1clust  =       cfg.plot1clust;
cfg_plot.clust2plot  =       cfg.clust2plot;

if strcmp(cfg.res, '10_norm')
    if strcmp(cfg.cut2, '4-4')
        cfg_plot.mlimE = 1:550; %1:85;
        cfg_plot.mlimR = 1:550; %1:85
        x = 15; %note that x is different from the original size of the matrix 
        cfg_plot.binsE = x; cfg_plot.binsR = x; 
        cfg_plot.labels_to_plotE = -.5:.5:1; 
        cfg_plot.labels_to_plotR = -.5:.5:1; 
        cfg_plot.limFE = 10.5; %cue onset line
        cfg_plot.limFR = 4.5;
        cfg_plot.l2excE = [1 3 5 7 9  ]; %[2 4 6 8 10 12 14 16 18 20 22 24 26 28]; 
        cfg_plot.l2excR = [1 3 5 7 9 11 13 15 17 19];%[2 4 6 8 10 12 14 16 18 20 22 24 26 28]; 
        cfg_plot.placeTX = [0 4.5 9.5 14.5]; 
        cfg_plot.placeTY = [0.5 5.5 10.5];  
        cfg_plot.plotCueOnset    =       1; %plot cue onset at encoding and retrieval
    end
end

if strcmp(cfg.res, '100_norm')
    if strcmp(cfg.cut2, '1-1')
        cfg_plot.mlimE = 3:17; %1:85;
        cfg_plot.mlimR = 3:17; %1:85
        x = 15; %note that x is different from the original size of the matrix 
        cfg_plot.binsE = x; cfg_plot.binsR = x; 
        cfg_plot.labels_to_plotE = -.5:.5:1; 
        cfg_plot.labels_to_plotR = -.5:.5:1; 
        cfg_plot.limFE = 10.5; %cue onset line
        cfg_plot.limFR = 4.5;
        cfg_plot.l2excE = [1 3 5 7 9  ]; %[2 4 6 8 10 12 14 16 18 20 22 24 26 28]; 
        cfg_plot.l2excR = [1 3 5 7 9 11 13 15 17 19];%[2 4 6 8 10 12 14 16 18 20 22 24 26 28]; 
        cfg_plot.placeTX = [0 4.5 9.5 14.5]; 
        cfg_plot.placeTY = [0.5 5.5 10.5];  
        cfg_plot.plotCueOnset    =       1; %plot cue onset at encoding and retrieval
    end

    if strcmp(cfg.cut2, '1-4')
        cfg_plot.mlimE = 3:17; %1:85;
        cfg_plot.mlimR = 3:47; %1:85
        x = 45; %note that x is different from the original size of the matrix 
        cfg_plot.binsE = 45; 
        cfg_plot.binsR = 15; 
        cfg_plot.labels_to_plotE = -.5:.5:1; 
        cfg_plot.labels_to_plotR = -.5:.5:4; 
        cfg_plot.limFE =  10.5;
        cfg_plot.limFR =  4.5;
        cfg_plot.l2excE = [1 3 ]; %[2 4 6 8 10 12 14 16 18 20 22 24 26 28]; 
        cfg_plot.l2excR = [1 3 5 7 9 11 13 15 17 19];%[2 4 6 8 10 12 14 16 18 20 22 24 26 28]; 
        cfg_plot.placeTX = [0.5 5.5 10.5 15.5 20.5 25.5 30.5 35.5 40.5 45.5]; 
        cfg_plot.placeTY = [0.5 5.5 10.5 15.5 20.5 25.5 30.5 35.5 40.5 45.5];  
        cfg_plot.plotCueOnset    =       1; %plot cue onset at encoding and retrieval
    end

    if strcmp(cfg.cut2, '4-4')
        cfg_plot.mlimE = 3:47 %1:45; %1:85;
        cfg_plot.mlimR = 3:47 %1:45; %1:85
        x = 45; %45 %note that x is different from the original size of the matrix 
        cfg_plot.binsE = x; cfg_plot.binsR = x; 
        cfg_plot.labels_to_plotE = -.5:.5:3.5; 
        cfg_plot.labels_to_plotR = -.5:.5:3.5; 
        cfg_plot.limFE = 40.5; %cue onset line
        cfg_plot.limFR = 4.5;
        cfg_plot.l2excE = [1 3 5 7 9  ]; %[2 4 6 8 10 12 14 16 18 20 22 24 26 28]; 
        cfg_plot.l2excR = [1 3 5 7 9 11 13 15 17 19];%[2 4 6 8 10 12 14 16 18 20 22 24 26 28]; 
        cfg_plot.placeTX = [0.5 4.5 9.5 14.5 19.5 24.5 29.5 34.5 39.5 44.5]; 
        cfg_plot.placeTY = [0.5 5.5 10.5 15.5 20.5 25.5 30.5 35.5 40.5 45.5];
        cfg_plot.plotCueOnset    =       1; %plot cue onset at encoding and retrieval
                   
    end
end 



if strcmp(cfg.res, '100_perm')
    if strcmp(cfg.cut2, '1-1')
        cfg_plot.mlimE = 8:17;
        cfg_plot.mlimR = 8:17;
        x = 10; %note that x is different from the original size of the matrix 
        cfg_plot.binsE = x; cfg_plot.binsR = x; 
        cfg_plot.labels_to_plotE = 0:.5:1; 
        cfg_plot.labels_to_plotR = 0:.5:1; 
        cfg_plot.limFE = 0; %cue onset line
        cfg_plot.limFR = 0;
        cfg_plot.l2excE = [ ]; %[2 4 6 8 10 12 14 16 18 20 22 24 26 28]; 
        cfg_plot.l2excR = [ ];%[2 4 6 8 10 12 14 16 18 20 22 24 26 28]; 
        cfg_plot.placeTX = [0.5 5.5 ]; %first is bottom left (horizontal axis)
        cfg_plot.placeTY = [0.5 5.5 ];  %left is top left (vertical axis) last is bottom left (vertical axis)
        cfg_plot.plotCueOnset    =       0; %plot cue onset at encoding and retrieval

    end

    if strcmp(cfg.cut2, '1-4')
        cfg_plot.mlimE = 8:17; %1:85;
        cfg_plot.mlimR = 8:47; %1:85
        cfg_plot.binsE = 40; 
        cfg_plot.binsR = 10; 
        cfg_plot.labels_to_plotE = 0:.5:1; 
        cfg_plot.labels_to_plotR = 0:.5:4; 
        cfg_plot.limFE =  0;
        cfg_plot.limFR =  0;
        cfg_plot.l2excE = [2 4 6 8 10 12 ]; %[2 4 6 8 10 12 14 16 18 20 22 24 26 28]; 
        cfg_plot.l2excR = [1 3 5 7 9 11 13 15 17 19];%[2 4 6 8 10 12 14 16 18 20 22 24 26 28]; 
        cfg_plot.placeTX = []; 
        cfg_plot.placeTY = [0.5 5.5 10.5 15.5 20.5 25.5 30.5 35.5 40.5 45.5];  
        cfg_plot.plotCueOnset    =       0; %plot cue onset at encoding and retrieval
        
    end

    if strcmp(cfg.cut2, '4-4')
        cfg_plot.mlimE = 8:42; %official
        cfg_plot.mlimR = 8:42; %official

        x = 40; %note that x is different from the original size of the matrix 
        cfg_plot.binsE = x; cfg_plot.binsR = x; 
        cfg_plot.labels_to_plotE = 0:.5:4; 
        cfg_plot.labels_to_plotR = 0:.5:4; 
        cfg_plot.limFE = 0; %cue onset line
        cfg_plot.limFR = 0;
        cfg_plot.l2excE = [1 3 5 7 9  ]; %[2 4 6 8 10 12 14 16 18 20 22 24 26 28]; 
        cfg_plot.l2excR = [1 3 5 7 9 11 13 15 17 19];%[2 4 6 8 10 12 14 16 18 20 22 24 26 28]; 
        cfg_plot.placeTX = [0.5 5.5 10.5 15.5 20.5 25.5 30.5 35.5 40.5 45.5]; 
        cfg_plot.placeTY = [0.5 5.5 10.5 15.5 20.5 25.5 30.5 35.5 40.5 45.5];
        cfg_plot.plotCueOnset    =       0; %plot cue onset at encoding and retrieval
        
           
    end
end 


    
    
% % % name of the images to be saved depending on conditions

if strcmp(cfg.cond1, 'SI') & strcmp(cfg.cond2,'DI')
    cfg_plot.lbls3plot = {['Same Item'] ['Different Item'] ['Same Item vs.' newline 'Different Item']};
    cfg_plot.imageName = ['_SI-DI-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) cfg_plot.imageName = ['_SI-DI-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'SIALL') & strcmp(cfg.cond2,'SI123')
    cfg_plot.lbls3plot = {['SIALL'] ['SI123'] ['SIALL vs.' newline 'SI123']};
    cfg_plot.imageName = ['_SIALL-SI123-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) cfg_plot.imageName = ['_SIALL-SI123-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'SI1') & strcmp(cfg.cond2,'SI2')
    cfg_plot.lbls3plot = {['SI1'] ['SI2'] ['SI1 vs.' newline 'SI2']};
    cfg_plot.imageName = ['_SI1-SI2-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) cfg_plot.imageName = ['_SI1-SI2-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'SIE') & strcmp(cfg.cond2,'DISC')
    cfg_plot.lbls3plot = {['SIE'] ['DISC'] ['SIE vs.' newline 'DISC']};
    cfg_plot.imageName = ['_SIE-DISC-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) cfg_plot.imageName = ['_SIE-DISC-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'SIE') & strcmp(cfg.cond2,'SIALL')
    cfg_plot.lbls3plot = {['SIE'] ['SIALL'] ['SIE vs.' newline 'SIALL']};
    cfg_plot.imageName = ['_SIE-SIALL-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) cfg_plot.imageName = ['_SIE-SIALL-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'SI1') & strcmp(cfg.cond2,'DISC1')
    cfg_plot.lbls3plot = {['SI1'] ['DISC1'] ['SI1 vs.' newline 'DISC1']};
    cfg_plot.imageName = ['_SI1-DISC1-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) cfg_plot.imageName = ['_SI1-DISC1-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'SI2') & strcmp(cfg.cond2,'DISC2')
    cfg_plot.lbls3plot = {['SI2'] ['DISC2'] ['SI2 vs.' newline 'DISC2']};
    cfg_plot.imageName = ['_SI2-DISC2-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) cfg_plot.imageName = ['_SI2-DISC2-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'SI3') & strcmp(cfg.cond2,'DISC3')
    cfg_plot.lbls3plot = {['SI3'] ['DISC3'] ['SI3 vs.' newline 'DISC3']};
    cfg_plot.imageName = ['_SI3-DISC3-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) cfg_plot.imageName = ['_SI3-DISC3-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'SISC') & strcmp(cfg.cond2,'DISC')
    cfg_plot.lbls3plot = {['SISC'] ['DISC'] ['SISC vs.' newline 'DISC']};
    cfg_plot.imageName = ['_SISC-DISC-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) cfg_plot.imageName = ['_SISC-DISC-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'DISC') & strcmp(cfg.cond2,'DIDC')
    cfg_plot.lbls3plot = {['DISC'] ['DIDC'] ['DISC vs.' newline 'DIDC']};
    cfg_plot.imageName = ['_DISC-DIDC-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) cfg_plot.imageName = ['_DISC-DIDC-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'SISC') & strcmp(cfg.cond2,'DIDC')
    cfg_plot.lbls3plot = {['SISC'] ['DIDC'] ['SISC vs.' newline 'DIDC']};
    cfg_plot.imageName = ['_SISC-DIDC-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) cfg_plot.imageName = ['_SISC-DIDC-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'SISC_EE') & strcmp(cfg.cond2,'DIDC_EE')
    cfg_plot.lbls3plot = {['SISC'] ['DIDC'] ['SISC vs.' newline 'DIDC']};
    cfg_plot.imageName = ['_SISC_EE-DIDC_EE-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) cfg_plot.imageName = ['_SISC_EE-DIDC_EE-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'SISC_EE') & strcmp(cfg.cond2,'DIDC_EE')
    cfg_plot.lbls3plot = {['SISC'] ['DIDC'] ['SISC vs.' newline 'DIDC']};
    cfg_plot.imageName = ['_SISC_EE-DIDC_EE-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) cfg_plot.imageName = ['_SISC_EE-DIDC_EE-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'SISC_EEC') & strcmp(cfg.cond2,'DISC_EEC')
    cfg_plot.lbls3plot = {['SISC'] ['DISC'] ['SISC vs.' newline 'DISC']};
    cfg_plot.imageName = ['_SISC_EEC-DISC_EEC-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) cfg_plot.imageName = ['_SISC_EEC-DISC_EEC-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'DISC_EE') & strcmp(cfg.cond2,'DIDC_EE')
    cfg_plot.lbls3plot = {['DISC'] ['DIDC'] ['DISC vs.' newline 'DIDC']};
    cfg_plot.imageName = ['_DISC_EE-DIDC_EE-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) cfg_plot.imageName = ['_DISC_EE-DIDC_EE-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end

if strcmp(cfg.cond1, 'SISC_M0') & strcmp(cfg.cond2,'DIDC_M0')
    cfg_plot.lbls3plot = {['SISC'] ['DIDC'] ['SISC vs.' newline 'DIDC']};
    cfg_plot.imageName = ['_SISC_M0-DIDC_M0-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) cfg_plot.imageName = ['_SISC_M0-DIDC_M0-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'SISC_M0') & strcmp(cfg.cond2,'DISC_M0')
    cfg_plot.lbls3plot = {['SISC'] ['DISC'] ['SISC vs.' newline 'DISC']};
    cfg_plot.imageName = ['_SISC_M0-DISC_M0-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) cfg_plot.imageName = ['_SISC_M0-DISC_M0-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'DISC_M0') & strcmp(cfg.cond2,'DIDC_M0')
    cfg_plot.lbls3plot = {['DISC'] ['DIDC'] ['DISC vs.' newline 'DIDC']};
    cfg_plot.imageName = ['_DISC_M0-DIDC_M0-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) cfg_plot.imageName = ['_DISC_M0-DIDC_M0-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end

if strcmp(cfg.cond1, 'SISC_EM2') & strcmp(cfg.cond2,'DISC_EM2')
    cfg_plot.lbls3plot = {['SISC'] ['DISC'] ['SISC vs.' newline 'DISC']};
    cfg_plot.imageName = ['_SISC_EM2-DISC_EM2-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) cfg_plot.imageName = ['_SISC_EM2-DISC_EM2-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'DISC_EM1') & strcmp(cfg.cond2,'DIDC_EM1')
    cfg_plot.lbls3plot = {['DISC'] ['DIDC'] ['DISC vs.' newline 'DIDC']};
    cfg_plot.imageName = ['_DISC_EM1-DIDC_EM1-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) cfg_plot.imageName = ['_DISC_EM1-DIDC_EM1-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'DISC_EM2') & strcmp(cfg.cond2,'DIDC_EM2')
    cfg_plot.lbls3plot = {['DISC'] ['DIDC'] ['DISC vs.' newline 'DIDC']};
    cfg_plot.imageName = ['_DISC_EM2-DIDC_EM2-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) cfg_plot.imageName = ['_DISC_EM2-DIDC_EM2-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end


if strcmp(cfg.cond1, 'SISC_EM2U') & strcmp(cfg.cond2,'DISC_EM2U')
    cfg_plot.lbls3plot = {['SISCU'] ['DISCU'] ['SISCU vs.' newline 'DISCU']};
    cfg_plot.imageName = ['_SISC_EM2U-DISC_EM2U-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) cfg_plot.imageName = ['_SISC_EM2U-DISC_EM2U-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'DISC_EM2U') & strcmp(cfg.cond2,'DIDC_EM2U')
    cfg_plot.lbls3plot = {['DISCU'] ['DIDCU'] ['DISCU vs.' newline 'DIDCU']};
    cfg_plot.imageName = ['_DISC_EM2U-DIDC_EM2U-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) cfg_plot.imageName = ['_DISC_EM2U-DIDC_EM2U-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end


if strcmp(cfg.cond1, 'SISC_M2M2') & strcmp(cfg.cond2,'DIDC_M2M2')
    cfg_plot.lbls3plot = {['SISC'] ['DIDC'] ['SISC vs.' newline 'DIDC']};
    cfg_plot.imageName = ['_SISC_M2M2-DIDC_M2M2-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) cfg_plot.imageName = ['_SISC_M2M2-DIDC_M2M2-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'SISC_M2M2') & strcmp(cfg.cond2,'DISC_M2M2')
    cfg_plot.lbls3plot = {['SISC'] ['DISC'] ['SISC vs.' newline 'DISC']};
    cfg_plot.imageName = ['_SISC_M2M2-DISC_M2M2-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) cfg_plot.imageName = ['_SISC_M2M2-DISC_M2M2-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'DISC_M2M2') & strcmp(cfg.cond2,'DIDC_M2M2')
    cfg_plot.lbls3plot = {['DISC'] ['DIDC'] ['DISC vs.' newline 'DIDC']};
    cfg_plot.imageName = ['_DISC_M2M2-DIDC_M2M2-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) cfg_plot.imageName = ['_DISC_M2M2-DIDC_M2M2-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'SISC_EECC') & strcmp(cfg.cond2,'DISC_EECC')
    cfg_plot.lbls3plot = {['SISC'] ['DISC'] ['SISC vs.' newline 'DISC']};
    cfg_plot.imageName = ['_SISC-DISC-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) cfg_plot.imageName = ['_SISC_EECC-DISC_EECC-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'DISC_EECC') & strcmp(cfg.cond2,'DIDC_EECC')
    cfg_plot.lbls3plot = {['DISC'] ['DIDC'] ['DISC vs.' newline 'DIDC']};
    cfg_plot.imageName = ['_DISC-DIDC-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) cfg_plot.imageName = ['_DISC_EECC-DIDC_EECC-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'SISC_EECI') & strcmp(cfg.cond2,'DISC_EECI')
    cfg_plot.lbls3plot = {['SISC'] ['DISC'] ['SISC vs.' newline 'DISC']};
    cfg_plot.imageName = ['_SISC-DISC-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) cfg_plot.imageName = ['_SISC_EECI-DISC_EECI-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'DISC_EECI') & strcmp(cfg.cond2,'DIDC_EECI')
    cfg_plot.lbls3plot = {['DISC'] ['DIDC'] ['DISC vs.' newline 'DIDC']};
    cfg_plot.imageName = ['_DISC-DIDC-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) cfg_plot.imageName = ['_DISC_EECI-DIDC_EECI-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'SISC_EEI') & strcmp(cfg.cond2,'DISC_EEI')
    cfg_plot.lbls3plot = {['SISC'] ['DISC'] ['SISC vs.' newline 'DISC']};
    cfg_plot.imageName = ['_SISC-DISC-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) cfg_plot.imageName = ['_SISC_EEI-DISC_EEI-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'DISC_EEI') & strcmp(cfg.cond2,'DIDC_EEI')
    cfg_plot.lbls3plot = {['DISC'] ['DIDC'] ['DISC vs.' newline 'DIDC']};
    cfg_plot.imageName = ['_DISC-DIDC-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) cfg_plot.imageName = ['_DISC_EEI-DIDC_EEI-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'SISC_EECC') & strcmp(cfg.cond2,'SISC_EEI')
    cfg_plot.lbls3plot = {['SISC_EECC'] ['SISC_EEI'] ['SISC_EECC vs.' newline 'SISC_EEI']};
    cfg_plot.imageName = ['_SISC_EECC-SISC_EEI-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) cfg_plot.imageName = ['_SISC_EECC-SISC_EEI-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'SISC_NCM2') & strcmp(cfg.cond2,'DISC_NCM2')
    cfg_plot.lbls3plot = {['SISC'] ['DISC'] ['SISC vs.' newline 'DISC']};
    cfg_plot.imageName = ['_SISC_NCM2-DISC_NCM2-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) cfg_plot.imageName = ['_SISC_NCM2-DISC_NCM2-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'DISC_NCM2') & strcmp(cfg.cond2,'DIDC_NCM2')
    cfg_plot.lbls3plot = {['DISC'] ['DIDC'] ['DISC vs.' newline 'DIDC']};
    cfg_plot.imageName = ['_DISC_NCM2-DIDC_NCM2-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) cfg_plot.imageName = ['_DISC_NCM2-DIDC_NCM2-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'DISC_M2A') & strcmp(cfg.cond2,'DIDC_M2A')
    cfg_plot.lbls3plot = {['DISC'] ['DIDC'] ['DISC vs.' newline 'DIDC']};
    cfg_plot.imageName = ['_DISC_M2A-DIDC_M2A-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) cfg_plot.imageName = ['_DISC_M2A-DIDC_M2A-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'CI') & strcmp(cfg.cond2,'I')
    cfg_plot.lbls3plot = {['Correct'] ['Incorrect'] ['Correct vs.' newline 'Incorrect']};
    cfg_plot.imageName = ['_CI-I-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) cfg_plot.imageName = ['_CI-I-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'CC') & strcmp(cfg.cond2,'I')
    cfg_plot.lbls3plot = {['Correct'] ['Incorrect'] ['Correct vs.' newline 'Incorrect']};
    cfg_plot.imageName = ['_CC-I-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) cfg_plot.imageName = ['_CC-I-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'DISC_123M2A') & strcmp(cfg.cond2,'DIDC_123M2A')
    cfg_plot.lbls3plot = {['DISC'] ['DIDC'] ['DISC vs.' newline 'DIDC']};
    cfg_plot.imageName = ['_DISC_123M2A-DIDC_123M2A-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) cfg_plot.imageName = ['_DISC_123M2A-DIDC_123M2A-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end

if strcmp(cfg.cond1, 'DISC_M2123V1') & strcmp(cfg.cond2,'DIDC_M2123V1')
    cfg_plot.lbls3plot = {['DISC'] ['DIDC'] ['DISC vs.' newline 'DIDC']};
    cfg_plot.imageName = ['_DISC_M2123V1-DIDC_M2123V1-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) cfg_plot.imageName = ['_DISC_M2123V1-DIDC_M2123V1-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'DISC_M2123V2') & strcmp(cfg.cond2,'DIDC_M2123V2')
    cfg_plot.lbls3plot = {['DISC'] ['DIDC'] ['DISC vs.' newline 'DIDC']};
    cfg_plot.imageName = ['_DISC_M2123V2-DIDC_M2123V2-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) cfg_plot.imageName = ['_DISC_M2123V2-DIDC_M2123V2-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'DISC_M2123CNCV1') & strcmp(cfg.cond2,'DIDC_M2123CNCV1')
    cfg_plot.lbls3plot = {['DISC'] ['DIDC'] ['DISC vs.' newline 'DIDC']};
    cfg_plot.imageName = ['_DISC_M2123CNCV1-DIDC_M2123CNCV1-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) cfg_plot.imageName = ['_DISC_M2123CNCV1-DIDC_M2123CNCV1-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'DISC_M2123CNCV2') & strcmp(cfg.cond2,'DIDC_M2123CNCV2')
    cfg_plot.lbls3plot = {['DISC'] ['DIDC'] ['DISC vs.' newline 'DIDC']};
    cfg_plot.imageName = ['_DISC_M2123CNCV2-DIDC_M2123CNCV2-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) cfg_plot.imageName = ['_DISC_M2123CNCV2-DIDC_M2123CNCV2-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'DISC_M2123NC') & strcmp(cfg.cond2,'DIDC_M2123NC')
    cfg_plot.lbls3plot = {['DISC'] ['DIDC'] ['DISC vs.' newline 'DIDC']};
    cfg_plot.imageName = ['_DISC_M2123NC-DIDC_M2123NC-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) cfg_plot.imageName = ['_DISC_M2123NC-DIDC_M2123NC-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'DISC_M2A123') & strcmp(cfg.cond2,'DIDC_M2A123')
    cfg_plot.lbls3plot = {['DISC'] ['DIDC'] ['DISC vs.' newline 'DIDC']};
    cfg_plot.imageName = ['_DISC_M2A123-DIDC_M2A123-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) cfg_plot.imageName = ['_DISC_M2A123-DIDC_M2A123-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end

if strcmp(cfg.cond1, 'DISC_EM2UV1') & strcmp(cfg.cond2,'DIDC_EM2UV1')
    cfg_plot.lbls3plot = {['DISC'] ['DIDC'] ['DISC vs.' newline 'DIDC']};
    cfg_plot.imageName = ['_DISC_EM2UV1-DIDC_EM2UV1-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) cfg_plot.imageName = ['_DISC_EM2UV1-DIDC_EM2UV1-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'SISC_EM2UV2') & strcmp(cfg.cond2,'DISC_EM2UV2')
    cfg_plot.lbls3plot = {['SISC'] ['DISC'] ['SISC vs.' newline 'DISC']};
    cfg_plot.imageName = ['_SISC_EM2UV2-DISC_EM2UV2-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) cfg_plot.imageName = ['_SISC_EM2UV2-_DISC_EM2UV2-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'DISC_EM2UV2') & strcmp(cfg.cond2,'DIDC_EM2UV2')
    cfg_plot.lbls3plot = {['DISC'] ['DIDC'] ['DISC vs.' newline 'DIDC']};
    cfg_plot.imageName = ['_DISC_EM2UV2-DIDC_EM2UV2-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) cfg_plot.imageName = ['_DISC_EM2UV2-DIDC_EM2UV2-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'DISC_M1A') & strcmp(cfg.cond2,'DIDC_M1A')
    cfg_plot.lbls3plot = {['DISC'] ['DIDC'] ['DISC vs.' newline 'DIDC']};
    cfg_plot.imageName = ['_DISC_M1A-DIDC_M1A-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) cfg_plot.imageName = ['_DISC_M1A-DIDC_M1A-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'SCSP_M2M2') & strcmp(cfg.cond2,'SCDP_M2M2')
    cfg_plot.lbls3plot = {['SCSP'] ['SCDP'] ['SCSP vs.' newline 'SCDP']};
    cfg_plot.imageName = ['_SCSP_M2M2-SCDP_M2M2-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) cfg_plot.imageName = ['_SCSP_M2M2-SCDP_M2M2-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end




cfg_plot = cfg_plot;
 
 
end

 

