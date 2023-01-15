function [cfg_plot] = set_reinst_plot_wm (cfg)
    
reinst_plot_cfg             =       [];
reinst_plot_cfg.remClust    =       cfg.remClust;
reinst_plot_cfg.clim        =       cfg.clim;
reinst_plot_cfg.climT       =       cfg.climT;
reinst_plot_cfg.square      =       cfg.square;
reinst_plot_cfg.saveimg     =       cfg.saveimg;
reinst_plot_cfg.plot1clust  =       cfg.plot1clust;
reinst_plot_cfg.clust2plot  =       cfg.clust2plot;

if strcmp(cfg.res, '10_norm')
    if strcmp(cfg.cut2, '4-4')
        reinst_plot_cfg.mlimE = 1:550; %1:85;
        reinst_plot_cfg.mlimR = 1:550; %1:85
        x = 15; %note that x is different from the original size of the matrix 
        reinst_plot_cfg.binsE = x; reinst_plot_cfg.binsR = x; 
        reinst_plot_cfg.labels_to_plotE = -.5:.5:1; 
        reinst_plot_cfg.labels_to_plotR = -.5:.5:1; 
        reinst_plot_cfg.limFE = 10.5; %cue onset line
        reinst_plot_cfg.limFR = 4.5;
        reinst_plot_cfg.l2excE = [1 3 5 7 9  ]; %[2 4 6 8 10 12 14 16 18 20 22 24 26 28]; 
        reinst_plot_cfg.l2excR = [1 3 5 7 9 11 13 15 17 19];%[2 4 6 8 10 12 14 16 18 20 22 24 26 28]; 
        reinst_plot_cfg.placeTX = [0 4.5 9.5 14.5]; 
        reinst_plot_cfg.placeTY = [0.5 5.5 10.5];  
        reinst_plot_cfg.plotCueOnset    =       1; %plot cue onset at encoding and retrieval
        reinst_plot_cfg.plotD           =       1;
    end
end

if strcmp(cfg.res, '100_norm')
    if strcmp(cfg.cut2, '1-1')
        reinst_plot_cfg.mlimE = 3:17; %1:85;
        reinst_plot_cfg.mlimR = 3:17; %1:85
        x = 15; %note that x is different from the original size of the matrix 
        reinst_plot_cfg.binsE = x; reinst_plot_cfg.binsR = x; 
        reinst_plot_cfg.labels_to_plotE = -.5:.5:1; 
        reinst_plot_cfg.labels_to_plotR = -.5:.5:1; 
        reinst_plot_cfg.limFE = 10.5; %cue onset line
        reinst_plot_cfg.limFR = 4.5;
        reinst_plot_cfg.l2excE = [1 3 5 7 9  ]; %[2 4 6 8 10 12 14 16 18 20 22 24 26 28]; 
        reinst_plot_cfg.l2excR = [1 3 5 7 9 11 13 15 17 19];%[2 4 6 8 10 12 14 16 18 20 22 24 26 28]; 
        reinst_plot_cfg.placeTX = [0 4.5 9.5 14.5]; 
        reinst_plot_cfg.placeTY = [0.5 5.5 10.5];  
        reinst_plot_cfg.plotCueOnset    =       1; %plot cue onset at encoding and retrieval
        reinst_plot_cfg.plotD           =       1;
    end

    if strcmp(cfg.cut2, '1-4')
        reinst_plot_cfg.mlimE = 3:17; %1:85;
        reinst_plot_cfg.mlimR = 3:47; %1:85
        x = 45; %note that x is different from the original size of the matrix 
        reinst_plot_cfg.binsE = 45; 
        reinst_plot_cfg.binsR = 15; 
        reinst_plot_cfg.labels_to_plotE = -.5:.5:1; 
        reinst_plot_cfg.labels_to_plotR = -.5:.5:4; 
        reinst_plot_cfg.limFE =  10.5;
        reinst_plot_cfg.limFR =  4.5;
        reinst_plot_cfg.l2excE = [1 3 ]; %[2 4 6 8 10 12 14 16 18 20 22 24 26 28]; 
        reinst_plot_cfg.l2excR = [1 3 5 7 9 11 13 15 17 19];%[2 4 6 8 10 12 14 16 18 20 22 24 26 28]; 
        reinst_plot_cfg.placeTX = [0.5 5.5 10.5 15.5 20.5 25.5 30.5 35.5 40.5 45.5]; 
        reinst_plot_cfg.placeTY = [0.5 5.5 10.5 15.5 20.5 25.5 30.5 35.5 40.5 45.5];  
        reinst_plot_cfg.plotCueOnset    =       1; %plot cue onset at encoding and retrieval
        reinst_plot_cfg.plotD           =       0;
    end

    if strcmp(cfg.cut2, '4-4')
        reinst_plot_cfg.mlimE = 3:47 %1:45; %1:85;
        reinst_plot_cfg.mlimR = 3:47 %1:45; %1:85
        x = 45; %45 %note that x is different from the original size of the matrix 
        reinst_plot_cfg.binsE = x; reinst_plot_cfg.binsR = x; 
        reinst_plot_cfg.labels_to_plotE = -.5:.5:4; 
        reinst_plot_cfg.labels_to_plotR = -.5:.5:4; 
        reinst_plot_cfg.limFE = 40.5; %cue onset line
        reinst_plot_cfg.limFR = 4.5;
        reinst_plot_cfg.l2excE = [1 3 5 7 9  ]; %[2 4 6 8 10 12 14 16 18 20 22 24 26 28]; 
        reinst_plot_cfg.l2excR = [1 3 5 7 9 11 13 15 17 19];%[2 4 6 8 10 12 14 16 18 20 22 24 26 28]; 
        reinst_plot_cfg.placeTX = [0.5 5.5 10.5 15.5 20.5 25.5 30.5 35.5 40.5 45.5]; 
        reinst_plot_cfg.placeTY = [0.5 5.5 10.5 15.5 20.5 25.5 30.5 35.5 40.5 45.5];
        reinst_plot_cfg.plotCueOnset    =       1; %plot cue onset at encoding and retrieval
        reinst_plot_cfg.plotD           =       1;
           
    end
end 



if strcmp(cfg.res, '100_perm')
    if strcmp(cfg.cut2, '1-1')
        reinst_plot_cfg.mlimE = 8:15; %1:85;
        reinst_plot_cfg.mlimR = 8:15; %1:85
        x = 10; %note that x is different from the original size of the matrix 
        reinst_plot_cfg.binsE = x; reinst_plot_cfg.binsR = x; 
        reinst_plot_cfg.labels_to_plotE = 0:.5:1; 
        reinst_plot_cfg.labels_to_plotR = 0:.5:1; 
        reinst_plot_cfg.limFE = 0; %cue onset line
        reinst_plot_cfg.limFR = 0;
        reinst_plot_cfg.l2excE = [ ]; %[2 4 6 8 10 12 14 16 18 20 22 24 26 28]; 
        reinst_plot_cfg.l2excR = [ ];%[2 4 6 8 10 12 14 16 18 20 22 24 26 28]; 
        reinst_plot_cfg.placeTX = [0.5 5.5 ]; %first is bottom left (horizontal axis)
        reinst_plot_cfg.placeTY = [0.5 5.5 ];  %left is top left (vertical axis) last is bottom left (vertical axis)
        reinst_plot_cfg.plotCueOnset    =       0; %plot cue onset at encoding and retrieval
        reinst_plot_cfg.plotD           =       0;
    end

    if strcmp(cfg.cut2, '1-4')
        reinst_plot_cfg.mlimE = 8:17; %1:85;
        reinst_plot_cfg.mlimR = 8:42; %1:85
        reinst_plot_cfg.binsE = 40; 
        reinst_plot_cfg.binsR = 8; 
        reinst_plot_cfg.labels_to_plotE = 0:.5:1; 
        reinst_plot_cfg.labels_to_plotR = 0:.5:4; 
        reinst_plot_cfg.limFE =  0;
        reinst_plot_cfg.limFR =  0;
        reinst_plot_cfg.l2excE = [2 4 6 8 10 12 ]; %[2 4 6 8 10 12 14 16 18 20 22 24 26 28]; 
        reinst_plot_cfg.l2excR = [1 3 5 7 9 11 13 15 17 19];%[2 4 6 8 10 12 14 16 18 20 22 24 26 28]; 
        reinst_plot_cfg.placeTX = []; 
        reinst_plot_cfg.placeTY = [0.5 5.5 10.5 15.5 20.5 25.5 30.5 35.5 40.5 45.5];  
        reinst_plot_cfg.plotCueOnset    =       0; %plot cue onset at encoding and retrieval
        reinst_plot_cfg.plotD           =       0;
    end

    if strcmp(cfg.cut2, '4-4')
        reinst_plot_cfg.mlimE = 6:45; %1:85;
        reinst_plot_cfg.mlimR = 6:45; %1:85
        x = 40; %note that x is different from the original size of the matrix 
        reinst_plot_cfg.binsE = x; reinst_plot_cfg.binsR = x; 
        reinst_plot_cfg.labels_to_plotE = 0:.5:4; 
        reinst_plot_cfg.labels_to_plotR = 0:.5:4; 
        reinst_plot_cfg.limFE = 0; %cue onset line
        reinst_plot_cfg.limFR = 0;
        reinst_plot_cfg.l2excE = [1 3 5 7 9  ]; %[2 4 6 8 10 12 14 16 18 20 22 24 26 28]; 
        reinst_plot_cfg.l2excR = [1 3 5 7 9 11 13 15 17 19];%[2 4 6 8 10 12 14 16 18 20 22 24 26 28]; 
        reinst_plot_cfg.placeTX = [0.5 5.5 10.5 15.5 20.5 25.5 30.5 35.5 40.5 45.5]; 
        reinst_plot_cfg.placeTY = [0.5 5.5 10.5 15.5 20.5 25.5 30.5 35.5 40.5 45.5];
        reinst_plot_cfg.plotCueOnset    =       0; %plot cue onset at encoding and retrieval
        reinst_plot_cfg.plotD           =       0;
           
    end
end 





%reinst_plot_cfg.binsE = length(reinst_plot_cfg.mlimE); 
%reinst_plot_cfg.binsR = length(reinst_plot_cfg.mlimR);

    
    

if strcmp(cfg.cond1, 'SI') & strcmp(cfg.cond2,'DI')
    reinst_plot_cfg.lbls3plot = {['Same Item'] ['Different Item'] ['Same Item vs.' newline 'Different Item']};
    reinst_plot_cfg.imageName = ['_SI-DI-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) reinst_plot_cfg.imageName = ['_SI-DI-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'SIALL') & strcmp(cfg.cond2,'SI123')
    reinst_plot_cfg.lbls3plot = {['SIALL'] ['SI123'] ['SIALL vs.' newline 'SI123']};
    reinst_plot_cfg.imageName = ['_SIALL-SI123-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) reinst_plot_cfg.imageName = ['_SIALL-SI123-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'SI1') & strcmp(cfg.cond2,'SI2')
    reinst_plot_cfg.lbls3plot = {['SI1'] ['SI2'] ['SI1 vs.' newline 'SI2']};
    reinst_plot_cfg.imageName = ['_SI1-SI2-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) reinst_plot_cfg.imageName = ['_SI1-SI2-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'SIE') & strcmp(cfg.cond2,'DISC')
    reinst_plot_cfg.lbls3plot = {['SIE'] ['DISC'] ['SIE vs.' newline 'DISC']};
    reinst_plot_cfg.imageName = ['_SIE-DISC-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) reinst_plot_cfg.imageName = ['_SIE-DISC-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'SIE') & strcmp(cfg.cond2,'SIALL')
    reinst_plot_cfg.lbls3plot = {['SIE'] ['SIALL'] ['SIE vs.' newline 'SIALL']};
    reinst_plot_cfg.imageName = ['_SIE-SIALL-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) reinst_plot_cfg.imageName = ['_SIE-SIALL-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'SI1') & strcmp(cfg.cond2,'DISC1')
    reinst_plot_cfg.lbls3plot = {['SI1'] ['DISC1'] ['SI1 vs.' newline 'DISC1']};
    reinst_plot_cfg.imageName = ['_SI1-DISC1-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) reinst_plot_cfg.imageName = ['_SI1-DISC1-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'SI2') & strcmp(cfg.cond2,'DISC2')
    reinst_plot_cfg.lbls3plot = {['SI2'] ['DISC2'] ['SI2 vs.' newline 'DISC2']};
    reinst_plot_cfg.imageName = ['_SI2-DISC2-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) reinst_plot_cfg.imageName = ['_SI2-DISC2-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'SI3') & strcmp(cfg.cond2,'DISC3')
    reinst_plot_cfg.lbls3plot = {['SI3'] ['DISC3'] ['SI3 vs.' newline 'DISC3']};
    reinst_plot_cfg.imageName = ['_SI3-DISC3-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) reinst_plot_cfg.imageName = ['_SI3-DISC3-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'SISC') & strcmp(cfg.cond2,'DISC')
    reinst_plot_cfg.lbls3plot = {['SISC'] ['DISC'] ['SISC vs.' newline 'DISC']};
    reinst_plot_cfg.imageName = ['_SISC-DISC-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) reinst_plot_cfg.imageName = ['_SISC-DISC-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'DISC') & strcmp(cfg.cond2,'DIDC')
    reinst_plot_cfg.lbls3plot = {['DISC'] ['DIDC'] ['DISC vs.' newline 'DIDC']};
    reinst_plot_cfg.imageName = ['_DISC-DIDC-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) reinst_plot_cfg.imageName = ['_DISC-DIDC-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'SISC') & strcmp(cfg.cond2,'DIDC')
    reinst_plot_cfg.lbls3plot = {['SISC'] ['DIDC'] ['SISC vs.' newline 'DIDC']};
    reinst_plot_cfg.imageName = ['_SISC-DIDC-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) reinst_plot_cfg.imageName = ['_SISC-DIDC-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'SISC_EE') & strcmp(cfg.cond2,'DIDC_EE')
    reinst_plot_cfg.lbls3plot = {['SISC'] ['DIDC'] ['SISC vs.' newline 'DIDC']};
    reinst_plot_cfg.imageName = ['_SISC_EE-DIDC_EE-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) reinst_plot_cfg.imageName = ['_SISC_EE-DIDC_EE-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'SISC_EE') & strcmp(cfg.cond2,'DIDC_EE')
    reinst_plot_cfg.lbls3plot = {['SISC'] ['DIDC'] ['SISC vs.' newline 'DIDC']};
    reinst_plot_cfg.imageName = ['_SISC_EE-DIDC_EE-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) reinst_plot_cfg.imageName = ['_SISC_EE-DIDC_EE-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'SISC_EEC') & strcmp(cfg.cond2,'DISC_EEC')
    reinst_plot_cfg.lbls3plot = {['SISC'] ['DISC'] ['SISC vs.' newline 'DISC']};
    reinst_plot_cfg.imageName = ['_SISC_EEC-DISC_EEC-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) reinst_plot_cfg.imageName = ['_SISC_EEC-DISC_EEC-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'DISC_EE') & strcmp(cfg.cond2,'DIDC_EE')
    reinst_plot_cfg.lbls3plot = {['DISC'] ['DIDC'] ['DISC vs.' newline 'DIDC']};
    reinst_plot_cfg.imageName = ['_DISC_EE-DIDC_EE-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) reinst_plot_cfg.imageName = ['_DISC_EE-DIDC_EE-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end

if strcmp(cfg.cond1, 'SISC_M0') & strcmp(cfg.cond2,'DIDC_M0')
    reinst_plot_cfg.lbls3plot = {['SISC'] ['DIDC'] ['SISC vs.' newline 'DIDC']};
    reinst_plot_cfg.imageName = ['_SISC_M0-DIDC_M0-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) reinst_plot_cfg.imageName = ['_SISC_M0-DIDC_M0-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'SISC_M0') & strcmp(cfg.cond2,'DISC_M0')
    reinst_plot_cfg.lbls3plot = {['SISC'] ['DISC'] ['SISC vs.' newline 'DISC']};
    reinst_plot_cfg.imageName = ['_SISC_M0-DISC_M0-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) reinst_plot_cfg.imageName = ['_SISC_M0-DISC_M0-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'DISC_M0') & strcmp(cfg.cond2,'DIDC_M0')
    reinst_plot_cfg.lbls3plot = {['DISC'] ['DIDC'] ['DISC vs.' newline 'DIDC']};
    reinst_plot_cfg.imageName = ['_DISC_M0-DIDC_M0-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) reinst_plot_cfg.imageName = ['_DISC_M0-DIDC_M0-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end

if strcmp(cfg.cond1, 'SISC_EM2') & strcmp(cfg.cond2,'DISC_EM2')
    reinst_plot_cfg.lbls3plot = {['SISC'] ['DISC'] ['SISC vs.' newline 'DISC']};
    reinst_plot_cfg.imageName = ['_SISC_EM2-DISC_EM2-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) reinst_plot_cfg.imageName = ['_SISC_EM2-DISC_EM2-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'DISC_EM2') & strcmp(cfg.cond2,'DIDC_EM2')
    reinst_plot_cfg.lbls3plot = {['DISC'] ['DIDC'] ['DISC vs.' newline 'DIDC']};
    reinst_plot_cfg.imageName = ['_DISC_EM2-DIDC_EM2-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) reinst_plot_cfg.imageName = ['_DISC_EM2-DIDC_EM2-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end


if strcmp(cfg.cond1, 'SISC_EM2U') & strcmp(cfg.cond2,'DISC_EM2U')
    reinst_plot_cfg.lbls3plot = {['SISCU'] ['DISCU'] ['SISCU vs.' newline 'DISCU']};
    reinst_plot_cfg.imageName = ['_SISC_EM2U-DISC_EM2U-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) reinst_plot_cfg.imageName = ['_SISC_EM2U-DISC_EM2U-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'DISC_EM2U') & strcmp(cfg.cond2,'DIDC_EM2U')
    reinst_plot_cfg.lbls3plot = {['DISCU'] ['DIDCU'] ['DISCU vs.' newline 'DIDCU']};
    reinst_plot_cfg.imageName = ['_DISC_EM2U-DIDC_EM2U-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) reinst_plot_cfg.imageName = ['_DISC_EM2U-DIDC_EM2U-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end


if strcmp(cfg.cond1, 'SISC_M2M2') & strcmp(cfg.cond2,'DIDC_M2M2')
    reinst_plot_cfg.lbls3plot = {['SISC'] ['DIDC'] ['SISC vs.' newline 'DIDC']};
    reinst_plot_cfg.imageName = ['_SISC_M2M2-DIDC_M2M2-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) reinst_plot_cfg.imageName = ['_SISC_M2M2-DIDC_M2M2-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'SISC_M2M2') & strcmp(cfg.cond2,'DISC_M2M2')
    reinst_plot_cfg.lbls3plot = {['SISC'] ['DISC'] ['SISC vs.' newline 'DISC']};
    reinst_plot_cfg.imageName = ['_SISC_M2M2-DISC_M2M2-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) reinst_plot_cfg.imageName = ['_SISC_M2M2-DISC_M2M2-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'DISC_M2M2') & strcmp(cfg.cond2,'DIDC_M2M2')
    reinst_plot_cfg.lbls3plot = {['DISC'] ['DIDC'] ['DISC vs.' newline 'DIDC']};
    reinst_plot_cfg.imageName = ['_DISC_M2M2-DIDC_M2M2-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) reinst_plot_cfg.imageName = ['_DISC_M2M2-DIDC_M2M2-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'SISC_EECC') & strcmp(cfg.cond2,'DISC_EECC')
    reinst_plot_cfg.lbls3plot = {['SISC'] ['DISC'] ['SISC vs.' newline 'DISC']};
    reinst_plot_cfg.imageName = ['_SISC-DISC-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) reinst_plot_cfg.imageName = ['_SISC_EECC-DISC_EECC-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'DISC_EECC') & strcmp(cfg.cond2,'DIDC_EECC')
    reinst_plot_cfg.lbls3plot = {['DISC'] ['DIDC'] ['DISC vs.' newline 'DIDC']};
    reinst_plot_cfg.imageName = ['_DISC-DIDC-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) reinst_plot_cfg.imageName = ['_DISC_EECC-DIDC_EECC-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'SISC_EECI') & strcmp(cfg.cond2,'DISC_EECI')
    reinst_plot_cfg.lbls3plot = {['SISC'] ['DISC'] ['SISC vs.' newline 'DISC']};
    reinst_plot_cfg.imageName = ['_SISC-DISC-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) reinst_plot_cfg.imageName = ['_SISC_EECI-DISC_EECI-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'DISC_EECI') & strcmp(cfg.cond2,'DIDC_EECI')
    reinst_plot_cfg.lbls3plot = {['DISC'] ['DIDC'] ['DISC vs.' newline 'DIDC']};
    reinst_plot_cfg.imageName = ['_DISC-DIDC-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) reinst_plot_cfg.imageName = ['_DISC_EECI-DIDC_EECI-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'SISC_EEI') & strcmp(cfg.cond2,'DISC_EEI')
    reinst_plot_cfg.lbls3plot = {['SISC'] ['DISC'] ['SISC vs.' newline 'DISC']};
    reinst_plot_cfg.imageName = ['_SISC-DISC-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) reinst_plot_cfg.imageName = ['_SISC_EEI-DISC_EEI-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'DISC_EEI') & strcmp(cfg.cond2,'DIDC_EEI')
    reinst_plot_cfg.lbls3plot = {['DISC'] ['DIDC'] ['DISC vs.' newline 'DIDC']};
    reinst_plot_cfg.imageName = ['_DISC-DIDC-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) reinst_plot_cfg.imageName = ['_DISC_EEI-DIDC_EEI-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'SISC_EECC') & strcmp(cfg.cond2,'SISC_EEI')
    reinst_plot_cfg.lbls3plot = {['SISC_EECC'] ['SISC_EEI'] ['SISC_EECC vs.' newline 'SISC_EEI']};
    reinst_plot_cfg.imageName = ['_SISC_EECC-SISC_EEI-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) reinst_plot_cfg.imageName = ['_SISC_EECC-SISC_EEI-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'SISC_NCM2') & strcmp(cfg.cond2,'DISC_NCM2')
    reinst_plot_cfg.lbls3plot = {['SISC'] ['DISC'] ['SISC vs.' newline 'DISC']};
    reinst_plot_cfg.imageName = ['_SISC_NCM2-DISC_NCM2-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) reinst_plot_cfg.imageName = ['_SISC_NCM2-DISC_NCM2-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'DISC_NCM2') & strcmp(cfg.cond2,'DIDC_NCM2')
    reinst_plot_cfg.lbls3plot = {['DISC'] ['DIDC'] ['DISC vs.' newline 'DIDC']};
    reinst_plot_cfg.imageName = ['_DISC_NCM2-DIDC_NCM2-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) reinst_plot_cfg.imageName = ['_DISC_NCM2-DIDC_NCM2-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'DISC_M2A') & strcmp(cfg.cond2,'DIDC_M2A')
    reinst_plot_cfg.lbls3plot = {['DISC'] ['DIDC'] ['DISC vs.' newline 'DIDC']};
    reinst_plot_cfg.imageName = ['_DISC_M2A-DIDC_M2A-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) reinst_plot_cfg.imageName = ['_DISC_M2A-DIDC_M2A-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'CI') & strcmp(cfg.cond2,'I')
    reinst_plot_cfg.lbls3plot = {['Correct'] ['Incorrect'] ['Correct vs.' newline 'Incorrect']};
    reinst_plot_cfg.imageName = ['_CI-I-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) reinst_plot_cfg.imageName = ['_CI-I-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'CC') & strcmp(cfg.cond2,'I')
    reinst_plot_cfg.lbls3plot = {['Correct'] ['Incorrect'] ['Correct vs.' newline 'Incorrect']};
    reinst_plot_cfg.imageName = ['_CC-I-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) reinst_plot_cfg.imageName = ['_CC-I-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'DISC_123M2A') & strcmp(cfg.cond2,'DIDC_123M2A')
    reinst_plot_cfg.lbls3plot = {['DISC'] ['DIDC'] ['DISC vs.' newline 'DIDC']};
    reinst_plot_cfg.imageName = ['_DISC_123M2A-DIDC_123M2A-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) reinst_plot_cfg.imageName = ['_DISC_123M2A-DIDC_123M2A-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end

if strcmp(cfg.cond1, 'DISC_M2123V1') & strcmp(cfg.cond2,'DIDC_M2123V1')
    reinst_plot_cfg.lbls3plot = {['DISC'] ['DIDC'] ['DISC vs.' newline 'DIDC']};
    reinst_plot_cfg.imageName = ['_DISC_M2123V1-DIDC_M2123V1-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) reinst_plot_cfg.imageName = ['_DISC_M2123V1-DIDC_M2123V1-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'DISC_M2123V2') & strcmp(cfg.cond2,'DIDC_M2123V2')
    reinst_plot_cfg.lbls3plot = {['DISC'] ['DIDC'] ['DISC vs.' newline 'DIDC']};
    reinst_plot_cfg.imageName = ['_DISC_M2123V2-DIDC_M2123V2-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) reinst_plot_cfg.imageName = ['_DISC_M2123V2-DIDC_M2123V2-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'DISC_M2123CNCV1') & strcmp(cfg.cond2,'DIDC_M2123CNCV1')
    reinst_plot_cfg.lbls3plot = {['DISC'] ['DIDC'] ['DISC vs.' newline 'DIDC']};
    reinst_plot_cfg.imageName = ['_DISC_M2123CNCV1-DIDC_M2123CNCV1-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) reinst_plot_cfg.imageName = ['_DISC_M2123CNCV1-DIDC_M2123CNCV1-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'DISC_M2123CNCV2') & strcmp(cfg.cond2,'DIDC_M2123CNCV2')
    reinst_plot_cfg.lbls3plot = {['DISC'] ['DIDC'] ['DISC vs.' newline 'DIDC']};
    reinst_plot_cfg.imageName = ['_DISC_M2123CNCV2-DIDC_M2123CNCV2-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) reinst_plot_cfg.imageName = ['_DISC_M2123CNCV2-DIDC_M2123CNCV2-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'DISC_M2123NC') & strcmp(cfg.cond2,'DIDC_M2123NC')
    reinst_plot_cfg.lbls3plot = {['DISC'] ['DIDC'] ['DISC vs.' newline 'DIDC']};
    reinst_plot_cfg.imageName = ['_DISC_M2123NC-DIDC_M2123NC-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) reinst_plot_cfg.imageName = ['_DISC_M2123NC-DIDC_M2123NC-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'DISC_M2A123') & strcmp(cfg.cond2,'DIDC_M2A123')
    reinst_plot_cfg.lbls3plot = {['DISC'] ['DIDC'] ['DISC vs.' newline 'DIDC']};
    reinst_plot_cfg.imageName = ['_DISC_M2A123-DIDC_M2A123-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) reinst_plot_cfg.imageName = ['_DISC_M2A123-DIDC_M2A123-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end

if strcmp(cfg.cond1, 'DISC_EM2UV1') & strcmp(cfg.cond2,'DIDC_EM2UV1')
    reinst_plot_cfg.lbls3plot = {['DISC'] ['DIDC'] ['DISC vs.' newline 'DIDC']};
    reinst_plot_cfg.imageName = ['_DISC_EM2UV1-DIDC_EM2UV1-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) reinst_plot_cfg.imageName = ['_DISC_EM2UV1-DIDC_EM2UV1-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'SISC_EM2UV2') & strcmp(cfg.cond2,'DISC_EM2UV2')
    reinst_plot_cfg.lbls3plot = {['SISC'] ['DISC'] ['SISC vs.' newline 'DISC']};
    reinst_plot_cfg.imageName = ['_SISC_EM2UV2-DISC_EM2UV2-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) reinst_plot_cfg.imageName = ['_SISC_EM2UV2-_DISC_EM2UV2-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end
if strcmp(cfg.cond1, 'DISC_EM2UV2') & strcmp(cfg.cond2,'DIDC_EM2UV2')
    reinst_plot_cfg.lbls3plot = {['DISC'] ['DIDC'] ['DISC vs.' newline 'DIDC']};
    reinst_plot_cfg.imageName = ['_DISC_EM2UV2-DIDC_EM2UV2-' num2str(cfg.subj2exc)  '.png']; 
    if (cfg.runperm) reinst_plot_cfg.imageName = ['_DISC_EM2UV2-DIDC_EM2UV2-' num2str(cfg.subj2exc) '.png_c.png'] ; end
end




cfg_plot = reinst_plot_cfg;
 
 
end

 

