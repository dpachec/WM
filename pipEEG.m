%%
clear
subjID = 's10';
path = '/Users/danielpacheco/Desktop/electrode_MNI/china/';
%path = '/Users/danielpacheco/Desktop/jings_data/';

cd ([path subjID]);

%% convert to nii
subjID = 's07';
dcmSource = [path subjID '/postCT'];
%niiFolder = [path subjID '/anat/' subjID '_T1w.nii.gz'];
niiFolder = [path subjID '/anat/' subjID '_CT.nii.gz'];
%niiFolder = [path subjID '/anat/again/' subjID '_CT.nii.gz'];

dicm2nii(dcmSource, niiFolder)


%% preprocessing of the anatomical MRI

mri = ft_read_mri([path subjID '/anat/' subjID '_T1w.nii.gz']);


%%
ft_determine_coordsys(mri, 'interactive', 'no')

% left-right orientation



%%
cfg = [];
cfg.method = 'interactive';
cfg.coordsys = 'acpc';
mri_acpc = ft_volumerealign(cfg, mri);

%% Write the preprocessed anatomical MRI out to a file as shown below.

cfg = [];
cfg.filename = [path '/' subjID '/anat/' subjID '_MR_acpc'];
cfg.filetype = 'nifti';
cfg.parameter = 'anatomy';
ft_volumewrite(cfg, mri_acpc);




%% Cortical surface extraction with FreeSurfer

fshome = '/Applications/freesurfer/7.1.1';
subdir = ['/Users/danielpacheco/Desktop/jings_data/' subjID '/anat'];
mrfile = ['/Users/danielpacheco/Desktop/jings_data/' subjID '/anat/' subjID '_MR_acpc.nii'];
system(['export FREESURFER_HOME=' fshome '; ' ...
'source $FREESURFER_HOME/SetUpFreeSurfer.sh; ' ...
'mri_convert -c -oc 0 0 0 ' mrfile ' ' [subdir '/tmp.nii'] '; ' ...2
'recon-all -i ' [subdir '/tmp.nii'] ' -s ' 'freesurfer' ' -sd ' subdir ' -all']);




%%
clear
subjID = 's10';
path = '/Users/danielpacheco/Desktop/china/';
%path = '/Users/danielpacheco/Desktop/jings_data/';

cd ([path subjID]);

%% Import the extracted cortical surfaces into the MATLAB workspace 
%%and examine their quality.

pial_lh = ft_read_headshape([path subjID '/freesurfer/surf/lh.pial']);
pial_lh.coordsys = 'acpc';
ft_plot_mesh(pial_lh);
lighting gouraud;
camlight;


%% Import the FreeSurfer-processed MRI
fsmri_acpc = ft_read_mri([path subjID '/freesurfer/mri/T1.mgz']);
%fsmri_acpc = ft_read_mri([path subjID '/freesurfer/mri/T1.mgz']);
fsmri_acpc.coordsys = 'acpc';



%% Import the anatomical CT
ct = ft_read_mri([path subjID '/anat/' subjID '_CT.nii.gz']);




%% determine the native orientation of the anatomical CT’s
ft_determine_coordsys(ct, 'interactive', 'no')

%% Align the anatomical CT to the CTF head surface coordinate system

cfg = [];
cfg.method = 'interactive';
cfg.coordsys = 'acpc'; % use acpc if not possible with ctf
ct_ctf = ft_volumerealign(cfg, ct);

%% convert the CT’s coordinate system into an approximation of the ACPC coordinate system
%Fuse the CT with the MRI using the below command. ~ aprox 40s
%Write the MRI-fused anatomical CT out to a file 

tic

ct_acpc = ft_convert_coordsys(ct_ctf, 'acpc');

%%Fuse the CT with the MRI using the below command. ~ aprox 40s
cfg = [];
cfg.method = 'spm';
cfg.spmversion = 'spm12';
cfg.coordsys = 'acpc';
cfg.viewresult = 'yes';
ct_acpc_f = ft_volumerealign(cfg, ct_acpc, fsmri_acpc);

%% Write the MRI-fused anatomical CT out to a file 
cfg = [];
cfg.filename = [subjID '_CT_acpc_f'];
cfg.filetype = 'nifti';
cfg.parameter = 'anatomy';
ft_volumewrite(cfg, ct_acpc_f);

toc

%% Import the header information from the recording file

hdr = ft_read_header([path subjID '/' subjID '.edf']);


%% Localize the electrodes in the post-implant CT 
%mri = ft_read_mri([path subjID '/anat/' subjID '_T1w.nii.gz']);
%ct_acpc_f = ft_read_mri([path subjID '/anat/' subjID '_CT_acpc_f.nii']);
%fsmri_acpc = ft_read_mri([path subjID '/freesurfer/mri/T1.mgz']);

ct_acpc_f  = ft_read_mri(['/Users/danielpacheco/Documents/iEEG_data_analysis/extinction/data/iEEG/electrode_localization/c_sub01/c_sub01_ct_acpc_f.nii']);
fsmri_acpc = ft_read_mri(['/Users/danielpacheco/Documents/iEEG_data_analysis/extinction/data/iEEG/electrode_localization/c_sub01/c_sub01_T1_fs.nii']);

fsmri_acpc.coordsys = 'acpc';
%hdr = ft_read_header([path subjID '/' subjID '.edf']);

cfg = [];
cfg.channel = hdr.label;
%cfg.channel = string(1:64);
elec_acpc_f = ft_electrodeplacement(cfg, ct_acpc_f, fsmri_acpc);



%% to edit electrode struct
% elec_acpc_f.label(18:25) = [];
% elec_acpc_f.elecpos(18:25,:) = [];
% elec_acpc_f.chanpos(18:25,:) = [];

%%

ft_plot_ortho(fsmri_acpc.anatomy, 'transform', fsmri_acpc.transform, 'style', 'intersect');
ft_plot_sens(elec_acpc_f, 'label', 'on', 'fontcolor', 'w');


%%

save([subjID '_elec_acpc_f.mat'], 'elec_acpc_f');



%% Register the subject’s brain to the standard MNI brain ~aprox 1:30min
tic

cfg = [];
cfg.nonlinear = 'yes';
cfg.spmversion = 'spm12';
cfg.spmmethod  = 'new';
fsmri_mni = ft_volumenormalise(cfg, fsmri_acpc);

toc

%% Use the resulting deformation parameters to obtain the electrode positions in standard MNI space
% aprox 20s

tic
elec_mni_frv = elec_acpc_f;
elec_mni_frv.elecpos = ft_warp_apply(fsmri_mni.params, elec_acpc_f.elecpos, 'individual2sn');
elec_mni_frv.chanpos = ft_warp_apply(fsmri_mni.params, elec_acpc_f.chanpos, 'individual2sn');
elec_mni_frv.coordsys = 'mni';
toc

%% Visualize the cortical mesh extracted from the standard MNI brain along with the spatially normalized electrode
pial_lh = ft_read_headshape([path subjID '/freesurfer/surf/lh.pial']);
pial_lh.coordsys = 'acpc';
ft_plot_mesh(pial_lh); 
ft_plot_sens(elec_mni_frv);
view([-55 10]);
material dull; 
lighting gouraud; 
camlight;

%% save normalized electrode to file

save([subjID '_elec_mni_frv.mat'], 'elec_mni_frv');
disp('done')

%% plot average + mni electrodes
tic

%atlas = ft_read_atlas([path 'average_MNI_brain/fsaverage6/mri/aseg.mgz']);
atlas = ft_read_atlas(['/Users/danielpacheco/Documents/iEEG_data_analysis/WM/electrode_MNI/china/average_MNI_brain/fsaverage6/mri/aseg.mgz']);
atlas.coordsys = 'mni';
cfg            = [];
cfg.inputcoord = 'mni';
cfg.atlas      = atlas;
%cfg.roi        = {'Left-Cerebral-Cortex'};
cfg.roi        = {'Right-Hippocampus'};
mask_rha = ft_volumelookup(cfg, atlas);

seg             = keepfields(atlas, {'dim', 'unit','coordsys','transform'});
seg.brain       = mask_rha;
cfg             = [];
cfg.method      = 'iso2mesh';
cfg.radbound    = 2;
cfg.maxsurf     = 0;
cfg.tissue      = 'brain';
cfg.numvertices = 1000;
cfg.smooth      = 3;
cfg.spmversion  = 'spm12';
mesh_rha_average = ft_prepare_mesh(cfg, seg);

%% plot subject specific 2 compare

atlas = ft_read_atlas([path subjID '/freesurfer/mri/aparc+aseg.mgz']);
atlas.coordsys = 'acpc';
cfg            = [];
cfg.inputcoord = 'acpc';
cfg.atlas      = atlas;
cfg.roi        = {'Right-Hippocampus'};
mask_rha     = ft_volumelookup(cfg, atlas);

seg = keepfields(atlas, {'dim', 'unit','coordsys','transform'});
seg.brain = mask_rha;
cfg             = [];
cfg.method      = 'iso2mesh';
cfg.radbound    = 2;
cfg.maxsurf     = 0;
cfg.tissue      = 'brain';
cfg.numvertices = 1000;
cfg.smooth      = 3;
cfg.spmversion  = 'spm12';
mesh_rha_subject = ft_prepare_mesh(cfg, seg);

toc

%% final figure

figure(1)
ft_plot_mesh(mesh_rha_average,  'facealpha', .5);
ft_plot_sens(elec_mni_frv);
title('Average hippocampus + subject MNI');

%%
figure(2)
ft_plot_mesh(mesh_rha_subject,  'facealpha', .5);
ft_plot_sens(elec_acpc_f);
title('Subject specific hippocampus');



%%

[vertices faces] = readObj('left_hippocampus.obj');

pL = patch('Faces',faces,'Vertices',vertices); hold on;
pL.LineStyle = 'none';      % remove the lines
l = light('Position',[-0.4 0.2 0.9],'Style','infinite')
material([.9 .7 .3]) %sets the ambient/diffuse/specular strength of the objects.
view(90,0)

%%
%cd 'D:\owncube\miniXIM\_WM\WM_datasets\china\WM01_m_20190826_hangzhou\caidabao_preMR'
projectdir = 'D:\owncube\miniXIM\_WM\WM_datasets\china\WM01_m_20190826_hangzhou\caidabao_preMR';
dicomFiles = dir( fullfile(projectdir, '*.dcm' ));
files = {dicomFiles.name};
y = length(dicomFiles)
X = zeros(512, 512, 1, y, 'uint8');
% Read the series of images.
for p=1:y
   filename = fullfile( projectdir, dicomFiles(p).name );
   X(:,:,1,p) = dicomread(filename);
end
% Display the image stack.
montage(X,[])


%%
projectdir = 'D:\owncube\miniXIM\_WM\WM_datasets\china\WM01_m_20190826_hangzhou\caidabao_preMR';
dicomFiles = dir( fullfile(projectdir, '*.dcm' ));
files = {dicomFiles.name};




dicom2nifti(files)



%% get all labels again

subjID = 'ASJ';
path = '/Users/danielpacheco/Desktop/agency_mni_electrodes/';


%%

atlas = ft_read_atlas([path subjID '/freesurfer/mri/aparc+aseg.mgz']);
atlas.coordsys = 'acpc';
cfg            = [];
cfg.inputcoord = 'acpc';
cfg.atlas      = atlas;
cfg.roi        = {'Right-Hippocampus'};
mask_rha     = ft_volumelookup(cfg, atlas);



%% 
clear 
myE = readtable('example_2_hippocampus.csv');
elec_mni_frv.unit = 'mm'; 
M = zeros(size(myE,1));
M(1:1+size(M,1):end) = 1;
elec_mni_frv.tra = M;
elec_mni_frv.coordsys = 'mni';

for i = 1:size(myE,1)
   elec_mni_frv.label(i,:) = myE{i, 1};
   elec_mni_frv.elecpos(i,1) = myE{i, 2};
   elec_mni_frv.elecpos(i,2) = myE{i, 3};
   elec_mni_frv.elecpos(i,3) = myE{i, 4};
    
end

elec_mni_frv.chanpos = elec_mni_frv.elecpos;

%% plot average + mni electrodes
tic
path = '/Users/danielpacheco/Desktop/electrode_MNI/china/';
atlas = ft_read_atlas([path 'average_MNI_brain/fsaverage6/mri/aseg.mgz']);
atlas.coordsys = 'mni';
cfg            = [];
cfg.inputcoord = 'mni';
cfg.atlas      = atlas;
%cfg.roi        = {'Left-Cerebral-Cortex'};
cfg.roi        = {'Right-Hippocampus'};
mask_rha = ft_volumelookup(cfg, atlas);

seg             = keepfields(atlas, {'dim', 'unit','coordsys','transform'});
seg.brain       = mask_rha;
cfg             = [];
cfg.method      = 'iso2mesh';
cfg.radbound    = 2;
cfg.maxsurf     = 0;
cfg.tissue      = 'brain';
cfg.numvertices = 1000;
cfg.smooth      = 3;
cfg.spmversion  = 'spm12';
mesh_rha_average = ft_prepare_mesh(cfg, seg);


figure(1)
ft_plot_mesh(mesh_rha_average,  'facealpha', .5);
ft_plot_sens(elec_mni_frv);
title('Average hippocampus + subject MNI');


































