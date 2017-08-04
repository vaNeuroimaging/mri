%Parameters to set for post_fc_processing_batch.m

%data list, in the same format required by FCPROCESS
datalist = '/data/cn4/dgreene/Patients/AdultTS/NoFieldMap/COHORTSELECT_CONTROLS/NEW_DATALIST.txt';

%tmask list, in the same format required by FCPROCESS
tmasklist = '/data/cn4/dgreene/Patients/AdultTS/NoFieldMap/COHORTSELECT_CONTROLS/NEW_TMASKLIST.txt';

%Location of final fc-processed data
fcprocessed_funcdata_dir = '/data/cn4/dgreene/Patients/AdultTS/NoFieldMap/ATS_controls_initial/';

%Location cifti files will be written to. Sub-folders will be created
%within this location.
outfolder = '/data/cn4/evan/Evan_brain/';

%Location of the volumetric subcortical mask label file. Data within this
%mask will be smoothed and included in the cifti files. This should contain
%different labels for each volumetric structure identified by freesurfer.
subcort_mask = '/data/cn4/laumannt/subcortical_mask/subcortical_mask_LR_333.nii';

%Location of subjects' fs_LR-registered surfaces. Underneath this folder,
%subject data should be in '[subjectnumber]/7112b_fs_LR/' ; this folder
%shoudl contain 'Native' and 'fsaverage_LR32k' subfolders with surfaces.
fs_LR_surfdir = '/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/';

%Location of atlas medial wall masks (where the cortical surfaces don't
%have cortical data). Data from these regions will not be included in the
%cifti files.
medial_mask_L = '/data/cn4/laumannt/subcortical_mask/L.atlasroi.32k_fs_LR.shape.gii';
medial_mask_R = '/data/cn4/laumannt/subcortical_mask/R.atlasroi.32k_fs_LR.shape.gii';

%Sigma of the smoothing kernel to be applied geodesically to surface data and
%volumetrically to the subcortical data
smoothnum = 2.55;
