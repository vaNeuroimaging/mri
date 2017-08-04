function params = post_fc_processing_batch_params
%
%Parameters to set for post_fc_processing_batch.m

%--------------------------------------------------------------------------

%data list, in the same format required by FCPROCESS
datalist = '/data/cn4/dgreene/Patients/AdultTS/NoFieldMap/COHORTSELECT_CONTROLS/NEW_DATALIST.txt';

%tmask list, in the same format required by FCPROCESS
tmasklist = '/data/cn4/dgreene/Patients/AdultTS/NoFieldMap/COHORTSELECT_CONTROLS/NEW_TMASKLIST.txt';

%Location of final UNSMOOTHED fc-processed data
fcprocessed_funcdata_dir = '/data/cn4/dgreene/Patients/AdultTS/NoFieldMap/ATS_controls_initial/';

%Generalied suffix of final UNSMOOTHED fc-processed data, such that the
%data files are named [subject number][suffix].4dfp.img
fcprocessed_funcdata_dir_suffix = '_333_zmdt_resid_ntrpl_bpss_zmdt';

%Location cifti files will be written to. This folder will be created, and
%sub-folders will be created within this location. 
outfolder = '/data/cn4/evan/Evan_brain/';

%Location of the volumetric subcortical mask label file. Data within this
%mask will be smoothed and included in the cifti files. This should contain
%different labels for each volumetric structure identified by freesurfer.
subcort_mask = '/data/cn/data1/scripts/CIFTI_RELATED/Resources/cifti_masks/subcortical_mask_LR_333.nii';

%Location of subjects' fs_LR-registered surfaces. Underneath this folder,
%subject data should be in '[subjectnumber]/7112b_fs_LR/' ; this folder
%should contain 'Native' and 'fsaverage_LR32k' subfolders with surfaces.
fs_LR_surfdir = '/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/';

%Location of atlas medial wall masks (where the cortical surfaces don't
%have cortical data). Data from these regions will not be included in the
%cifti files.
medial_mask_L = '/data/cn/data1/scripts/CIFTI_RELATED/Resources/cifti_masks/L.atlasroi.32k_fs_LR.shape.gii';
medial_mask_R = '/data/cn/data1/scripts/CIFTI_RELATED/Resources/cifti_masks/R.atlasroi.32k_fs_LR.shape.gii';

%Additional cifti files will be produced using the specified "smallwall"
%medial wall masks, which are highly eroded versions of the atlas medial
%wall masks. These are dataset-specific, and are usually built as the
%cortical locations where at least one subject in the dataset has no data
%projected. Ciftis made using "smallwall" masks are primarily used for
%parcellation. Leave empty ([]) to omit generation of ciftis with smallwall
%masks.
sw_medial_mask_L = '/data/pruett/CPD/ImagingStudies/BabySibs_fcMRI//Zeran/voxel_community/mode_49_V24_wmparc/L.atlasroi_noproj.func.gii';
sw_medial_mask_R = '/data/pruett/CPD/ImagingStudies/BabySibs_fcMRI//Zeran/voxel_community/mode_49_V24_wmparc/R.atlasroi_noproj.func.gii';

%Sigma of the smoothing kernel to be applied geodesically to surface data and
%volumetrically to the subcortical data
smoothnum = 2.55;



%-------------------------------------------------------------------------
%DO NOT CHANGE

%load variables into an output structure
varnames = who;
for i = 1:length(varnames)
    evalc(['params.' varnames{i} ' = ' varnames{i} ';']);
end

