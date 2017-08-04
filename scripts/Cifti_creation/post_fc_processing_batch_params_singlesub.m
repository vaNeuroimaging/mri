function params = post_fc_processing_batch_params_singlesub
%
%Parameters to set for post_fc_processing_batch.m

%--------------------------------------------------------------------------

%Location of the volumetric subcortical mask label file. Data within this
%mask will be smoothed and included in the cifti files. This should contain
%different labels for each volumetric structure identified by freesurfer.
subcort_mask = '/home/data/scripts/Resources/cifti_masks/subcortical_mask_LR_333_MNI.nii.gz';

%Location of atlas medial wall masks (where the cortical surfaces don't
%have cortical data). Data from these regions will not be included in the
%cifti files.
medial_mask_L = '/home/data/scripts/Resources/cifti_masks/L.atlasroi.32k_fs_LR.shape.gii';
medial_mask_R = '/home/data/scripts/Resources/cifti_masks/R.atlasroi.32k_fs_LR.shape.gii';

%Additional cifti files will be produced using the specified "smallwall"
%medial wall masks, which are highly eroded versions of the atlas medial
%wall masks. These are dataset-specific, and are usually built as the
%cortical locations where at least one subject in the dataset has no data
%projected. Ciftis made using "smallwall" masks are primarily used for
%parcellation. Leave empty ([]) to omit generation of ciftis with smallwall
%masks.
%sw_medial_mask_L = '/data/pruett/CPD/ImagingStudies/BabySibs_fcMRI//Zeran/voxel_community/mode_49_V24_wmparc/L.atlasroi_noproj.func.gii';
%sw_medial_mask_R = '/data/pruett/CPD/ImagingStudies/BabySibs_fcMRI//Zeran/voxel_community/mode_49_V24_wmparc/R.atlasroi_noproj.func.gii';

%Sigma of the smoothing kernel to be applied geodesically to surface data and
%volumetrically to the subcortical data




%-------------------------------------------------------------------------
%DO NOT CHANGE

%load variables into an output structure
varnames = who;
for i = 1:length(varnames)
    evalc(['params.' varnames{i} ' = ' varnames{i} ';']);
end

