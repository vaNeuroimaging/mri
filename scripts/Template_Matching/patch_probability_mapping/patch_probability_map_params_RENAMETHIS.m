function params = patch_probability_map_params
%params = patch_probability_map_params
%
% params to be specified for Patch_probability_map.m

%Location of a file containing subject-specific system maps (one subject
%per column)
subnetworksfile = ['Templatematch_dice_bysubject_kden0.05.dscalar.nii'];

%Location of a file describing the surface area of each cortical vertex
surfaceareasfile = '/data/cn/data1/scripts/CIFTI_RELATED/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.LR.midthickness.32k_fs_LR_surfaceareas.dtseries.nii';

%Maximum distance metric (mean closest point-to-point distance) beyond
%which two subjects' patches are not allowed to be matched.
distthresh = 10;

%Minimum subjects who must have a cluster of patches for it to survive in
%the final thresholded output.
pctsubthresh = 1/3;

%Subjects' system maps must be clustered into contiguous pieces before
%matching. Set to 1 if this needs to be done; set to 0 if this has already
%happened and you want to re-use the previous clustering (saves a lot of time). 
new_cluster = 0;

%Name of file containing data from subjects' discrete contiguous system
%pieces. Will be created or loaded.
system_clustering_names = 'System_clustering.mat';

%Labels of systems. The index of each name in this variable should match
%the numeric value of the system in the subjects' maps.
networklabels = {'Default','Visual','FrontoPar','FrontoPar 2','DorsalAttn','DorsalAttn2','VentAttn','Salience','CingOperc','MortorHand','MotorMouth','Auditory','MTL1','MTL2','MedPar','ParOccip'}; %Standard Power colors
%networklabels = {'Default','Visual','FrontoPar','FrontoPar 2','DorsalAttn','DorsalAttn2','VentAttn','Salience','CingOperc','MortorHand','MotorMouth','Auditory','MTL1','STS','MedPar','ParOccip','Orbitofrontal','LatVisual'};

%A file describing vertex-to-vertex geodesic distances. Use a version
%with large cross-hemispheric distances (this prevents an unlikely scenario
%of matching clusters across the medial wall). Using a uint8 version is
%recommended as it will be faster to load, and the script reduces it to
%uint8 anyway to save time.
distances_file = '/data/cn/data1/scripts/CIFTI_RELATED/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Cifti_geo_distances_xhemisphere_large_uint8.mat';

%Folder where results will be written to.
outputfolder = '/data/cn4/evan/RestingState/SystemClusters/';

%-------------------------------------------------------------------------
%DO NOT CHANGE

%load variables into an output structure
varnames = who;
for i = 1:length(varnames)
    evalc(['params.' varnames{i} ' = ' varnames{i} ';']);
end