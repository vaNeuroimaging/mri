function params = template_matching_params_file
%Parameters to set for Template_matching_cortexonly.m


%List of subject timecourses in cifti format
surfdatafile = '/data/cn/data1/scripts/CIFTI_RELATED/Template_Matching/Surfdatalist_120_108.txt'; %should be formatted as (subject number)  (cifti timecourse location)


%List of subject tmasks
tmaskfile = '/data/cn/data1/scripts/CIFTI_RELATED/Template_Matching/Finaltmasklist_120_108.txt'; %should be formatted as (subject number)  (tmask location)


%Exclusion distance, in mm. Connectivity map overlap will not be assessed
%within this distance of a seed vertex
xdistance = 30;


%Minimum size of contiguous system patches. Patches smaller than this will
%be removed and filled in with neighboring system identities
minclustersizemm2 = 30;


%Density threshold for template system connectivity maps and subject
%vertex-seeded connectivity maps
kden_thresh = .05;


%Location of template maps. Should contain a 'templates' variable that is
%(#vertices) X (#systems) and contains unthresholded connectivity maps for
%each system. Should also contain 'IDs', a 1 X (#systems) vector of numeric
%IDs for each system.
templatesfile = '/data/cn/data1/scripts/CIFTI_RELATED/Template_Matching/Templates_consensus.mat'; 
%'/data/cn4/evan/RestingState/Ind_variability/Templates_Yeo.mat'


%Matrix containing vertex-to-vertex geodesic distances.
distancesfile = '/data/cn/data1/scripts/CIFTI_RELATED/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Cifti_geo_distances_xhemisphere_euclidean_uint8.mat';


%Location of subject fs_LR surfaces (will be in /subjectname/7112b_fs_LR/fsaverage_LR32k/ within this folder)
surffolder = '/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/';

%-------------------------------------------------------------------------
%DO NOT CHANGE

%load variables into an output structure
varnames = who;
for i = 1:length(varnames)
    evalc(['params.' varnames{i} ' = ' varnames{i} ';']);
end
