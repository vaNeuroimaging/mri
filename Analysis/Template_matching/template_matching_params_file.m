function params = template_matching_params_file
%Parameters to set for Template_matching_cortexonly.m


%List of subject timecourses in cifti format
surfdatafile = '/home/data/Analysis/connectome_datalist_10min.txt'; %should be formatted as (subject number)  (cifti timecourse location)


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
templatesfile = '/home/data/scripts/Template_Matching/Templates_consensus.mat'; 


%Matrix containing vertex-to-vertex geodesic distances.
distancesfile = '/home/data/scripts/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.LR.32k_fs_LR.distances_surfgeo_voleuc_normalwall_standardsubcort_uint8.mat';


%-------------------------------------------------------------------------
%DO NOT CHANGE

%load variables into an output structure
varnames = who;
for i = 1:length(varnames)
    evalc(['params.' varnames{i} ' = ' varnames{i} ';']);
end
