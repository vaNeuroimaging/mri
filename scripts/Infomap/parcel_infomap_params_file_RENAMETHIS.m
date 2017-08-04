function params = parcel_infomap_params_file
%
%Parameters to set for parcel_infomap.m

%Parcels being evaluated
parcelfile = '/data/cn4/evan/Published_parcels/Parcels_LR.dtseries.nii';

%Folder results will be written to
outputfolder = '/data/cn4/evan/Temp/partials/partial_corr/';



%Density thresholds to run
thresholdarray = [.005 : .001 : .05];

%Geodesic distance exclusion, in mm (i.e. parcels cannnot be "connected" if
%they are closer than this distance)
xdistance = 30;

%Smallest collection of parcels that is allowed to exist as a separate
%community
networksizeminimum = 5;



%Do you want to calculate the average parcel-to-parcel correlation matrix?
%1 = yes; 0 = no (i.e. it already exists and will be loaded)
calc_corrmat = 0;

%Do you want to calculate the average parcel-to-parcel geodesic distance matrix?
%1 = yes; 0 = no (i.e. it already exists and will be loaded)
calc_distances = 0;

%Name of correlation matrix to be calculated and saved or loaded
corrmatname = 'corrmat.mat';

%Name of distance matrix to be calculated and saved or loaded
parceldistancesname = 'parcel_distances.mat';


%List of data files and tmask files. 
%ONLY USED if calc_corrmat = 1;
cohortfile = '/data/cn5/selfRegulation/V4Process_nosmooth/Surfdatalist_120_108.txt'; %Should be formatted as (IDnumber)  (cifti timeseries location)
tmaskfile = '/data/cn5/selfRegulation/V4Process_nosmooth/Finaltmasklist_120_108.txt'; %Should be formatted as (IDnumber)  (tmask location)

%Matrix containing vertex-to-vertex geodesic distances.
%ONLY USED if calc_distances = 1;
distancesfile = '/data/cn/data1/scripts/CIFTI_RELATED/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Cifti_geo_distances_xhemisphere_euclidean_uint8.mat';


%-------------------------------------------------------------------------
%DO NOT CHANGE

%load variables into an output structure
varnames = who;
for i = 1:length(varnames)
    evalc(['params.' varnames{i} ' = ' varnames{i} ';']);
end