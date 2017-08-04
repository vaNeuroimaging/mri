function  surface_parcellation_lowRAM_singlesub(subject,cifti_file,surfdir,tmaskfile,subsample,dosubcort,outputdir)
%surface_parcellation_lowRAM_singlesub(subject,cifti_file,surfdir,tmaskfile,subsample,dosubcort,outputdir)
%
% Generate gradient-based parcellation on a single subject's surface
% registered data in cifti format. 
%
% 'subject' is a string containing the subject's ID
%
% 'cifti_file' is a string containing the path to the subject's cifti
% timecourse
%
% 'surfdir' is a string containing the path to the subject's 32k fs_LR
% directory
%
% 'tmaskfile' is a string containing the path to the subject's tmask file,
% a text file of ones and zeros indicating which timepoints to exclude
%
% 'subsample' indicates the subsampling of maps to run the gradient on.
% Omit or leave empty ([]) to use the default of subsample=100, which means
% that the gradient will be run on 1/100th of similarity maps, randomly
% selected. The primary purpose of this is to greatly speed up the
% parcellation, which can take weeks otherwise. Testing indicates that
% subsampling every 100th map produces final parcellations which correlate
% at r>.99 with non-subsampled maps. 
%
% 'dosubcort' should be a binary (1 or zero) indicating whether subcortical
% structures will also be parcellated. Parcellating subcortical structures
% approximately doubles total parcellation time. Note that subcortical
% parcellations have not been validated. Omit or leave empty ([]) to use
% dosubcort=1;
%
% 'outputdir' specifies the folder to write results into. Omit or leave
% empty ([]) to write into the working directory. 
%
%
% Requires the Cifti Resources scripts to be in your path (e.g.,
% /data/cn/data1/scripts/CIFTI_RELATED/Resources/ and subfolders)
%
% TOL 01/25/13, modified by EMG 03/3/16


% PARAMETERS TO SET
smooth = 2.55; % sigma for geodesic smoothing applied to gradient maps
workbenchdir = 'nice /data/cn/data1/scripts/CIFTI_RELATED/Resources/workbench/bin_linux64/'; %location of workbench
divisions = 10;
%-----------------------------------------------------------------------


if ~exist('subsample') || isempty(subsample)
    subsample = 100;
end

if ~exist('dosubcort') || isempty(dosubcort)
    dosubcort = 1;
end

if ~exist('outputdir') || isempty(outputdir)
    outputdir = pwd;
end



% Make output folder
mkdir(outputdir)
cd(outputdir)



tmask = load(tmaskfile);

%Name of this subject's midthickness surface
midsurf_32k_sub = {[surfdir '/' subject '.L.midthickness.32k_fs_LR.surf.gii'],[surfdir '/' subject '.R.midthickness.32k_fs_LR.surf.gii']};


disp(['Processing subject ' subject])

% Load cifti file
ciftistruct = ft_read_cifti_mod(cifti_file);
cifti_timecourse = single(ciftistruct.data);
ciftistruct.data = [];
ciftistruct_orig = ciftistruct;

% Set up some values to be used across subjects

num_samples = round(size(cifti_timecourse,1)/subsample);

randinds = randperm(size(cifti_timecourse,1));
randinds = randinds(1:num_samples);

% Calculate correlation maps
disp('Calculating correlation map')
%cifti_corrmap = paircorr_mod(cifti_timecourse(:,logical(tmask))');

subsampled_corrmap = paircorr_mod(cifti_timecourse(:,logical(tmask))',cifti_timecourse(randinds,logical(tmask))');
subsampled_corrmap(isnan(subsampled_corrmap)) = 0;
subsampled_corrmap = FisherTransform(subsampled_corrmap);


% Remove NaNs (produced if vertices have no data)
%cifti_corrmap(isnan(cifti_corrmap)) = 0;
% Apply the Fisher tranformation
%cifti_corrmap = single(FisherTransform(cifti_corrmap));

corrofcorr = zeros(size(cifti_timecourse,1),num_samples,'single');

divisionsize = ceil(size(corrofcorr,1)/divisions);

% Calculate correlation similarity
disp('Calculating similarity map')
for division = 1:divisions
    divisioninds = [((division-1) * divisionsize +1) : min([size(corrofcorr,1), (division * divisionsize)])];
    
    division_corrmap = paircorr_mod(cifti_timecourse(:,logical(tmask))',cifti_timecourse(divisioninds,logical(tmask))');
    division_corrmap(isnan(division_corrmap)) = 0;
    division_corrmap = FisherTransform(division_corrmap);
    
    corrofcorr(divisioninds,:) = paircorr_mod(division_corrmap,subsampled_corrmap);
    
    clear division_corrmap
    
end
clear cifti_timecourse division_corrmap subsampled_corrmap

% Remove NaNs
corrofcorr(isnan(corrofcorr)) = 0;
% Apply the Fisher tranformation
corrofcorr = FisherTransform(corrofcorr);

% Write out corr of corr cifti file for gradient calculation
ciftistruct_orig.data = corrofcorr;
if ~dosubcort
    ncortverts = nnz((ciftistruct_orig.brainstructure > 0) & (ciftistruct_orig.brainstructure < 3));
    ciftistruct_orig.data = ciftistruct_orig.data(1:ncortverts,:);
    ciftistruct_orig.brainstructure(ciftistruct_orig.brainstructure > 2) = [];
    ciftistruct_orig.pos(ciftistruct_orig.brainstructure > 2) = [];
    ciftistruct_orig.brainstructurelabel(3:end) = [];
    ciftistruct_orig = rmfield(ciftistruct_orig,'dim');
    ciftistruct_orig = rmfield(ciftistruct_orig,'transform');
end
clear corrofcorr
ft_write_cifti_mod([outputdir '/corrofcorr_LR_subcort'],ciftistruct_orig);
ciftistruct_orig.data = [];


% Calculate gradients
disp('Calculating gradient')
gradsname = [subject '_corrofcorr_allgrad_LR_subcort'];
system([ workbenchdir '/wb_command -cifti-gradient ' outputdir '/corrofcorr_LR_subcort.dtseries.nii COLUMN ' outputdir '/' gradsname '.dtseries.nii -left-surface ' midsurf_32k_sub{1} ' -right-surface ' midsurf_32k_sub{2}]);


delete([outputdir '/corrofcorr_LR_subcort.dtseries.nii'])



% Smooth gradients before edge detection
disp('Smoothing gradient')
system([workbenchdir '/wb_command -cifti-smoothing ' outputdir '/' gradsname '.dtseries.nii ' num2str(smooth) ' ' num2str(smooth) ' COLUMN ' outputdir '/' gradsname '_smooth' num2str(smooth) '.dtseries.nii -left-surface ' midsurf_32k_sub{1} ' -right-surface ' midsurf_32k_sub{2}]);

neighbors = cifti_neighbors([outputdir '/' gradsname '_smooth' num2str(smooth) '.dtseries.nii']);

% Load smoothed gradients
ciftistruct = ft_read_cifti_mod([outputdir '/' gradsname '_smooth' num2str(smooth) '.dtseries.nii']);

%Average smoothed gradients across maps (unused but useful to look at)
ciftistruct_orig.data = mean(ciftistruct.data,2);
ft_write_cifti_mod([outputdir '/' gradsname '_smooth' num2str(smooth) 'avg'],ciftistruct_orig);

fullgrads_smooth = ciftistruct.data;
ciftistruct.data = [];

delete([outputdir '/' gradsname '.dtseries.nii'])

disp('Calculating edges')


% Get local minima of each smoothed gradient map
minimametrics = metric_minima_all_cifti(fullgrads_smooth,3,neighbors);

% Run watershed-by-flooding algorithm on each gradient map to generate edges
edges = watershed_algorithm_all_par_cifti(fullgrads_smooth,minimametrics,200,1,neighbors);
clear fullgrads_smooth

% Average across gradient maps and save
edge_density = mean(edges==0,2);
ciftistruct.data = edge_density;
ft_write_cifti_mod([outputdir '/' gradsname '_smooth' num2str(smooth) '_wateredge_avg'],ciftistruct);

% Save watershed edges from all gradient maps
ciftistruct.data = edges;
ft_write_cifti_mod([outputdir '/' gradsname '_smooth' num2str(smooth) '_wateredge_all'],ciftistruct);
clear edges ciftistruct_orig

