function  surface_parcellation_concat_cifti_wateredge(cohortfile,outputdir,tmasktype,hems,edges)
%Generate gradient-based parcellation on surface registered subject data
% This script requires a parcellation cohortfile file. The cohortfile includes
% a list of the subject names, the cifti timeseries to be parcellated including
% its path, and the directory where the subject's surface is (as generated by 
% Freesurfer and registered to fs_LR space), formatted as:
% 
% subjectname funcpath/cifti_file surfdir
% 
% e.g. 
%
% vc33416 /data/cn4/laumannt/vc33416__BOLD_LR_surf_subcort_smooth2.55_32k_fsLR.dtseries.nii /data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/vc33416/7112b_fs_LR
%
% This function will concatenate all input cifti timeseries into a single
% cifti timeseries on which all subsequent operations are performed.
% With 'tmasktype', users can enter a tmasklist.txt as created by
% fcimage_analysis to mask out timepoints in their data. Enter 'none' if no
% temporal mask is needed.
% If 'smooth' is set to  'yes' data will be smoothed along the cortical
% surface after surface projection. Enter 'no' if no surface smoothing is
% desired.
% With 'hems', users can specify that either the 'LEFT', 'RIGHT', or 'BOTH'
% hemispheres be parcellated.
% If 'edges' is set to 'yes' a final non-maxima suppression step will be
% applied to the gradient maps and an edge frequency map will be output in
% addition to an average gradient map. Enter 'no' if edge detection is not
% desired.
% This version of surface parcellation will always average correlation data
% from all subjects, generating a single gradient and/or edge map for the
% group.
% TOL 08/06/13

workbenchdir = 'env nice -n 6 /data/cn4/laumannt/workbench/bin_linux64/';
roidir = '/data/cn4/laumannt/subcortical_mask';
HEMS = {'L';'R'};
hemname = {'LEFT';'RIGHT'};
hemname_low = {'left';'right'};

[ign subnum] = system(['cat ' cohortfile ' | wc -l']);
subnum = str2num(subnum);

% Read in subject names, functional volume locations, and surface directory
for s = 1:subnum
[ign subjects{s,1}] = system(['awk -F '' '' ''NR==' num2str(s) '{printf $1}'' ' cohortfile]);
[ign cifti_files{s,1}] = system(['awk -F '' '' ''NR==' num2str(s) '{printf $2}'' ' cohortfile]);
%[ign funcpaths{s,1}] = system(['awk -F '' '' ''NR==' num2str(s) '{printf $2}'' ' cohortfile ' | awk -F ''/'' ''sub(FS $NF,x)''']);
%[ign cifti_files{s,1}] = system(['awk -F '' '' ''NR==' num2str(s) '{printf $2}'' ' cohortfile ' | awk -F ''/'' ''{print $NF}'' | awk -F ''.'' ''{print $1}''']);
%[ign cifti_files{s,1}] = system(['awk -F '' '' ''NR==' num2str(s) '{printf $2}'' ' cohortfile ' | awk -F ''/'' ''{print $NF}''']);
%[ign surfdirs{s,1}] = system(['awk -F '' '' ''NR==' num2str(s) '{printf $3}'' ' cohortfile]);
end

% Check masktype and the existence of custom masks
switch tmasktype
    case 'none'
    otherwise
        [tmasksubjects tmaskfiles]=textread(tmasktype,'%s%s');
        if ~isequal(tmasksubjects,subjects)
            error('masklist subjects do not match cohortfile subjects');
        end
        for i=1:numel(tmaskfiles)
            fprintf('\t%d\t%s\n',i,subjects{i,1});
            if ~exist(tmaskfiles{i,1})
                error('\t%s is not found\n',tmaskfiles{i,1});
            end
        end
end

switch hems
    case 'LEFT'
        h = 1;
    case 'RIGHT'
        h = 2;
    case 'BOTH'
        h = [1:2];
end

% Make output folder
system(['mkdir ' outputdir])
 
% Concatenate data
system(['rm ' outputdir '/allsubs_total_tmask.txt'])
system(['touch ' outputdir '/allsubs_total_tmask.txt'])
time_concat = 'allsubs_LR_timeseries.dtseries.nii';
mergestring = [workbenchdir '/wb_command -cifti-merge ' outputdir '/' time_concat];
for s = 1:length(subjects)
    disp(['Processing subject #' num2str(s) ': ' subjects{s}])
    subject = strcat(subjects{s});
    
%     %create cifti with eroded medial wall
%     volume = ['/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/FCPROCESS_bandpass_interp_nosmooth/' subject '/' subject '_333_zmdt_resid_ntrpl_bpss_zmdt_sd1.nii'];
%     cifti_file = [outputdir '/' subject '_32k_fsLR_eroded.dtseries.nii'];
%     create_cifti_from_volume_timeseries(volume,cifti_file);

    %cifti_file = [funcpaths{s} '/' cifti_files{s}];
    cifti_file = cifti_files{s};
    
    % Concatenate tmasks
    switch tmasktype
        case 'none'
        otherwise
        system(['cat ' tmaskfiles{s} ' >> ' outputdir '/allsubs_total_tmask.txt'])
    end
    mergestring = strcat(mergestring,[' -cifti ' cifti_file]);
  
end
disp('Concatenating timeseries')
system(mergestring)    

for hem = h   
    
    concat_corr = ['allsubs_concat_corr_' HEMS{hem} '.dconn.nii'];
    
    % Calculate correlation maps of concatenated timeseries
    disp('Calculating correlation maps')
    switch tmasktype
        case 'none'
            system([workbenchdir '/wb_command -cifti-correlation ' outputdir '/' time_concat ' ' outputdir '/' concat_corr ' -roi-override -' hemname_low{hem} '-roi ' roidir '/' HEMS{hem} '.atlasroi_erode3.32k_fs_LR.shape.gii -fisher-z'])
 
        otherwise
            system([workbenchdir '/wb_command -cifti-correlation ' outputdir '/' time_concat ' ' outputdir '/' concat_corr ' -roi-override -' hemname_low{hem} '-roi ' roidir '/' HEMS{hem} '.atlasroi_erode3.32k_fs_LR.shape.gii -weights ' outputdir '/allsubs_total_tmask.txt -fisher-z'])
            
    end
    disp('done')
    
    disp('Removing NaNs')
    system([workbenchdir '/wb_command -cifti-math "a" ' outputdir '/allsubs_concat_corr_temp.dconn.nii -fixnan 0 -var a ' outputdir '/' concat_corr])
    disp('done')
    system(['mv ' outputdir '/allsubs_concat_corr_temp.dconn.nii ' outputdir '/' concat_corr])
        
    % Calculate corr of corr
    disp('Calculating corr of corr')
    system([workbenchdir '/wb_command -cifti-correlation ' outputdir '/' concat_corr ' ' outputdir '/avg_corrofcorr_' HEMS{hem} '.dconn.nii'])
    disp('done')
 
    % Calculate gradients       
    atlasdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser';
    midsurf_32k = [atlasdir '/fsaverage_LR32k/Conte69.' HEMS{hem} '.midthickness.32k_fs_LR.surf.gii'];

    disp('Calculating gradient')
    gradsname = ['avgcorrofcorr_allgrad_' HEMS{hem}];  
    system([workbenchdir '/wb_command -cifti-gradient ' outputdir '/avg_corrofcorr_' HEMS{hem} '.dconn.nii COLUMN ' outputdir '/' gradsname '.dconn.nii -' hemname_low{hem} '-surface ' midsurf_32k ' -surface-presmooth 2.55']);
    
    % Average gradients
    system([workbenchdir '/wb_command -cifti-reduce ' outputdir '/' gradsname '.dconn.nii MEAN ' outputdir '/' gradsname '_avg.dconn.nii'])

    % Smooth gradients before edge detection
    disp('Smoothing gradient')
    smooth = 2.55;
    system([workbenchdir '/wb_command -cifti-smoothing ' outputdir '/' gradsname '.dconn.nii ' num2str(smooth) ' 0 COLUMN ' outputdir '/' gradsname '_smooth' num2str(smooth) '.dconn.nii -' hemname_low{hem} '-surface ' midsurf_32k]);
   
    % Average smooth gradients
    system([workbenchdir '/wb_command -cifti-reduce ' outputdir '/' gradsname '_smooth' num2str(smooth) '.dconn.nii MEAN ' outputdir '/' gradsname '_smooth' num2str(smooth) '_avg.dconn.nii'])

    % Load medial wall mask
    medial_wall = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall_erode3.' HEMS{hem} '.32k_fs_LR.func.gii']);   
    medial_wall = medial_wall.cdata;  
    
    % Convert cifti to load into matlab
    system([workbenchdir '/wb_command -cifti-convert -to-gifti-ext ' outputdir '/' gradsname '_smooth' num2str(smooth) '.dconn.nii ' outputdir '/' gradsname '_smooth' num2str(smooth) '.func.gii']);
    
    grads = gifti([outputdir '/' gradsname '_smooth' num2str(smooth) '.func.gii']);
    fullgrads = zeros(32492,nnz(medial_wall==0));
    fullgrads(~logical(medial_wall),:) = grads.cdata; 
    clear grads
       
    % Calculate edge maps
    switch edges
        case 'no'
        case 'yes'   
            disp('Calculating edges')
            fullgrads_medial = fullgrads;
            fullgrads_medial(logical(medial_wall),:) = 1000;
            minimametrics = metric_minima_all(fullgrads_medial,3); 
            clear fullgrads_medial
            labels = watershed_algorithm_all_par(fullgrads,minimametrics,200,1);
            labels_avg = mean(labels==0,2);
            save(gifti(single(labels_avg)),[outputdir '/' gradsname '_smooth' num2str(smooth) '_wateredge_avg.func.gii']);
            
            replabels_avg = repmat(labels_avg,1,nnz(medial_wall));
            replabels_avg = replabels_avg > .2;
            sumedges = sum((alledges .* replabels_avg),1);
            ignore_verts_thresh = mean(sumedges) - std(sumedges);
            goodverts_label_avg = mean(labels(:,sumedges>ignore_verts_thresh)==0,2);
            save(gifti(single(goodverts_label_avg)),[outputdir '/' gradsname '_smooth' num2str(smooth) '_wateredge_goodverts_avg.func.gii']);
            
    end
    %system(['rm ' outputdir '/' concat_corr])
    system(['rm ' outputdir '/avg_corrofcorr_' HEMS{hem} '.dconn.nii'])
    system(['rm ' outputdir '/' gradsname '_smooth' num2str(smooth) '.func.*']);
    system(['rm ' outputdir '/' gradsname '.dconn.nii']);
    system(['rm ' outputdir '/' gradsname '_smooth' num2str(smooth) '.dconn.nii'])
end
%hems
%exit