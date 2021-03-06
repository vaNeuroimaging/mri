function  surface_parcellation_cifti(cohortfile,outputdir,tmasktype,smooth,hems,edges)
%Generate gradient-based parcellation on surface registered subject data
% This script requires a parcellation cohortfile file. The cohortfile includes
% a list of the subject names, the functional volume to be parcellated including
% its path, and the directory where the subject's surface is (as generated by 
% Freesurfer and registered to fs_LR space), formatted as:
% 
% subjectname funcpath/funcvol surfdir
% 
% e.g. 
%
% vc33416 /data/cn4/laumannt/vc33416_rest1.4dfp.img /data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/vc33416/7112b_fs_LR
%
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
% TOL 01/25/13

smoothnum = 2.55;
mask4dfp='/data/cn4/laumannt/Standard/glm_atlas_mask_333.4dfp.img';
workbenchdir = '/data/cn4/laumannt/workbench/bin_linux64/';
roidir = '/data/cn4/laumannt/subcortical_mask';
HEMS = {'L';'R'};
hemname = {'LEFT';'RIGHT'};

[ign subnum] = system(['cat ' cohortfile ' | wc -l']);
subnum = str2num(subnum);

% Read in subject names, functional volume locations, and surface directory
for s = 1:subnum
[ign subjects{s,1}] = system(['awk -F '' '' ''NR==' num2str(s) '{printf $1}'' ' cohortfile]);
[ign funcpaths{s,1}] = system(['awk -F '' '' ''NR==' num2str(s) '{printf $2}'' ' cohortfile ' | awk -F ''/'' ''sub(FS $NF,x)''']);
[ign cifti_files{s,1}] = system(['awk -F '' '' ''NR==' num2str(s) '{printf $2}'' ' cohortfile ' | awk -F ''/'' ''{print $NF}'' | awk -F ''.'' ''{print $1}''']);
[ign surfdirs{s,1}] = system(['awk -F '' '' ''NR==' num2str(s) '{printf $3}'' ' cohortfile]);
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

for hem = h
    
    maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' HEMS{hem} '.32k_fs_LR.func.gii'];
    mask = gifti(maskname);
    mask = mask.cdata;
    
    for s = 1:length(subjects)
        disp(['Processing subject #' num2str(s) ': ' subjects{s}])
        subject = strcat(subjects{s});
        funcpath = strcat(funcpaths{s});
        cifti_file = strcat(cifti_files{s});
        surfdir = strcat(surfdirs{s});     
        
        % Tmask data
        disp('Calculating correlation maps')
        switch tmasktype
            case 'none'
            otherwise
            system([workbenchdir '/wb_command -cifti-correlation ' funcdir '/' cifti_file '.dtseries.nii ' outputdir '/' cifti_file '_corr.dconn.nii -roi-override -' hemname{hem} '-roi ' roidir '/' HEMS{hem} '.atlasroi.32k_fs_LR.shape.gii -weights ' tmaskfiles{s} ' -fisher-z'])
            
        end
        disp('done')
      
        disp('Removing NaNs')
        system([workbenchdir '/wb_command -cifti-math "a" ' outputdir '/' cifti_file '_corr_temp.dconn.nii -fixnan 0 -var a ' outputdir '/' cifti_file '_corr.dconn.nii'])
        disp('done')
        system(['mv ' outputdir '/' cifti_file '_corr_temp.dconn.nii ' outputdir '/' cifti_file '_corr.dconn.nii'])
        
        %Initialize average correlation matrix on sub 1, then add matrices
        %for subsequent subs/runs
        if  (s == 1)
            disp('Initializing average correlation matrix')
            system(['mv ' outputdir '/' cifti_file '_corr.dconn.nii ' outputdir '/avg_corr_' HEMS{hem} '.dconn.nii'])
        else
            disp('Adding to running average')
            system([workbenchdir '/wb_command -cifti-math "a+b" ' outputdir '/avg_temp.dconn.nii -var a ' outputdir '/avg_corr_' HEMS{hem} '.dconn.nii -var b ' outputdir '/' cifti_file '_corr.dconn.nii']);
            disp('done')
            system(['mv ' outputdir '/avg_temp.dconn.nii ' outputdir '/avg_corr_' HEMS{hem} '.dconn.nii'])
            system(['rm ' outputdir '/' cifti_file '_corr.dconn.nii'])
        end
    end   
    %Take average
    if numel(subjects)>1
        disp('Calculating average correlation map')
        system([workbenchdir '/wb_command -cifti-math "a/' num2str(length(subjects)) '" ' outputdir '/avg_temp.dconn.nii -var a ' outputdir '/avg_corr_' HEMS{hem} '.dconn.nii '])
        disp('done')
        system(['mv ' outputdir '/avg_temp.dconn.nii ' outputdir '/avg_corr_' HEMS{hem} '.dconn.nii'])
    end
        
    %Calculate corr of corr
    disp('Calculating corr of corr')
    system([workbenchdir '/wb_command -cifti-correlation ' outputdir '/avg_corr_' HEMS{hem} '.dconn.nii ' outputdir '/avg_corrofcorr_' HEMS{hem} '.dconn.nii'])
    disp('done')
 
    %Calculate gradient       
    atlasdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser';
    specfile = [atlasdir '/fsaverage_LR32k/Conte69.' HEMS{hem} '.32k_fs_LR.c5.spec'];
    midsurf_32k = [atlasdir '/fsaverage_LR32k/Conte69.' HEMS{hem} '.midthickness.32k_fs_LR.surf.gii'];
    disp('Calculating gradient')
    system([workbenchdir '/wb_command -cifti-gradient ' outputdir '/avgcorrofcorr_' HEMS{hem} '.dconn.nii COLUMN ' outputdir '/avgcorrofcorr_allgrad_' HEMS{hem} '.dconn.nii -' hemname{hem} '-surface ' midsurf_32k ' -surface-presmooth 2.55']);
    system([workbenchdir '/wb_command -cifti-convert -to-gifti-ext ' outputdir '/avgcorrofcorr_allgrad_' HEMS{hem} '.dconn.nii ' outputdir '/avgcorrofcorr_allgrad_' HEMS{hem} '.func.gii']);
    
    grads = gifti([outputdir '/avgcorrofcorr_' HEMS{hem} '_allgrad.func.gii']);
    fullgrads = zeros(32492,29696);
    fullgrads(~logical(mask),:) = grads.cdata; 
    clear grads
    meangrad = mean(fullgrads,2);
    save(gifti(meangrad),[outputdir '/avgcorrofcorr_allgrad_' HEMS{hem} '_avg.func.gii'])  
    disp('done')
   
    %Smooth gradient before edge detection
    disp('Smoothing gradient')
    smooth = 2.55;
    system([workbenchdir '/wb_command -cifti-smoothing ' outputdir '/avgcorrofcorr_allgrad_' HEMS{hem} '.dconn.nii ' num2str(smooth) ' 0 COLUMN ' outputdir '/avgcorrofcorr_allgrad_' HEMS{hem} '_smooth' num2str(smooth) '.dconn.nii -' hemname{hem} '-surface ' midsurf_32k]);
    system([workbenchdir '/wb_command -cifti-convert -to-gifti-ext ' outputdir '/avgcorrofcorr_allgrad_' HEMS{hem} '_smooth' num2str(smooth) '.dconn.nii ' outputdir '/avgcorrofcorr_allgrad_' HEMS{hem} '_smooth' num2str(smooth) '.func.gii']);

    grads = gifti([outputdir '/avgcorrofcorr_allgrad_' HEMS{hem} '_smooth' num2str(smooth) '.func.gii']);
    fullgrads = zeros(32492,29696);
    fullgrads(~logical(mask),:) = grads.cdata; 
    clear grads
    meangrad = mean(fullgrads,2);
    save(gifti(meangrad),[outputdir '/avgcorrofcorr_allgrad_' HEMS{hem} '_smooth' num2str(smooth) '_avg.func.gii']);  
        
    % Calculate edge maps
    switch edges
        case 'no'
        case 'yes'   
            disp('Calculating edges')
            surface_edges_all_test_faster(fullgrads,specfile,outputdir,[outputdir '/avgcorrofcorr_allgrad_' HEMS{hem} '_smooth' num2str(smooth) '_edges_avg'],1);
    end
    

end
exit