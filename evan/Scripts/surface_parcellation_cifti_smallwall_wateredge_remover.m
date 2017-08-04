function  surface_parcellation_cifti_smallwall_wateredge_remover(avgcorrfilename,outputdir,hems,edges)

%function  surface_parcellation_cifti_smallwall_wateredge_remover(cohortfile,outputdir,tmasktype,hems,edges)
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

% smoothnum = 2.55;
% mask4dfp='/data/cn4/laumannt/Standard/glm_atlas_mask_333.4dfp.img';
 workbenchdir = 'env nice -n 4 /data/cn4/laumannt/workbench/bin_linux64/';
% %roidir = '/data/cn4/laumannt/subcortical_mask';
 roidir = '/data/hcp-zfs/home/laumannt/120_parcellation';
 HEMS = {'L';'R'};
 hemname = {'LEFT';'RIGHT'};
 hemname_low = {'left';'right'};
% 
% [ign subnum] = system(['cat ' cohortfile ' | wc -l']);
% subnum = str2num(subnum);
% 
% % Read in subject names, functional volume locations, and surface directory
% for s = 1:subnum
% [ign subjects{s,1}] = system(['awk -F '' '' ''NR==' num2str(s) '{printf $1}'' ' cohortfile]);
% [ign funcpaths{s,1}] = system(['awk -F '' '' ''NR==' num2str(s) '{printf $2}'' ' cohortfile ' | awk -F ''/'' ''sub(FS $NF,x)''']);
% [ign cifti_files{s,1}] = system(['awk -F '' '' ''NR==' num2str(s) '{printf $2}'' ' cohortfile ' | awk -F ''/'' ''{print $NF}'' | awk -F ''.'' ''{print $1"."$2}''']);
% %[ign cifti_files{s,1}] = system(['awk -F '' '' ''NR==' num2str(s) '{printf $2}'' ' cohortfile ' | awk -F ''/'' ''{print $NF}''']);
% [ign surfdirs{s,1}] = system(['awk -F '' '' ''NR==' num2str(s) '{printf $3}'' ' cohortfile]);
% end
% 
% % Check masktype and the existence of custom masks
% switch tmasktype
%     case 'none'
%     otherwise
%         [tmasksubjects tmaskfiles]=textread(tmasktype,'%s%s');
%         if ~isequal(tmasksubjects,subjects)
%             error('masklist subjects do not match cohortfile subjects');
%         end
%         for i=1:numel(tmaskfiles)
%             fprintf('\t%d\t%s\n',i,subjects{i,1});
%             if ~exist(tmaskfiles{i,1})
%                 error('\t%s is not found\n',tmaskfiles{i,1});
%             end
%         end
% end
% 
switch hems
    case 'LEFT'
        h = 1;
    case 'RIGHT'
        h = 2;
%     case 'BOTH'
%         h = [1:2];
end
% 
% % Make output folder
% system(['mkdir ' outputdir])

for hem = h
    
    
%     for s = 1:length(subjects)
%         disp(['Processing subject #' num2str(s) ': ' subjects{s}])
%         subject = strcat(subjects{s});
%         funcpath = strcat(funcpaths{s});
%         cifti_file = strcat(cifti_files{s});
%         surfdir = strcat(surfdirs{s});     
%         
%         % Tmask data
%         disp('Calculating correlation maps')
%         switch tmasktype
%             case 'none'
%             system([workbenchdir '/wb_command -cifti-correlation ' funcpath '/' cifti_file '.dtseries.nii ' outputdir '/' cifti_file '_corr.dconn.nii -roi-override -' hemname_low{hem} '-roi ' roidir '/' HEMS{hem} '.atlasroi_group_noproj.func.gii -fisher-z'])
%             otherwise
%             system([workbenchdir '/wb_command -cifti-correlation ' funcpath '/' cifti_file '.dtseries.nii ' outputdir '/' cifti_file '_corr.dconn.nii -roi-override -' hemname_low{hem} '-roi ' roidir '/' HEMS{hem} '.atlasroi_group_noproj.func.gii -weights ' tmaskfiles{s} ' -fisher-z'])
%             
%         end
%         disp('done')
%       
%         disp('Removing NaNs')
%         system([workbenchdir '/wb_command -cifti-math "a" ' outputdir '/' cifti_file '_corr_temp.dconn.nii -fixnan 0 -var a ' outputdir '/' cifti_file '_corr.dconn.nii'])
%         disp('done')
%         system(['mv ' outputdir '/' cifti_file '_corr_temp.dconn.nii ' outputdir '/' cifti_file '_corr.dconn.nii'])
%         
%         %Initialize average correlation matrix on sub 1, then add matrices
%         %for subsequent subs/runs
%         if  (s == 1)
%             disp('Initializing average correlation matrix')
%             system(['mv ' outputdir '/' cifti_file '_corr.dconn.nii ' outputdir '/avg_corr_' HEMS{hem} '.dconn.nii'])
%         else
%             disp('Adding to running average')
%             system([workbenchdir '/wb_command -cifti-math "a+b" ' outputdir '/avg_temp.dconn.nii -var a ' outputdir '/avg_corr_' HEMS{hem} '.dconn.nii -var b ' outputdir '/' cifti_file '_corr.dconn.nii']);
%             disp('done')
%             system(['mv ' outputdir '/avg_temp.dconn.nii ' outputdir '/avg_corr_' HEMS{hem} '.dconn.nii'])
%             system(['rm ' outputdir '/' cifti_file '_corr.dconn.nii'])
%         end
%     end   
%     %Take average
%     if numel(subjects)>1
%         disp('Calculating average correlation map')
%         system([workbenchdir '/wb_command -cifti-math "a/' num2str(length(subjects)) '" ' outputdir '/avg_temp.dconn.nii -var a ' outputdir '/avg_corr_' HEMS{hem} '.dconn.nii '])
%         disp('done')
%         system(['mv ' outputdir '/avg_temp.dconn.nii ' outputdir '/avg_corr_' HEMS{hem} '.dconn.nii'])
%     end
       
% mask = gifti([roidir '/' HEMS{hem} '.atlasroi_group_noproj.func.gii']); mask = mask.cdata;
% maskL = gifti([roidir '/L.atlasroi_group_noproj.func.gii']); maskL = maskL.cdata;
% 
% disp('Selecting verts for this hemisphere')
% 
% system([workbenchdir '/wb_command -cifti-convert -to-gifti-ext ' avgcorrfilename ' ' outputdir '/Temp.func.gii']); 
% avgcorr = gifti([outputdir '/Temp.func.gii']); avgcorr = avgcorr.cdata;
% avgcorr = avgcorr([1:nnz(mask)] + (nnz(maskL)*strcmp('R',HEMS{hem})) , :);
% save(gifti(single(avgcorr)),[outputdir '/Temp.func.gii'],'ExternalFileBinary')
% system([workbenchdir '/wb_command -cifti-convert -from-gifti-ext ' outputdir '/Temp.func.gii ' outputdir '/avg_corr_' HEMS{hem} '.dconn.nii']); 
% delete([outputdir '/Temp.func.gii*'])


    %Calculate corr of corr
    disp('Calculating corr of corr')
    %system([workbenchdir '/wb_command -cifti-correlation ' outputdir '/avg_corr_' HEMS{hem} '.dconn.nii ' outputdir '/avg_corrofcorr_' HEMS{hem} '.dconn.nii'])
    system([workbenchdir '/wb_command -cifti-correlation ' avgcorrfilename ' ' outputdir '/avg_corrofcorr_' HEMS{hem} '.dconn.nii'])
    disp('done')
 
    % Calculate gradients       
    atlasdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser';
    midsurf_32k = [atlasdir '/fsaverage_LR32k/Conte69.' HEMS{hem} '.midthickness.32k_fs_LR.surf.gii'];

    disp('Calculating gradient')
      gradsname = ['avgcorrofcorr_allgrad_' HEMS{hem}];  
    system([workbenchdir '/wb_command -cifti-gradient ' outputdir '/avg_corrofcorr_' HEMS{hem} '.dconn.nii COLUMN ' outputdir '/' gradsname '.dconn.nii -' hemname_low{hem} '-surface ' midsurf_32k ' -surface-presmooth 2.55']);
    
    % Average gradients
    system([workbenchdir '/wb_command -cifti-reduce ' outputdir '/' gradsname '.dconn.nii MEAN ' outputdir '/' gradsname '_avg.dtseries.nii'])

    % Smooth gradients before edge detection
    disp('Smoothing gradient')
      smooth = 2.55;
    system([workbenchdir '/wb_command -cifti-smoothing ' outputdir '/' gradsname '.dconn.nii ' num2str(smooth) ' 0 COLUMN ' outputdir '/' gradsname '_smooth' num2str(smooth) '.dconn.nii -' hemname_low{hem} '-surface ' midsurf_32k]);
   
    % Average smooth gradients
    system([workbenchdir '/wb_command -cifti-reduce ' outputdir '/' gradsname '_smooth' num2str(smooth) '.dconn.nii MEAN ' outputdir '/' gradsname '_smooth' num2str(smooth) '_avg.dtseries.nii'])
    
    % Load medial wall mask
    medial_wall = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall_erode3.' HEMS{hem} '.32k_fs_LR.func.gii']);   
    %medial_wall = gifti(['/data/hcp-zfs/home/laumannt/120_parcellation/medial_exclude_' HEMS{hem} '.func.gii']);
    medial_wall = medial_wall.cdata;  
    
    % Convert cifti to load into matlab
    system([workbenchdir '/wb_command -cifti-convert -to-gifti-ext ' outputdir '/' gradsname '_smooth' num2str(smooth) '.dconn.nii ' outputdir '/' gradsname '_smooth' num2str(smooth) '.func.gii']);

    grads = gifti([outputdir '/avgcorrofcorr_allgrad_' HEMS{hem} '_smooth' num2str(smooth) '.func.gii']);
    fullgrads = zeros(32492,nnz(medial_wall==0));
    fullgrads(~logical(medial_wall),:) = grads.cdata; 
    clear grads
    
    % For nonmax edge calculation
    atlasdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser';
    specfile = [atlasdir '/fsaverage_LR32k/Conte69.' HEMS{hem} '.32k_fs_LR.c5.spec'];
    
    % Calculate edge maps
    switch edges
        case 'no'
        case 'yes'   
            disp('Calculating edges')
            fullgrads_medial = fullgrads;
            fullgrads_medial(logical(medial_wall),:) = 1000;
            minimametrics = metric_minima_all(fullgrads_medial,3); 
            clear fullgrads_medial
            %Original water edge
            labels = watershed_algorithm_all_par(fullgrads,minimametrics,200,1);
            labels_avg = mean(labels==0,2);
            save(gifti(single(labels_avg)),[outputdir '/' gradsname '_smooth' num2str(smooth) '_wateredge_avg.func.gii']);
            
            save([outputdir 'labels_' HEMS{hem}],'labels','-v7.3');
            
            %Water edge with below threshold edges removed
            labels = watershed_edgeremover(labels,fullgrads,0.01);
            labels_avg = mean(labels,2);
            save(gifti(single(labels_avg)),[outputdir '/' gradsname '_smooth' num2str(smooth) '_wateredge_thresh0.01_avg.func.gii']);
            
            %Nonmaxima suppression-based edges
%            surface_edges_all_test_faster(fullgrads,specfile,outputdir,[gradsname '_smooth' num2str(smooth) '_edge_avg'],1);
    end
    
    system(['rm ' outputdir '/avg_corrofcorr_' HEMS{hem} '.dconn.nii'])
   % system(['rm ' outputdir '/' gradsname '_smooth' num2str(smooth) '.func.*']);
    system(['rm ' outputdir '/' gradsname '.dconn.nii']);
    system(['rm ' outputdir '/' gradsname '_smooth' num2str(smooth) '.dconn.nii'])

end
%exit