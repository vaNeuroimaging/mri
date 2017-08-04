function  surface_parcellation_subsurf_nosmooth_allgrad(tmasklist)

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

outputdir = pwd;

workbenchdir = 'nice /data/cn4/evan/workbench/bin_linux64/';
orig_labelfile = '/data/cn4/laumannt/subcortical_mask/mode_subcortical_label.nii';

medial_walls = {'/data/hcp-zfs/home/laumannt/120_parcellation/L.atlasroi_group_noproj.func.gii','/data/hcp-zfs/home/laumannt/120_parcellation/R.atlasroi_group_noproj.func.gii'};

%  system(['nifti_4dfp -4 ' orig_labelfile ' mode_subcortical.4dfp.img'])
%  system(['t4img_4dfp none mode_subcortical.4dfp.img mode_subcortical_' targetres '.4dfp.img -O' targetres])
%  system(['nifti_4dfp -n mode_subcortical_' targetres '.4dfp.img mode_subcortical_' targetres '.nii'])
%  system(['wb_command -volume-label-import mode_subcortical_' targetres '.nii /data/cn4/laumannt/subcortical_mask/FreeSurferSubcorticalLabelTableLut_nobrainstem_LR.txt mode_subcortical_' targetres '_label.nii'])
%  
%  save(gifti(single(zeros(32492,1))),'Temp_surface.func.gii')
%  system(['wb_command -cifti-create-dense-timeseries Template_' targetres '.dtseries.nii -volume mode_subcortical_' targetres '.nii mode_subcortical_' targetres '_label.nii -left-metric Temp_surface.func.gii -roi-left ' medial_walls{1} ' -right-metric Temp_surface.func.gii -roi-right ' medial_walls{1}])
%  system(['wb_command -cifti-convert -to-gifti-ext Template_' targetres '.dtseries.nii Template_forwriting' targetres '.func.gii'])
%  system(['wb_command -cifti-correlation Template_' targetres '.dtseries.nii Template_' targetres '_corr.dconn.nii'])
% % system(['wb_command -cifti-correlation Template_' targetres '_corr.dconn.nii Template_' targetres '_corrofcorr.dconn.nii'])
%  system(['wb_command -cifti-convert -to-gifti-ext Template_' targetres '_corr.dconn.nii Template_' targetres '_corr.func.gii'])
% delete(['Template_' targetres '_corr.dconn.nii']);
% delete(['Template_' targetres '_corr.func.gii.data']);

corrtemplatefile = ['/data/hcp-zfs/home/laumannt/120_parcellation/cifti_correlation_smallwall/60sub_C1_avg_corr_LR.func.gii'];

template_file = ['/data/hcp-zfs/home/laumannt/120_parcellation/cifti_timeseries_smallwall/vc35175_BOLD_LR_surf_subcort_32k_fsLR_smooth2.55.func.gii'];
%subcort_ind = logical(cifti_read(subcort_mask));

global neighbors
neighbors = cifti_neighbors(template_file);

[subjects, tmaskfiles]=textread(tmasklist,'%s%s');

% Make output folder
system(['mkdir ' outputdir])

cd(outputdir)
    
%    templatedata = cifti_read(template_file);
%    fullgrads = zeros(size(templatedata,1));
    
gradsname = ['avg_corrofcorr_allgrad'];
    tic
    for s = 1:length(subjects)
        
        subject = subjects{s};
        %cifti_file = cifti_files{s}(1:(end-13));
        cifti_file = ['/data/hcp-zfs/home/laumannt/120_parcellation/cifti_timeseries_smallwall/' subject '_BOLD_LR_surf_subcort_32k_fsLR_smooth2.55.dtseries.nii'];
        %system(['wb_command -cifti-resample ' cifti_file ' COLUMN ' template_file ' COLUMN BARYCENTRIC TRILINEAR Temp_subdata.dtseries.nii -volume-predilate 4'])
%        surfdir = surfdirs{s};
        tmask = load(tmaskfiles{s});
        surf_L = ['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' subject '/7112b_fs_LR/fsaverage_LR32k/' subject '.L.midthickness.32k_fs_LR.surf.gii'];
        surf_R = ['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' subject '/7112b_fs_LR/fsaverage_LR32k/' subject '.R.midthickness.32k_fs_LR.surf.gii'];
        disp(['Processing subject #' num2str(s) ': ' subject])
        
        % Convert cifti to gifti and load timecourse
        disp('Calculating correlation maps')
%        evalc(['!' workbenchdir '/wb_command -cifti-convert -to-gifti-ext ' cifti_file '.dtseries.nii ' cifti_file '.func.gii']);
        
%        cifti_timecourse = gifti([cifti_file '.func.gii']);
%        cifti_timecourse = cifti_timecourse.cdata;
        
cifti_timecourse = cifti_read(cifti_file);

%    system(['wb_command -cifti-correlation Temp_subdata.dtseries.nii Temp_subdata_corr.dconn.nii -fisher-z -weights ' tmaskfiles{s}])
        
        % Calculate correlation maps
        %cifti_corrmap = corrcoef(cifti_timecourse(cort_ind,logical(tmask))');
         cifti_corrmap = paircorr_mod(cifti_timecourse(:,logical(tmask))');
         disp('done')
        
        
        disp('Removing NaNs')
%        system([workbenchdir '/wb_command -cifti-math "a" Temp_subdata_corr_nonan.dconn.nii -fixnan 0 -var a Temp_subdata_corr.dconn.nii'])
%        delete('Temp_subdata_corr.dconn.nii')
        toc
         cifti_corrmap(isnan(cifti_corrmap)) = 0;
         cifti_corrmap = FisherTransform(cifti_corrmap);
        
        % Calculate corr of corr
        disp('Calculating corr of corr')
%        system(['wb_command -cifti-correlation Temp_subdata_corr_nonan.dconn.nii Temp_subdata_corrofcorr.dconn.nii -fisher-z'])
%        disp('Removing NaNs')
%        system([workbenchdir '/wb_command -cifti-math "a" Temp_subdata_corrofcorr_nonan.dconn.nii -fixnan 0 -var a Temp_subdata_corrofcorr.dconn.nii'])
%        delete('Temp_subdata_corrofcorr.dconn.nii')
%        delete('Temp_subdata_corr_nonan.dconn.nii')
         corrofcorr = corrcoef(cifti_corrmap);
         clear cifti_corrmap
         corrofcorr(isnan(corrofcorr)) = 0;
         corrofcorr = FisherTransform(corrofcorr);
        disp('done')
        
        toc
        
        % Write out cifti file
        
        cifti_write_wHDR(corrofcorr,corrtemplatefile,['Temp_subdata_corrofcorr'],'dconn')
        clear corrofcorr
        toc
        % Calculate gradients
        disp('Calculating gradient')
%        gradsname = [];
        evalc(['!' workbenchdir '/wb_command -cifti-gradient Temp_subdata_corrofcorr.dconn.nii COLUMN Temp_subdata_grads.dconn.nii -left-surface ' surf_L ' -right-surface ' surf_R])%' -surface-presmooth 2.55']);
        delete('Temp_subdata_corrofcorr.dconn.nii')
        toc
        if s==1
            copyfile('Temp_subdata_grads.dconn.nii',[gradsname '_runningavg.dconn.nii'])
        else
            system([workbenchdir '/wb_command -cifti-math "a+b" ' gradsname '_runningavg.dconn.nii -var a Temp_subdata_grads.dconn.nii -var b ' gradsname '_runningavg.dconn.nii']);
        end
        delete('Temp_subdata_grads.dconn.nii');
%         
%         %         % Average gradients
%         %         system([workbenchdir '/wb_command -cifti-reduce ' outputdir '/' gradsname '.dconn.nii MEAN ' outputdir '/' gradsname '_avg.dtseries.nii'])
%         %
%         
%         % Convert gradients and load
%         evalc(['!' workbenchdir '/wb_command -cifti-convert -to-gifti-ext ' outputdir '/' gradsname '.dconn.nii ' outputdir '/' gradsname '.func.gii'])
%         grads = gifti([outputdir '/' gradsname '.func.gii']);
%         toc
%         
%         % Add subject gradients to running average
%         fullgrads = [fullgrads + grads.cdata];
        toc
    end
    
    
    % Average gradients across subjects
    %fullgrads = fullgrads./length(subjects);
    
    disp('Calculating average correlation map')
        system([workbenchdir '/wb_command -cifti-math "a/' num2str(length(subjects)) '" ' gradsname '.dconn.nii -var a ' gradsname '_runningavg.dconn.nii'])
        delete([gradsname '_runningavg.dconn.nii'])
        disp('done')
    
    % Save out average gradients
    %templatefile = ['/data/cn4/evan/ROIs/subcort_corrofcorr_template.func.gii'];
    
    
    %cifti_write_wHDR(fullgrads,corrtemplatefile,gradsname,'dconn')
    %clear fullgrads
    
    % Smooth gradients before edge detection
    atlasdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser'; % Use atlas surface now that data is averaged
    midsurf_32k_L = [atlasdir '/fsaverage_LR32k/Conte69.L.midthickness.32k_fs_LR.surf.gii'];
    midsurf_32k_R = [atlasdir '/fsaverage_LR32k/Conte69.R.midthickness.32k_fs_LR.surf.gii'];
    disp('Smoothing gradient')
    smooth = 2.55;
    system([workbenchdir '/wb_command -cifti-smoothing ' outputdir '/' gradsname '.dconn.nii ' num2str(smooth) ' ' num2str(smooth) ' COLUMN ' outputdir '/' gradsname '_smooth' num2str(smooth) '.dconn.nii -left-surface ' midsurf_32k_L ' -right-surface ' midsurf_32k_R]);
%     
    % Average smooth gradients
    system([workbenchdir '/wb_command -cifti-reduce ' outputdir '/' gradsname '_smooth' num2str(smooth) '.dconn.nii MEAN ' outputdir '/' gradsname '_smooth' num2str(smooth) 'avg.dtseries.nii'])
    
    % Convert gradients and load
    %evalc(['!' workbenchdir '/wb_command -cifti-convert -to-gifti-ext ' outputdir '/' gradsname '_smooth' num2str(smooth) '.dconn.nii ' outputdir '/' gradsname '_smooth' num2str(smooth) '.func.gii'])
    %grads = gifti([outputdir '/' gradsname '_smooth' num2str(smooth) '.func.gii']); grads = grads.cdata;
    grads = cifti_read([outputdir '/' gradsname '_smooth' num2str(smooth) '.dconn.nii']);
    
    
    % Calculate edge maps
    
            disp('Calculating edges')
            minimametrics = metric_minima_all_cifti(grads,3);
            clear fullgrads_medial
            %Original water edge
            labels = watershed_algorithm_all_par_cifti(grads,minimametrics,200,1);
            labels_avg = mean(labels==0,2);
            %save(gifti(single(labels_avg)),[outputdir '/' gradsname '_smooth' num2str(smooth) '_wateredge_avg.func.gii']);
            %cifti_write_wHDR(single(labels_avg),'/data/cn4/evan/ROIs/subcort_corr_template.func.gii',[outputdir '/' gradsname '_smooth' num2str(smooth) '_wateredge_avg'])
%             output = zeros(59412 + size(labels_avg),1);
%             output(59412+1:end) = labels_avg;

            cifti_write_wHDR(labels_avg,template_file,[outputdir '/' gradsname '_smooth' num2str(smooth) '_wateredge_avg'])
            
            save([outputdir '/labels_subcort'],'labels','-v7.3');
            
%             % Water edge with below threshold edges removed
%             thresh = 0.01;%:.005:.04;
%             for t = 1:length(thresh)
%                 newlabels = watershed_edgeremover(labels,fullgrads_smooth,thresh(t));
%                 labels_avg = mean(newlabels,2);
%                 save(gifti(single(labels_avg)),[outputdir '/' gradsname '_smooth' num2str(smooth) '_wateredge_thresh' num2str(thresh(t)) '_avg.func.gii']);
%             end
%             
%             % Nonmaxima suppression-based edges
%             surface_edges_all_test_faster(fullgrads_smooth,specfile,outputdir,[gradsname '_edge_avg'],1);
   
    
    % system(['rm ' outputdir '/avg_corrofcorr_' HEMS{hem} '.dconn.nii'])
    % system(['rm ' outputdir '/' gradsname '_smooth' num2str(smooth) '.func.*']);
    % system(['rm ' outputdir '/' gradsname '.dconn.nii']);
    % system(['rm ' outputdir '/' gradsname '_smooth' num2str(smooth) '.dconn.nii'])
    

end

function minimametric = metric_minima_all_cifti(metric,neighdist)

global neighbors
         

minimametric = zeros(size(metric));
for i = 1:length(metric)
    
    nodeneigh = i;
    newneigh = nodeneigh;
    curneigh = newneigh;
    
    for n = 1:neighdist
        
        for t = 1:length(curneigh)
            newneigh = [newneigh neighbors(curneigh(t),2:end)];
            newneigh(isnan(newneigh)) = [];
        end
        curneigh = setdiff(newneigh,nodeneigh);
        nodeneigh = union(nodeneigh,newneigh,'stable');
        
    end
    
    nodeneighval_all = metric(nodeneigh(2:end),:);
    origval_all = repmat(metric(nodeneigh(1),:),[length(nodeneigh)-1 1]);
    
    minval = sum(origval_all<nodeneighval_all,1);
    %[minval mini] = min(nodeneighval_all);
    minimametric(i,:) = single(minval==(length(nodeneigh)-1));
    
end
end

function labels = watershed_algorithm_all_par_cifti(edgemetrics,minimametrics,stepnum,fracmaxh)

global neighbors

%Label initial markers with unique value
labels = zeros(size(minimametrics));
    
minh = min(edgemetrics(:));
maxh = max(edgemetrics(:));

stoph = maxh*fracmaxh;
step = (maxh-minh)/stepnum;
hiter = minh:step:stoph;
system('mkdir countingDir')
system('touch countingDir/countfile')
matlabpool open 8
numlabels = size(labels,2);
for j = 1:2
    
    if j == 1
        labelstodo = 1:floor(numlabels/2);
    elseif j == 2
        labelstodo = (floor(numlabels/2)+1):numlabels;
    end
    parfor l = labelstodo%size(labels,2)
        system(['echo ' num2str(l) ' >> countingDir/countfile']);
        [ign count] = system('cat countingDir/countfile | wc -l');
        disp(['Calculating edges on metric #:' num2str(l) ', Count = ' num2str(count)])
        label = labels(:,l);
        edgemetric = edgemetrics(:,l);
        
        labelpos = find(minimametrics(:,l)==1);
        randval = randn(length(labelpos),1);
        [randign randind] = sort(randval);
        temp = 1:length(labelpos);
        labelnums = temp(randind);
        label(labelpos) = labelnums;
        watershed_zones = zeros(size(label));
        for i = 1:length(hiter);
            %disp(['Number of iterations will be ' num2str(length(hiter)) ', Iteration = ' num2str(i)])
            %maskpos = find(edgemetric<sortedge(i)); % Take values in metric less than current iteration
            maskmetrics = edgemetric<hiter(i); % Take values in metric less than current iteration
            maskmetrics = maskmetrics & ~label>0 & ~watershed_zones;
            
            maskpos = find(sum(maskmetrics,2)>0);
            
            randnum = randperm(length(maskpos));
            maskpos = maskpos(randnum);
            
            for m = 1:length(maskpos) %For all nodes at this threshold
                
                nodeneigh = neighbors(maskpos(m),2:end);
                maskinthismetric = maskmetrics(maskpos(m),:);
                %nodeneigh = neighbors(sortedgepos(i),2:7);
                
                nodeneigh(isnan(nodeneigh)) = [];
                nodeneighlab = label(nodeneigh,:);
                
                %Find minimum value other than 0 among neighbors
                minfindnodeneighlab = nodeneighlab;
                minfindnodeneighlab(nodeneighlab==0) = 100000;
                minnodeneighlab = min(minfindnodeneighlab,[],1);
                
                %Find maximum value other than 0 among neighbors
                maxfindnodeneighlab = nodeneighlab;
                maxfindnodeneighlab(nodeneighlab==0) = -100000;
                maxnodeneighlab = max(maxfindnodeneighlab,[],1);
                
                %If min and max differ (i.e. two or more neighbor water), watershed
                %zone
                watershed_nodes = (minnodeneighlab~=maxnodeneighlab) & (minnodeneighlab~=100000) & maskinthismetric;
                watershed_zones(maskpos(m),watershed_nodes) = 1;
                label(maskpos(m),watershed_nodes) = 0;
                
                %If min and max the same but different from 0, add to neighbor
                %water
                next_to_water = (minnodeneighlab==maxnodeneighlab) & (minnodeneighlab~=100000) & maskinthismetric;
                label(maskpos(m),next_to_water) = minnodeneighlab(next_to_water);
                
                
            end
            
        end
        labels(:,l) = label;
    end
end
system('rm -r countingDir')
matlabpool close

%   
%  save(gifti(label),[outputdir '/' filestem 'watershed.gii'])
%  system(['mv ' outputdir '/' filestem 'watershed.gii ' outputdir '/' filestem 'watershed.func.gii'])
end