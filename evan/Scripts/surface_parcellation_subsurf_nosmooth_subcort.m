function  surface_parcellation_subsurf_nosmooth_subcort(tmasklist,outputdir)
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
%workbenchdir = 'env nice -n 4 /data/cn4/laumannt/workbench/bin_linux64/';
workbenchdir = 'nice /data/cn4/evan/workbench/bin_linux64/';
%roidir = '/data/cn4/laumannt/subcortical_mask';
%roidir = '/data/hcp-zfs/home/laumannt/120_parcellation';
roidir = '/data/cn4/laumannt/subcortical_mask';
HEMS = {'L';'R'};
hemname = {'LEFT';'RIGHT'};
hemname_low = {'left';'right'};

templatefile = ['/data/cn4/evan/ROIs/subcort_corrofcorr_template.func.gii'];

subcort_mask = '/data/cn4/evan/ROIs/cifti_subcort_mask.dtseries.nii';
subcort_ind = logical(cifti_read(subcort_mask));

[subjects tmaskfiles]=textread(tmasklist,'%s%s');

% Make output folder
system(['mkdir ' outputdir])

cd(outputdir)
    
    
    fullgrads = zeros(nnz(subcort_ind));
    

    tic
    for s = 1:length(subjects)
        
        subject = subjects{s};
        %cifti_file = cifti_files{s}(1:(end-13));
        cifti_file = ['/data/hcp-zfs/home/laumannt/120_parcellation/cifti_timeseries_normalwall/' subject '_BOLD_LR_surf_subcort_32k_fsLR_smooth2.55.dtseries.nii'];
%        surfdir = surfdirs{s};
        tmask = load(tmaskfiles{s});
        surfdir = ['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' subject '/7112b_fs_LR/'];
        disp(['Processing subject #' num2str(s) ': ' subject])
        
        % Convert cifti to gifti and load timecourse
        disp('Calculating correlation maps')
%        evalc(['!' workbenchdir '/wb_command -cifti-convert -to-gifti-ext ' cifti_file '.dtseries.nii ' cifti_file '.func.gii']);
        
%        cifti_timecourse = gifti([cifti_file '.func.gii']);
%        cifti_timecourse = cifti_timecourse.cdata;
        cifti_timecourse = cifti_read(cifti_file);
        
        % Calculate correlation maps
        %cifti_corrmap = corrcoef(cifti_timecourse(cort_ind,logical(tmask))');
        cifti_corrmap = paircorr_mod(cifti_timecourse(subcort_ind,logical(tmask))',cifti_timecourse(:,logical(tmask))');
        disp('done')
        toc
        
        disp('Removing NaNs')
        cifti_corrmap(isnan(cifti_corrmap)) = 0;
        
        % Calculate corr of corr
        disp('Calculating corr of corr')
        corrofcorr = corrcoef(cifti_corrmap');
        disp('done')
        toc
        
        % Write out cifti file
        
        cifti_write_wHDR(corrofcorr,templatefile,['corrofcorr_subcort'],'dconn')
        toc
        % Calculate gradients
        disp('Calculating gradient')
        gradsname = ['corrofcorr_allgrad_subcort'];
        evalc(['!' workbenchdir '/wb_command -cifti-gradient ' outputdir '/corrofcorr_subcort.dconn.nii COLUMN ' outputdir '/' gradsname '.dconn.nii'])%' -surface-presmooth 2.55']);
        
        %         % Average gradients
        %         system([workbenchdir '/wb_command -cifti-reduce ' outputdir '/' gradsname '.dconn.nii MEAN ' outputdir '/' gradsname '_avg.dtseries.nii'])
        %
        
        % Convert gradients and load
        evalc(['!' workbenchdir '/wb_command -cifti-convert -to-gifti-ext ' outputdir '/' gradsname '.dconn.nii ' outputdir '/' gradsname '.func.gii'])
        grads = gifti([outputdir '/' gradsname '.func.gii']);
        toc
        
        % Add subject gradients to running average
        fullgrads = [fullgrads + grads.cdata];
        toc
    end
    
    
    % Average gradients across subjects
    fullgrads = fullgrads./length(subjects);
    
    % Save out average gradients
    templatefile = ['/data/cn4/evan/ROIs/subcort_corrofcorr_template.func.gii'];
    
    gradsname = ['avg_corrofcorr_allgrad_subcort'];
    cifti_write_wHDR(fullgrads,templatefile,gradsname,'dconn')
    clear fullgrads
    
    % Smooth gradients before edge detection
    atlasdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser'; % Use atlas surface now that data is averaged
    %midsurf_32k = [atlasdir '/fsaverage_LR32k/Conte69.' HEMS{hem} '.midthickness.32k_fs_LR.surf.gii'];
    disp('Smoothing gradient')
      smooth = 2.55;
     system([workbenchdir '/wb_command -cifti-smoothing ' outputdir '/' gradsname '.dconn.nii 0 ' num2str(smooth) ' COLUMN ' outputdir '/' gradsname '_smooth' num2str(smooth) '.dconn.nii']);
%     
    % Average smooth gradients
    system([workbenchdir '/wb_command -cifti-reduce ' outputdir '/' gradsname '_smooth' num2str(smooth) '.dconn.nii MEAN ' outputdir '/' gradsname '_smooth' num2str(smooth) 'avg.dtseries.nii'])
    
    % Convert gradients and load
    evalc(['!' workbenchdir '/wb_command -cifti-convert -to-gifti-ext ' outputdir '/' gradsname '_smooth' num2str(smooth) '.dconn.nii ' outputdir '/' gradsname '_smooth' num2str(smooth) '.func.gii'])
    grads = gifti([outputdir '/' gradsname '_smooth' num2str(smooth) '.func.gii']); grads = grads.cdata;
    
    
    % Calculate edge maps
    
            disp('Calculating edges')
            minimametrics = metric_minima_all_subcort(grads,3);
            clear fullgrads_medial
            %Original water edge
            labels = watershed_algorithm_all_par_subcort(grads,minimametrics,200,1);
            labels_avg = mean(labels==0,2);
            %save(gifti(single(labels_avg)),[outputdir '/' gradsname '_smooth' num2str(smooth) '_wateredge_avg.func.gii']);
            %cifti_write_wHDR(single(labels_avg),'/data/cn4/evan/ROIs/subcort_corr_template.func.gii',[outputdir '/' gradsname '_smooth' num2str(smooth) '_wateredge_avg'])
            output = zeros(66697,1);
            output(end-(size(labels_avg)-1):end) = labels_avg;
            cifti_write_wHDR(output,[],[outputdir '/' gradsname '_smooth' num2str(smooth) '_wateredge_avg'])
            
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

function minimametric = metric_minima_all_subcort(metric,neighdist)

neighbors = cifti_subcort_neighbors('/data/hcp-zfs/home/laumannt/120_parcellation/cifti_timeseries_normalwall/vc25125_BOLD_LR_surf_subcort_32k_fsLR_smooth2.55.func.gii');
         

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

function labels = watershed_algorithm_all_par_subcort(edgemetrics,minimametrics,stepnum,fracmaxh)

neighbors = cifti_subcort_neighbors('/data/hcp-zfs/home/laumannt/120_parcellation/cifti_timeseries_normalwall/vc25125_BOLD_LR_surf_subcort_32k_fsLR_smooth2.55.func.gii');

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