function Cluster_correct_ttestpaired_cifti(data1,data2,height_threshs,alphas,outputstem,surfacearea_cifti,template_cifti)
% Cluster_correct_ttest2_cifti(data1,data2,height_threshs,alphas,outputstem,[surfacearea_cifti],[template_cifti])
%
% Conducts a paired t-test between the data in the matrix data1 and the
% matrix data2. Corrects for multiple comparisons at the cluster level by
% randomly permuting the condition identities 10,000 times, conducting
% paired t-tests, and assessing how often clusters of a given size /
% t-statistic threshold emerge in the random data.
%
% Inputs:
% data1: a [#vertices X #subjects] matrix, or a dtseries containing such a matrix.
% data2: a [#vertices X #subjects] matrix, or a dtseries containing such a matrix.
% height_threshs : a vector of a priori t-statistic thresholds to use.
% alphas: a vector containing desired alphas to correct to.
% outputstem: the base name for the cluster-corrected output t-image. The
%   height threshold, cluster extent threshold, and alpha will be appended.
% surfacearea_cifti: an optional argument containing the surface areas of
%   the surfaces to be used. This should be a surface-only cifti file in
%   the same space as data1 and data2. Leave empty or omit to use the
%   default, which is the surface area of the fs_LR 32k Conte69 atlas with
%   the standard medial wall.
% template_cifti: an optional argument containing any cifti file in the
%   same space as data1 and data2, used to get the vertex neighbor
%   adjacencies and to write the output. Leave empty or omit to use the
%   default, which is the fs_LR 32k Conte69 atlas witht he standard medial
%   wall.
%
% E.Gordon 12/2014


iterations = 10000;
default_surfaceareas = '/data/cn4/evan/fsaverage_LR32k/Conte69.LR.midthickness.32k_fs_LR_surfaceareas_normalwall.dtseries.nii';
default_neighbors = '/data/cn4/evan/fsaverage_LR32k/Cifti_surf_neighbors_LR_normalwall.mat';

if exist('surfacearea_cifti') && ~isempty(surfacearea_cifti)
    surfaceareas = cifti_read(surfacearea_cifti);
else
    surfaceareas = cifti_read(default_surfaceareas);
end

if exist('template_cifti') && ~isempty(template_cifti)
    neighbors = cifti_neighbors(template_cifti);
else
    neighbors = smartload(default_neighbors);
    template_cifti = default_surfaceareas;
end

if ischar(data1)
    data1 = cifti_read(data1);
end
if ischar(data2)
    data2 = cifti_read(data2);
end

concatenated = [data1 data2];

nanspresent = 0;
if nnz(isnan(concatenated)) > 0
    disp('NaNs present in data--iterating t-tests will run slower')
    nanspresent = 1;
end

groups = [ones(size(data1,2),1) ; ones(size(data2,2),1)*2];

max_random_clustersizes = zeros(iterations,1);



for height_thresh = height_threshs
    string = [];
    for iternum = 1:iterations
        
        prevstring = string;
        string = ['Conducting permuted t-tests with a t-threshold of ' num2str(height_thresh) ': iteration ' num2str(iternum) ' of ' num2str(iterations)];
        fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
        
        permuted_groups = groups(randperm(length(groups)));
        if nanspresent
            [ign,ign2,ign3,STATS] = ttest(concatenated(:,permuted_groups==1)',concatenated(:,permuted_groups==2)');
        else
            [ign,ign2,ign3,STATS] = ttest_faster(concatenated(:,permuted_groups==1)',concatenated(:,permuted_groups==2)');
        end
        permuted_clusteredTs = cifti_cluster_surface(abs(STATS.tstat'),height_thresh,inf,0,surfaceareas,neighbors);
        permuted_clustersizes = zeros(size(permuted_clusteredTs,2),1);
        for i = 1:size(permuted_clusteredTs,2)
            permuted_clustersizes(i) = sum(permuted_clusteredTs(:,i) .* surfaceareas,1);
        end
        max_random_clustersizes(iternum) = max(permuted_clustersizes);
    end
    disp('')
    
    sorted_max_random_clustersizes = sort(max_random_clustersizes,'descend');
    
    [ign,ign2,ign3,STATS] = ttest(data1',data2');
        
    for alpha = alphas
        
        clustersize_thresh = sorted_max_random_clustersizes(round(iterations * alpha));
        
        clusteredTs = cifti_cluster_surface(abs(STATS.tstat'),height_thresh,inf,clustersize_thresh,surfaceareas,neighbors);
        output = STATS.tstat' .* sum(clusteredTs,2);
        
        cifti_write_wHDR(output,template_cifti,[outputstem '_T' num2str(height_thresh) '_k' num2str(clustersize_thresh) '_alpha' num2str(alpha)])
        
    end
end


  
    
