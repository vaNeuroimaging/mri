function Cluster_correct_ttest2_cifti(data1file,data2file,height_threshs,alphas,outputstem,surfacearea_cifti)
% output_matrix = Cluster_correct_ttest2_cifti(data1file,data2file,height_threshs,alphas,outputstem,[surfacearea_cifti])
%
% Conducts a two-sample t-test between the data in the matrix data1 and the
% matrix data2. Corrects for multiple comparisons at the cluster level by
% randomly permuting the group identities 1,000 times, conducting
% two-sample t-tests, and assessing how often clusters of a given size /
% t-statistic threshold emerge in the random data.
%
% Inputs:
% data1file: path to a cifti containing a [#vertices X #group1 subjects] matrix.
% data2file: path to a cifti containing a [#vertices X #group2 subjects] matrix.
% height_threshs : a vector of a priori t-statistic thresholds to use.
% alphas: a vector containing desired alphas to correct to.
% outputstem: the base name for the uncorrected and cluster-corrected output t-images.
% surfacearea_cifti: an optional argument containing the surface areas of
%   the surfaces to be used. This should be a surface-only cifti file in
%   the same space as data1 and data2. Leave empty or omit to use the
%   default, which is the surface area of the fs_LR 32k Conte69 atlas with
%   the standard medial wall.
%
% E.Gordon 12/2014


iterations = 1000;
default_surfaceareas = '/home/data/scripts/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.LR.midthickness.32k_fs_LR_surfaceareas_cortexonly.dtseries.nii';
%default_neighbors = '/data/cn4/evan/fsaverage_LR32k/Cifti_surf_neighbors_LR_normalwall.mat';

if exist('surfacearea_cifti') && ~isempty(surfacearea_cifti)
    surfaceareas = ft_read_cifti_mod(surfacearea_cifti);
else
    surfaceareas = ft_read_cifti_mod(default_surfaceareas);
end
surfaceareas = surfaceareas.data;


template_cifti = ft_read_cifti_mod(data1file);
neighbors = cifti_neighbors(data1file);

data1 = ft_read_cifti_mod(data1file); data1 = data1.data;
data2 = ft_read_cifti_mod(data2file); data2 = data2.data;

if size(data1,1) ~= size(data2,1)
    error('Different number of data points in data1 cifti and data2 cifti!')
end

if size(data1,1) ~= size(surfaceareas,1)
    error('Different number of data points in data ciftis and surface area cifti!')
end


[~,~,~,TRUESTATS] = ttest2(data1',data2');
    template_cifti.data = TRUESTATS.tstat';
    template_cifti.dimord = 'pos_time';
    ft_write_cifti_mod([outputstem '_uncorrected'],template_cifti);


template_cifti.dimord = 'pos_scalar';
template_cifti.mapname = cell(1,0);

concatenated = [data1 data2];

nanspresent = 0;
if nnz(isnan(concatenated)) > 0
    disp('NaNs present in data--iterating t-tests will run slower')
    nanspresent = 1;
end

groups = [ones(size(data1,2),1) ; ones(size(data2,2),1)*2];

max_random_clustersizes = zeros(iterations,1);

output_matrix = zeros(size(data1,1),0);

for height_thresh_num = 1:length(height_threshs)
    height_thresh = height_threshs(height_thresh_num);
    string = [];
    for iternum = 1:iterations
        
        prevstring = string;
        string = ['Conducting permuted t-tests with a t-threshold of ' num2str(height_thresh) ': iteration ' num2str(iternum) ' of ' num2str(iterations)];
        fprintf([repmat('\b',1,length(prevstring)) '%s'],string); 
        
        permuted_groups = groups(randperm(length(groups)));
        if nanspresent
            [~,~,~,STATS] = ttest2(concatenated(:,permuted_groups==1)',concatenated(:,permuted_groups==2)');
        else
            [~,~,~,STATS] = ttest2_faster(concatenated(:,permuted_groups==1)',concatenated(:,permuted_groups==2)');
        end
        permuted_clusteredTs = cifti_cluster_surface(abs(STATS.tstat'),height_thresh,inf,0,surfaceareas,neighbors);
        permuted_clustersizes = zeros(size(permuted_clusteredTs,2),1);
        for i = 1:size(permuted_clusteredTs,2)
            permuted_clustersizes(i) = sum(permuted_clusteredTs(:,i) .* surfaceareas,1);
        end
        max_random_clustersizes(iternum) = max(permuted_clustersizes);
    end
    disp(' ')
    
    sorted_max_random_clustersizes = sort(max_random_clustersizes,'descend');
    
    
        
    for alphanum = 1:length(alphas)
        
        alpha = alphas(alphanum);
        
        clustersize_thresh = sorted_max_random_clustersizes(round(iterations * alpha));
        disp(['Height threshold = ' num2str(height_thresh) ', Corrected alpha = ' num2str(alpha) ': cluster size threshold = ' num2str(clustersize_thresh) 'mm'])
        
        clusteredTs = cifti_cluster_surface(abs(TRUESTATS.tstat'),height_thresh,inf,clustersize_thresh,surfaceareas,neighbors);
        output = TRUESTATS.tstat' .* sum(clusteredTs,2);
        
        %if ~isempty(outputstem)
        %    cifti_write_wHDR(output,template_cifti,[outputstem '_corrected_T' num2str(height_thresh) '_k' num2str(clustersize_thresh) '_alpha' num2str(alpha)])
        %end
        
        output_matrix(:,end+1) = output;
        template_cifti.mapname{1,end+1} = ['Height = ' num2str(height_thresh) ', Alpha = ' num2str(alpha) ': k = ' num2str(clustersize_thresh) 'mm'];
        
    end
end

template_cifti.data = output_matrix;
ft_write_cifti_mod([outputstem '_corrected'],template_cifti);


  
    
