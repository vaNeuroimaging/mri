function Cluster_correct_correlation_cifti(datafile,regressor,height_threshs,alphas,outputstem,surfacearea_cifti)
%Cluster_correct_correlation_cifti(datafile,regressor,height_threshs,alphas,outputstem,[surfacearea_cifti])
%
% Conducts a correlation between the data in the matrix datafile and the
% regressor. Corrects for multiple comparisons at the cluster level by
% randomly permuting the group identities 1,000 times, conducting
% correlations, and assessing how often clusters of a given size /
% R threshold emerge in the random data.
%
% Inputs:
% datafile: path to a cifti containing a [#vertices X # subjects] matrix.
% regressor: a vector of length [#subjects] containing the values to be
%   regressed against the data. Can also be the path to a text file
%   containing this vector.
% height_threshs : a vector of a priori t-statistic thresholds to use.
% alphas: a vector containing desired alphas to correct to.
% outputstem: the base name for the uncorrected and cluster-corrected output t-images.
% surfacearea_cifti: an optional argument containing the surface areas of
%   the surfaces to be used. This should be a surface-only cifti file in
%   the same space as data1 and data2. Leave empty or omit to use the
%   default, which is the surface area of the fs_LR 32k Conte69 atlas with
%   the standard medial wall.
%
% E.Gordon 3/2016


iterations = 100;
default_surfaceareas = '/home/data/scripts/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.LR.midthickness.32k_fs_LR_surfaceareas_cortexonly.dtseries.nii';
%default_neighbors = '/data/cn4/evan/fsaverage_LR32k/Cifti_surf_neighbors_LR_normalwall.mat';

if exist('surfacearea_cifti') && ~isempty(surfacearea_cifti)
    surfaceareas = ft_read_cifti_mod(surfacearea_cifti);
else
    surfaceareas = ft_read_cifti_mod(default_surfaceareas);
end
surfaceareas = surfaceareas.data;

if ischar(regressor)
    regressor = load(regressor);
end


template_cifti = ft_read_cifti_mod(datafile);
data1 = template_cifti.data;
neighbors = cifti_neighbors(datafile);

if size(data1,1) ~= size(surfaceareas,1)
    error('Different number of data points in data ciftis and surface area cifti!')
end

if size(data1,2) ~= length(regressor)
    error('Different number of subjects in data ciftis and regressor!')
end

data1(isnan(data1)) = 0;

true_corrvals = paircorr_mod(data1',regressor(:));
template_cifti.data = true_corrvals;
template_cifti.dimord = 'pos_time';
ft_write_cifti_mod([outputstem '_uncorrected'],template_cifti);


template_cifti.dimord = 'pos_scalar';
template_cifti.mapname = cell(1,0);

max_random_clustersizes = zeros(iterations,1);

output_matrix = zeros(size(data1,1),0);

for height_thresh_num = 1:length(height_threshs)
    height_thresh = height_threshs(height_thresh_num);
    string = [];
    for iternum = 1:iterations
        
        prevstring = string;
        string = ['Conducting permuted correlations with an R threshold of ' num2str(height_thresh) ': iteration ' num2str(iternum) ' of ' num2str(iterations)];
        fprintf([repmat('\b',1,length(prevstring)) '%s'],string); 
        permuted_data = data1(:,randperm(size(data1,2)));
        
        
        corrvals = paircorr_mod(permuted_data',regressor(:));
        
        permuted_clusteredRs = cifti_cluster_surface(abs(corrvals),height_thresh,inf,0,surfaceareas,neighbors);
        permuted_clustersizes = zeros(size(permuted_clusteredRs,2),1);
        for i = 1:size(permuted_clusteredRs,2)
            permuted_clustersizes(i) = sum(permuted_clusteredRs(:,i) .* surfaceareas,1);
        end
        max_random_clustersizes(iternum) = max(permuted_clustersizes);
    end
    disp(' ')
    
    sorted_max_random_clustersizes = sort(max_random_clustersizes,'descend');
    
    
        
    for alphanum = 1:length(alphas)
        
        alpha = alphas(alphanum);
        
        clustersize_thresh = sorted_max_random_clustersizes(round(iterations * alpha));
        disp(['Height threshold = ' num2str(height_thresh) ', Corrected alpha = ' num2str(alpha) ': cluster size threshold = ' num2str(clustersize_thresh) 'mm'])
        
        clusteredRs = cifti_cluster_surface(abs(true_corrvals),height_thresh,inf,clustersize_thresh,surfaceareas,neighbors);
        output = true_corrvals .* sum(clusteredRs,2);
        
        output_matrix(:,end+1) = output;
        template_cifti.mapname{1,end+1} = ['Height = ' num2str(height_thresh) ', Alpha = ' num2str(alpha) ': k = ' num2str(clustersize_thresh) 'mm'];
        
    end
end

template_cifti.data = output_matrix;
ft_write_cifti_mod([outputstem '_corrected'],template_cifti);


  
    
