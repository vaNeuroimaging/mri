function Cluster_correct_correlation_nifti(datafile,regressor,height_threshs,alphas,outputstem)
%Cluster_correct_correlation_cifti(datafile,regressor,height_threshs,alphas,outputstem)
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

warning off
iterations = 100;

if ischar(regressor)
    regressor = load(regressor);
end


template_nifti = load_untouch_nii_2D(datafile);
data1 = template_nifti.img;

if size(data1,2) ~= length(regressor)
    error('Different number of subjects in data ciftis and regressor!')
end

data1(isnan(data1)) = 0;

true_corrvals = paircorr_mod(data1',regressor(:));
true_corrvals(isnan(true_corrvals)) = 0;
template_nifti.img = true_corrvals;
template_nifti.hdr.dime.dim(5) = 1;
save_untouch_nii_2D(template_nifti,[outputstem '_uncorrected.nii.gz'])


max_random_clustersizes = zeros(iterations,1);


for height_thresh_num = 1:length(height_threshs)
    height_thresh = height_threshs(height_thresh_num);
    string = [];
    for iternum = 1:iterations
        
        prevstring = string;
        string = ['Conducting permuted correlations with an R threshold of ' num2str(height_thresh) ': iteration ' num2str(iternum) ' of ' num2str(iterations)];
        fprintf([repmat('\b',1,length(prevstring)) '%s'],string); 
        permuted_regressor = regressor(randperm(length(regressor)));
        %permuted_data = data1(:,randperm(size(data1,2)));
        
        
        corrvals = paircorr_mod(data1',permuted_regressor(:));
        corrvals(isnan(corrvals)) = 0;
        
        template_nifti.img = corrvals;
        save_untouch_nii_2D(template_nifti,[outputstem '_temp.nii.gz'])
        delete('temp.txt')
        [~,~] = system(['cluster -i ' outputstem '_temp.nii.gz -t ' num2str(height_thresh) ' >>temp.txt']);
        [~,permuted_clustersizes,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~] = textread('temp.txt','%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','headerlines',1,'emptyvalue',0);
        
        template_nifti.img = -corrvals;
        save_untouch_nii_2D(template_nifti,[outputstem '_temp.nii.gz'])
        delete('temp.txt')
        [~,~] = system(['cluster -i ' outputstem '_temp.nii.gz -t ' num2str(height_thresh) ' >>temp.txt']);
        [~,neg_permuted_clustersizes,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~,~] = textread('temp.txt','%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f','headerlines',1,'emptyvalue',0);
        delete('temp.txt')
        
%         permuted_clusteredRs = cifti_cluster_surface(abs(corrvals),height_thresh,inf,0,surfaceareas,neighbors);
%         permuted_clustersizes = zeros(size(permuted_clusteredRs,2),1);
%         for i = 1:size(permuted_clusteredRs,2)
%             permuted_clustersizes(i) = sum(permuted_clusteredRs(:,i) .* surfaceareas,1);
%         end
        max_random_clustersizes(iternum) = max([permuted_clustersizes ; neg_permuted_clustersizes]);
        if max_random_clustersizes(iternum) > 200
            1;
        end
    end
    disp(' ')
    
    sorted_max_random_clustersizes = sort(max_random_clustersizes,'descend');
    
    
        
    for alphanum = 1:length(alphas)
        
        alpha = alphas(alphanum);
        
        clustersize_thresh = sorted_max_random_clustersizes(round(iterations * alpha));
        disp(['Height threshold = ' num2str(height_thresh) ', Corrected alpha = ' num2str(alpha) ': cluster size threshold = ' num2str(clustersize_thresh) 'voxels'])
        
        [~,~] = system(['cluster -i ' outputstem '_uncorrected.nii.gz -t ' num2str(height_thresh) ' --minextent=' num2str(clustersize_thresh) ' --othresh=' outputstem '_corrected_pos.nii.gz']);
        [~,~] = system(['fslmaths ' outputstem '_uncorrected.nii.gz -mul -1 ' outputstem '_uncorrected_neg.nii.gz']);
        [~,~] = system(['cluster -i ' outputstem '_uncorrected_neg.nii.gz -t ' num2str(height_thresh) ' --minextent=' num2str(clustersize_thresh) ' --othresh=' outputstem '_corrected_neg.nii.gz']);
        [~,~] = system(['fslmaths ' outputstem '_corrected_neg.nii.gz -mul -1 -add ' outputstem '_corrected_pos.nii.gz ' outputstem '_corrected_r' num2str(height_thresh) '_alpha' num2str(alpha) '_k' num2str(clustersize_thresh) '.nii.gz']);
        
        
    end
end
delete([outputstem '_temp.nii.gz'])



  
    
