function make_cov_corr(dconn,outstem)
%make_cov_corr(dconn,outstem)
%
% Make single-hemisphere covariance matrices for use in homogeneity
% testing. These matrices are each [hemisphere verts x hemisphere verts]
% and represent the covariance of the correlation patterns.
%
% Inputs:
%  dconn - a .dconn.nii file containing the (Fisher-transformed)
%   correlations of all points against all points
%  outstem - name of output files, which will be called [outstem]_L.mat and
%   [outstem]_R.mat
%
%EMG 08/26/15

data = ft_read_cifti_mod(dconn);
Linds = 1:nnz(data.brainstructure==1);
Rinds = (nnz(data.brainstructure==1)+1):(nnz(data.brainstructure==1) + nnz(data.brainstructure==2));

cov_corr_L = cov(data.data(:,Linds));
save([outstem '_L.mat'],'cov_corr_L','-v7.3')
clear cov_corr_L

cov_corr_R = cov(data.data(:,Rinds));
save([outstem '_L.mat'],'cov_corr_L','-v7.3')

