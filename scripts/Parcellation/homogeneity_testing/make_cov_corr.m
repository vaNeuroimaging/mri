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

if isstruct(dconn)
    data = dconn;
    clear dconn
else
    data = ft_read_cifti_mod(dconn);
end
Linds = 1:nnz(data.brainstructure==1);
Rinds = (nnz(data.brainstructure==1)+1):(nnz(data.brainstructure==1) + nnz(data.brainstructure==2));



cov_corr_L = cov(data.data(:,Linds));

out = data;
out.data = cov_corr_L;
out.brainstructure = data.brainstructure(1:32492);
out.pos = data.pos(1:32492,1);
out.brainstructurelabel = data.brainstructurelabel(1);
try out = rmfield(out,{'dim','transform'}); out = rmfield(out,'time');
catch; end
ft_write_cifti_mod([outstem '_L.dconn.nii'],out);

%save([outstem '_L.mat'],'cov_corr_L','-v7.3')
clear cov_corr_L out



cov_corr_R = cov(data.data(:,Rinds));

out = data;
out.data = cov_corr_R;
out.brainstructure = data.brainstructure(32493:64984);
out.brainstructure(out.brainstructure>0) = 1;
out.pos = data.pos(32493:64984,1);
out.brainstructurelabel = data.brainstructurelabel(2);
try out = rmfield(out,{'dim','transform'}); out = rmfield(out,'time');
catch; end
ft_write_cifti_mod([outstem '_R.dconn.nii'],out);

%save([outstem '_R.mat'],'cov_corr_R','-v7.3')
clear cov_corr_R out








