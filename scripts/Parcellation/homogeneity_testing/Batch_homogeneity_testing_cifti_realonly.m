function Batch_homogeneity_testing_cifti(cifti_parcellations_totest,cov_corr_files,outdir,rotationsfile,iterations,baddataname)
%Batch_homogeneity_testing_cifti(cifti_parcellations_totest,cov_corr_files,[outdir],[rotationsfile],[iterations],[baddataname])
%
% Tests the homogeneity of connectivity values in each parcel of multiple
% parcellations, generates rotated versions of those same parcellations,
% and tests homogeneity in the rotated versions
%
% Inputs:
% cifti_parcellations_totest - a string containing the location of a
%   parcellation, or a cell array of strings containing multiple
%   parcellations.
% cov_corr_files - a cell array with two strings: one with the location of
%   the .mat file containing the left hemisphere covariance matrix, and the
%   other with the location of the right hemisphere covaraince matrix
% outdir - location results will be written into. Omit or leave empty ([])
%   to use the current working directory.
% rotationsfile - an optional string input describing the location of a
%   .mat file containing random rotations around the x,y,and z axes. Useful
%   to specify e.g. to use the same random rotations as previous tests. Omit
%   or leave empty to generate a new set of rotations.
% iterations - an optional argument listing the number of random rotations
%   to generate and test. Omit or leave empty to use the default of 1000.
% baddataname - an optional argument identifying a binary cifti file
%   describing regions of the cortex in which parcel homogenities should not
%   be tested (e.g. because they have low SNR). Omit or leave empty to use a
%   default file: /data/cn/data1/scripts/CIFTI_RELATED/Resources/Baddata_bigcluster_LR.dtseries.nii.
%   Alternately, specify 'none' to test homogeneity everywhere.
%
% Requires the Cifti Resources scripts to be in your path (e.g.,
% /data/cn/data1/scripts/CIFTI_RELATED/Resources/ and subfolders)
%
%EMG 08/26/15


%Set defaults for missing inputs
if ~exist('outdir','var') || isempty(outdir)
    outdir = pwd;
end

if ~exist('baddataname','var') || isempty(baddataname)
    baddataname = '/home/data/scripts/Resources/Baddata_bigcluster_LR.dtseries.nii';
end

if ~exist('iterations','var') || isempty(iterations)
    iterations = 1000;
end

%Make a new set of rotations
if ~exist('rotationsfile','var') || isempty(rotationsfile)
    for iternum = 1:iterations
        
        rotations.xrot(iternum) = rand * 2*pi;
        rotations.yrot(iternum) = rand * 2*pi;
        rotations.zrot(iternum) = rand * 2*pi;
        
    end
    save('rotations.mat','rotations')
else
    rotations = smartload(rotationsfile);
end


nsurfverts = 32492;


L_baddata = zeros(nsurfverts,1);
R_baddata = zeros(nsurfverts,1);

%Load the baddata and convert to hemisphere surfaces
if ~strcmp(baddataname,'none')
    
    baddata = ft_read_cifti_mod(baddataname);
    
    L_baddatainds = baddata.brainstructure(1:nsurfverts)==1;
    L_baddata(L_baddatainds) = baddata.data(1:nnz(baddata.brainstructure==1));
    
    R_baddatainds = [baddata.brainstructure((nsurfverts+1) : (nsurfverts * 2))==2];
    R_baddata(R_baddatainds) = baddata.data((nnz(baddata.brainstructure==1)+1) : (nnz(baddata.brainstructure==1) + nnz(baddata.brainstructure==2)));
end


if ischar(cifti_parcellations_totest)
    cifti_parcellations_totest = {cifti_parcellations_totest};
end

%Loop through parcellations
for i = 1:length(cifti_parcellations_totest)
    this_cifti = cifti_parcellations_totest{i};
    [~, temp, ~] = fileparts(this_cifti);
    [~,filestem,~] = fileparts(temp);
    
    %Load the cifti
    data = ft_read_cifti_mod(this_cifti);
    
    %Save the cifti on the left hemisphere
    L_data = zeros(nsurfverts,1);
    L_datainds = data.brainstructure(1:nsurfverts)==1;
    L_data(L_datainds) = data.data(1:nnz(data.brainstructure==1));
    L_name = [outdir '/' filestem '_L'];
    save(gifti(single(L_data)),[L_name '.func.gii']);
    
    %Load the covariance
    cov_corr_L = ft_read_cifti_mod(cov_corr_files{1});
    cov_corr_L = cov_corr_L.data;
    
    if size(cov_corr_L,1) ~= nnz(L_datainds)
        error('Covariance matrices must be from the same space as the parcellation!')
    end
    
    %Run the homogeneity testing
    parcel_homogeneity_nobaddata_generic([L_name '.func.gii'],iterations,cov_corr_L,rotations,L_datainds,'L',L_name,L_baddata);
    
    clear cov_corr_L
    
    
    
    
    %Save the cifti on the right hemisphere
    R_data = zeros(nsurfverts,1);
    R_datainds = [data.brainstructure((nsurfverts+1) : (nsurfverts * 2))==2];
    R_data(R_datainds) = data.data((nnz(data.brainstructure==1)+1) : (nnz(data.brainstructure==1) + nnz(data.brainstructure==2)));
    R_name = [outdir '/' filestem '_R'];
    save(gifti(single(R_data)),[R_name '.func.gii']);
    
    %Load the covariance
    cov_corr_R = ft_read_cifti_mod(cov_corr_files{2});
    cov_corr_R = cov_corr_R.data;
    
    if size(cov_corr_R,1) ~= nnz(R_datainds)
        error('Covariance matrices must be from the same space as the parcellation!')
    end
    
    %Run the homogeneity testing
    parcel_homogeneity_nobaddata_generic([R_name '.func.gii'],iterations,cov_corr_R,rotations,R_datainds,'R',R_name,R_baddata);
    
    clear cov_corr_R
    
    
    
    %Combine results from two hemispheres
    
%     L_results = load([L_name '.mat']);
%     R_results = load([R_name '.mat']);
%     
%     realsizes = [L_results.realsizes  R_results.realsizes];
%     realeigvals_per_first = [L_results.realeigvals_per_first  R_results.realeigvals_per_first];
%     rotated_eigvals = [L_results.rotated_eigvals ; R_results.rotated_eigvals];
%     
%     save([filestem '.mat'],'realsizes','realeigvals_per_first','rotated_eigvals');
    
    delete([L_name '.mat']);
    delete([R_name '.mat']);
    
    
    L_homogeneity = gifti([outdir '/PCA_eigval_per_first_' filestem '_L.func.gii']);
    R_homogeneity = gifti([outdir '/PCA_eigval_per_first_' filestem '_R.func.gii']);
    data.data(1:nnz(L_datainds),1) = L_homogeneity.cdata(L_datainds);
    data.data((nnz(L_datainds)+1) : (nnz(L_datainds) + nnz(R_datainds))) = R_homogeneity.cdata(R_datainds);
    ft_write_cifti_mod([outdir '/PCA_eigval_per_first_' filestem '.dtseries.nii'],data)
    
    delete([outdir '/PCA_eigval_per_first_' filestem '_L.func.gii']);
    delete([outdir '/PCA_eigval_per_first_' filestem '_R.func.gii']);
    
    
end
    

