function Batch_homogeneity_testing_cifti(cifti_parcellations_totest,cov_corr_files,outdir,rotationsfile,iterations,nsurfverts,baddataname)
%Batch_homogeneity_testing_cifti(cifti_parcellations_totest,cov_corr_files,[rotationsfile],[iterations],[nsurfverts],baddatanames)

if ~exist('baddataname') || isempty(baddataname)
    baddataname = '/data/cn/data1/scripts/CIFTI_RELATED/Resources/Baddata_bigcluster_LR.dtseries.nii';
end

if ~exist('nsurfverts','var') || isempty(nsurfverts)
    nsurfverts = 32492;
end

if ~exist('iterations','var') || isempty(iterations)
    iterations = 1000;
end

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

baddata = ft_read_cifti_mod(baddataname);

for i = 1:length(cifti_parcellations_totest)
    this_cifti = cifti_parcellations_totest{i};
    [path, temp, ~] = fileparts(this_cifti);
    [~,filestem,~] = fileparts(temp);
    data = ft_read_cifti_mod(this_cifti);
    
    L_data = zeros(nsurfverts,1);
    L_datainds = data.brainstructure(1:nsurfverts)==1;
    L_data(L_datainds) = data.data(data.brainstructure==1);
    L_name = [outdir '/' filestem '_L'];
    save(gifti(single(L_data)),[L_name '.func.gii']);
    
    L_baddata = zeros(nsurfverts,1);
    L_baddatainds = baddata.brainstructure(1:nsurfverts)==1;
    L_baddata(L_baddatainds) = baddata.data(data.brainstructure==1);
    
    cov_corr_L = smartload(cov_corr_files{1});
    
    if size(cov_corr_L,1) ~= nnz(L_datainds)
        error('Covariance matrices must be from the same space as the parcellation!')
    end
    
    generate_rotated_parcellation_andPCA_nobaddata([L_name '.func.gii'],iterations,cov_corr_L,rotations,L_datainds,'L',L_name,L_baddata);
    
    clear cov_corr_L
    
    
    
    R_data = zeros(nsurfverts,1);
    R_datainds = [data.brainstructure((nsurfverts+1) : (nsurfverts * 2))==2];
    R_data(R_datainds) = data.data(data.brainstructure==2);
    R_name = [outdir '/' filestem '_R'];
    save(gifti(single(R_data)),[R_name '.func.gii']);
    
    R_baddata = zeros(nsurfverts,1);
    R_baddatainds = [baddata.brainstructure((nsurfverts+1) : (nsurfverts * 2))==2];
    R_baddata(R_baddatainds) = baddata.data(data.brainstructure==2);
    
    cov_corr_R = smartload(cov_corr_files{2});
    
    if size(cov_corr_R,1) ~= nnz(R_datainds)
        error('Covariance matrices must be from the same space as the parcellation!')
    end
    
    generate_rotated_parcellation_andPCA_nobaddata([R_name '.func.gii'],iterations,cov_corr_R,rotations,R_datainds,'R',R_name,R_baddata);
    
    clear cov_corr_R
    
end
    

