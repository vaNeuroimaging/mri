

for i = 4:42
    
    allparcelfilenames(i-3,:) = {['/data/cn4/evan/ROIs/Craddock/Craddock' sprintf('%04i',i) '_L.func.gii'],['/data/cn4/evan/ROIs/Craddock/Craddock' sprintf('%04i',i) '_R.func.gii']};
end
allparcelfilenames(end+1,:) = {'/data/cn4/evan/ROIs/Shen/Shen_200_L.func.gii','/data/cn4/evan/ROIs/Shen/Shen_200_R.func.gii'};
allparcelfilenames(end+1,:) = {'/data/cn4/evan/ROIs/Shen/Shen_100_L.func.gii','/data/cn4/evan/ROIs/Shen/Shen_100_R.func.gii'};



iterations = 100;



load('rotations.mat')



load('/data/cn5/selfRegulation/V4Process_nosmooth/cifti_covariance_normalwall/cov_corr_L_normalwall.mat')
for i = 1:size(allparcelfilenames,1)
    slashloc = strfind(allparcelfilenames{i,1},'/');
    outfileprefix = [allparcelfilenames{i,1}(slashloc(end)+1 : end-9)];
    
    generate_rotated_parcellation_andPCA_nobaddata(allparcelfilenames{i,1},iterations,cov_corr_L,rotations,1,'L',outfileprefix);
    
end
clear cov_corr_L


load('/data/cn5/selfRegulation/V4Process_nosmooth/cifti_covariance_normalwall/cov_corr_R_normalwall.mat')
for i = 1:size(allparcelfilenames,1)
    slashloc = strfind(allparcelfilenames{i,2},'/');
    outfileprefix = [allparcelfilenames{i,2}(slashloc(end)+1 : end-9)];
    generate_rotated_parcellation_andPCA_nobaddata(allparcelfilenames{i,2},iterations,cov_corr_R,rotations,1,'R',outfileprefix);
end
clear cov_corr_R
