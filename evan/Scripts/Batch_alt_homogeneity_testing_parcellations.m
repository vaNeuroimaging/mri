allparcelfilenames = {'/data/cn4/evan/RestingState/FC_Mapping_120/subsurf/nosmooth/120_subsurf_L_nosmooth_watershedmerge_0.4_tweaked.func.gii','/data/cn4/evan/RestingState/FC_Mapping_120/subsurf/nosmooth/120_subsurf_R_nosmooth_watershedmerge_0.4_tweaked.func.gii';...
'/data/cn4/evan/ROIs/Shen/Shen_300_L.func.gii','/data/cn4/evan/ROIs/Shen/Shen_300_R.func.gii';...
    '/data/cn4/evan/ROIs/Craddock/Craddock0032_L.func.gii','/data/cn4/evan/ROIs/Craddock/Craddock0032_R.func.gii';...
    '/data/cn4/evan/ROIs/264_surfvert_ROIs_L.func.gii','/data/cn4/evan/ROIs/264_surfvert_ROIs_R.func.gii';...
    '/data/cn4/evan/ROIs/Power_Neuron_2011_L_parcels.func.gii','/data/cn4/evan/ROIs/Power_Neuron_2011_R_parcels.func.gii';...'/data/cn4/evan/RestingState/FC_Mapping_120/ConsensusParcels_L.func.gii','/data/cn4/evan/RestingState/FC_Mapping_120/ConsensusParcels_R.func.gii';...
    '/data/cn4/evan/RestingState/FC_Mapping_120/Yeo_17_L_parcels.func.gii','/data/cn4/evan/RestingState/FC_Mapping_120/Yeo_17_R_parcels.func.gii';...
    '/data/cn4/evan/ROIs/Brodmann_L.func.gii','/data/cn4/evan/ROIs/Brodmann_R.func.gii';...
    '/data/cn4/evan/ROIs/AAL222_MNI_L.func.gii','/data/cn4/evan/ROIs/AAL222_MNI_R.func.gii'};

%allparcelfilenames = {'/data/cn5/selfRegulation/V4Process_nosmooth/gradients_120_108_combined_subsurf_nosmooth/120_108_combined_watershedmerge_0.35_tweaked.func.gii','/data/cn5/selfRegulation/V4Process_nosmooth/gradients_120_108_combined_subsurf_nosmooth/120_108_combined_R_watershedmerge_0.35_tweaked.func.gii'};


cd /data/cn4/evan/RestingState/FC_Mapping_120/AlternateHomogeneity_testing_in_108/

%allparcelfilenames = {'/data/cn4/evan/RestingState/FC_Mapping_120/120_L_wateredgethresh_watershedmerge_0.45.func.gii','/data/cn4/evan/RestingState/FC_Mapping_120/120_R_wateredgethresh_watershedmerge_0.45.func.gii'};

iterations = 100;

makenewrotations = 0;


if makenewrotations
    
    for iternum = 1:iterations
        
        rotations.xrot(iternum) = rand * 2*pi;
        rotations.yrot(iternum) = rand * 2*pi;
        rotations.zrot(iternum) = rand * 2*pi;
        
    end
    
    save('rotations.mat','rotations')
    
else
    
    load('rotations_all.mat')
    %load('rotations.mat')
    
    rotations.xrot = rotations.xrot(1:iterations);
    rotations.yrot = rotations.yrot(1:iterations);
    rotations.zrot = rotations.zrot(1:iterations);
    
end


    load('/data/cn5/selfRegulation/V4Process_nosmooth/cifti_correlation_normalwall/108_avg_corr_L_normalwall.mat')
    crosscorr = FisherTransform(paircorr_mod(corr_L));
    clear corr_L

for i = 1:size(allparcelfilenames,1)
    
    
    
%     load('/data/cn4/evan/RestingState/FC_Mapping_120/cov_corr_normalwall_L.mat')
%     
    slashloc = strfind(allparcelfilenames{i,1},'/');
    outfileprefix = [allparcelfilenames{i,1}(slashloc(end)+1 : end-9)];
    
    generate_rotated_parcellation_andmeancrosscorr_nobaddata(allparcelfilenames{i,1},iterations,crosscorr,rotations,1,'L',outfileprefix);
        
    
end
    
load('/data/cn5/selfRegulation/V4Process_nosmooth/cifti_correlation_normalwall/108_avg_corr_R_normalwall.mat')
    crosscorr = FisherTransform(paircorr_mod(corr_R));
    clear corr_R

for i = 1:size(allparcelfilenames,1)
    
    
%    load('/data/cn4/evan/RestingState/FC_Mapping_120/cov_corr_normalwall_R.mat')
    slashloc = strfind(allparcelfilenames{i,2},'/');
    outfileprefix = [allparcelfilenames{i,2}(slashloc(end)+1 : end-9)];
    generate_rotated_parcellation_andmeancrosscorr_nobaddata(allparcelfilenames{i,2},iterations,crosscorr,rotations,1,'R',outfileprefix);
    
    
    
end



