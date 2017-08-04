allparcelfilenames ={'/data/cn4/evan/Temp/bothhem_parcellation_test/avg_corrmat/Parcels_avg_corrmat_edgethresh_0.45_L.func.gii','/data/cn4/evan/Temp/bothhem_parcellation_test/avg_corrmat/Parcels_avg_corrmat_edgethresh_0.45_R.func.gii'};
%};
%allparcelfilenames = {'/data/cn4/evan/RestingState/Ind_variability/Poldrome/Poldrome_84_subsurf_edge_L_watershedmerge_0.45_tweaked.func.gii','/data/cn4/evan/RestingState/Ind_variability/Poldrome/Poldrome_84_subsurf_edge_R_watershedmerge_0.45_tweaked.func.gii'};
%{'/data/cn4/evan/ROIs/AAL222_MNI_L.func.gii','/data/cn4/evan/ROIs/AAL222_MNI_R.func.gii'};
%allparcelfilenames = {'/data/cn4/evan/RestingState/FC_Mapping_120/subsurf/120_subsurf_L_watershedmerge_0.5_tweaked.func.gii','/data/cn4/evan/RestingState/FC_Mapping_120/subsurf/120_subsurf_R_watershedmerge_0.5_tweaked.func.gii'};
 %allparcelfilenames = {'/data/cn4/evan/ROIs/264_surfvert_ROIs_L_removeedge.func.gii','/data/cn4/evan/ROIs/264_surfvert_ROIs_R_removeedge.func.gii';...
%     '/data/cn4/evan/RestingState/FC_Mapping_120/ConsensusParcels_L_removeedge.func.gii','/data/cn4/evan/RestingState/FC_Mapping_120/ConsensusParcels_R_removeedge.func.gii';...
% '/data/cn4/evan/ROIs/Craddock/Craddock_350_L_removeedge.func.gii','/data/cn4/evan/ROIs/Craddock/Craddock_350_R_removeedge.func.gii';...
% '/data/cn4/evan/RestingState/FC_Mapping_120/Yeo_17_L_parcels_removeedge.func.gii','/data/cn4/evan/RestingState/FC_Mapping_120/Yeo_17_R_parcels_removeedge.func.gii';...
% '/data/cn4/evan/ROIs/mode.L.aparc.a2009s.32k_fs_LR_removeedge.func.gii','/data/cn4/evan/ROIs/mode.R.aparc.a2009s.32k_fs_LR_removeedge.func.gii'};

division = '';

iterations = 1000;

makenewrotations = 0;


if makenewrotations

for iternum = 1:iterations
    
       rotations.xrot(iternum) = rand * 2*pi;
       rotations.yrot(iternum) = rand * 2*pi;
       rotations.zrot(iternum) = rand * 2*pi;
       
end

save('rotations.mat','rotations')

else
    
    load(['rotations' division '.mat'])
    
end

    
    
    for i = 1:size(allparcelfilenames,1)
        load('/data/cn5/selfRegulation/V4Process_nosmooth/cifti_covariance_normalwall/cov_corr_L_normalwall.mat')
        %load /data/cn4/evan/RestingState/Ind_variability/Poldrome/cov_corr_L.mat
            slashloc = strfind(allparcelfilenames{i,1},'/');
            outfileprefix = [allparcelfilenames{i,1}(slashloc(end)+1 : end-9) division];
            
            %if ~exist([outfileprefix '.mat'])
            
                generate_rotated_parcellation_andPCA_nobaddata(allparcelfilenames{i,1},iterations,cov_corr_L,rotations,1,'L',outfileprefix);
                %generate_rotated_parcellation_andPCA_nobaddata(allparcelfilenames{i,1},iterations,cov_corr_L,rotations,1,'L',outfileprefix);
            %end
            
        clear cov_corr_L
    
    
    load('/data/cn5/selfRegulation/V4Process_nosmooth/cifti_covariance_normalwall/cov_corr_R_normalwall.mat')
    %load /data/cn4/evan/RestingState/Ind_variability/Poldrome/cov_corr_R.mat
    %for i = 1:size(allparcelfilenames,1)
        slashloc = strfind(allparcelfilenames{i,2},'/');
    outfileprefix = [allparcelfilenames{i,2}(slashloc(end)+1 : end-9) division];
    generate_rotated_parcellation_andPCA_nobaddata(allparcelfilenames{i,2},iterations,cov_corr_R,rotations,1,'R',outfileprefix);
    %generate_rotated_parcellation_andPCA_nobaddata(allparcelfilenames{i,2},iterations,cov_corr_R,rotations,1,'R',outfileprefix);
    %end
    clear cov_corr_R
    end
