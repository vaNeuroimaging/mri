% allparcelfilenames = {'/data/cn4/evan/RestingState/FC_Mapping_120/Cohorts/C1_Lwatershedmerge_0.45.func.gii','/data/cn4/evan/RestingState/FC_Mapping_120/Cohorts/C1_Rwatershedmerge_0.45.func.gii';...
% '/data/cn4/evan/ROIs/264_surfvert_ROIs_L.func.gii','/data/cn4/evan/ROIs/264_surfvert_ROIs_R.func.gii';...
% '/data/cn4/evan/RestingState/FC_Mapping_120/ConsensusParcels_L.func.gii','/data/cn4/evan/RestingState/FC_Mapping_120/ConsensusParcels_R.func.gii';...
% '/data/cn4/evan/ROIs/Craddock/Craddock_350_L.func.gii','/data/cn4/evan/ROIs/Craddock/Craddock_350_R.func.gii';...
% '/data/cn4/evan/RestingState/FC_Mapping_120/Yeo_17_L_parcels.func.gii','/data/cn4/evan/RestingState/FC_Mapping_120/Yeo_17_R_parcels.func.gii';...
% '/data/cn4/evan/ROIs/mode.L.aparc.a2009s.32k_fs_LR.func.gii','/data/cn4/evan/ROIs/mode.R.aparc.a2009s.32k_fs_LR.func.gii'};

% allparcelfilenames = {'/data/cn4/evan/ROIs/264_surfvert_ROIs_L_removeedge.func.gii','/data/cn4/evan/ROIs/264_surfvert_ROIs_R_removeedge.func.gii';...
%     '/data/cn4/evan/RestingState/FC_Mapping_120/ConsensusParcels_L_removeedge.func.gii','/data/cn4/evan/RestingState/FC_Mapping_120/ConsensusParcels_R_removeedge.func.gii';...
% '/data/cn4/evan/ROIs/Craddock/Craddock_350_L_removeedge.func.gii','/data/cn4/evan/ROIs/Craddock/Craddock_350_R_removeedge.func.gii';...
% '/data/cn4/evan/RestingState/FC_Mapping_120/Yeo_17_L_parcels_removeedge.func.gii','/data/cn4/evan/RestingState/FC_Mapping_120/Yeo_17_R_parcels_removeedge.func.gii';...
% '/data/cn4/evan/ROIs/mode.L.aparc.a2009s.32k_fs_LR_removeedge.func.gii','/data/cn4/evan/ROIs/mode.R.aparc.a2009s.32k_fs_LR_removeedge.func.gii'};

allparcelfilenames = {'/data/cn4/evan/RestingState/FC_Mapping_120/Cohorts/C1_Consensus_L_parcels.func.gii','/data/cn4/evan/RestingState/FC_Mapping_120/Cohorts/C1_Consensus_R_parcels.func.gii'};

iterations = 500;


    load('rotations.mat')

    load('/data/cn4/evan/RestingState/FC_Mapping_120/Homogeneity_testing_in_60/C2_cov_corr_L_normalwall.mat')
    for i = 1:size(allparcelfilenames,1)
            slashloc = strfind(allparcelfilenames{i,1},'/');
            outfileprefix = [allparcelfilenames{i,1}(slashloc(end)+1 : end-9)];
            
            generate_rotated_parcellation_andPCA_nobaddata(allparcelfilenames{i,1},iterations,cov_corr_L,rotations,1,'L',outfileprefix);
            
    end
    clear cov_corr_L
    
    
    load('/data/cn4/evan/RestingState/FC_Mapping_120/Homogeneity_testing_in_60/C2_cov_corr_R_normalwall.mat')
    for i = 1:size(allparcelfilenames,1)
        slashloc = strfind(allparcelfilenames{i,2},'/');
    outfileprefix = [allparcelfilenames{i,2}(slashloc(end)+1 : end-9)];
    generate_rotated_parcellation_andPCA_nobaddata(allparcelfilenames{i,2},iterations,cov_corr_R,rotations,1,'R',outfileprefix);
    end
    clear cov_corr_R


    load('rotations2.mat')

    load('/data/cn4/evan/RestingState/FC_Mapping_120/Homogeneity_testing_in_60/C2_cov_corr_L_normalwall.mat')
    for i = 1:size(allparcelfilenames,1)
            slashloc = strfind(allparcelfilenames{i,1},'/');
            outfileprefix = [allparcelfilenames{i,1}(slashloc(end)+1 : end-9) '2'];
            
            generate_rotated_parcellation_andPCA_nobaddata(allparcelfilenames{i,1},iterations,cov_corr_L,rotations,1,'L',outfileprefix);
            
    end
    clear cov_corr_L
    
    
    load('/data/cn4/evan/RestingState/FC_Mapping_120/Homogeneity_testing_in_60/C2_cov_corr_R_normalwall.mat')
    for i = 1:size(allparcelfilenames,1)
        slashloc = strfind(allparcelfilenames{i,2},'/');
    outfileprefix = [allparcelfilenames{i,2}(slashloc(end)+1 : end-9) '2'];
    generate_rotated_parcellation_andPCA_nobaddata(allparcelfilenames{i,2},iterations,cov_corr_R,rotations,1,'R',outfileprefix);
    end
    clear cov_corr_R
