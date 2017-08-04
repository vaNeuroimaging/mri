% allparcelfilenames = {'/data/cn4/evan/RestingState/FC_Mapping_120/Cohorts/C1_Lwatershedmerge_0.45.func.gii','/data/cn4/evan/RestingState/FC_Mapping_120/Cohorts/C1_Rwatershedmerge_0.45.func.gii';...
% '/data/cn4/evan/ROIs/264_surfvert_ROIs_L.func.gii','/data/cn4/evan/ROIs/264_surfvert_ROIs_R.func.gii';...
% '/data/cn4/evan/RestingState/FC_Mapping_120/ConsensusParcels_L.func.gii','/data/cn4/evan/RestingState/FC_Mapping_120/ConsensusParcels_R.func.gii';...
% '/data/cn4/evan/ROIs/Craddock/Craddock_350_L.func.gii','/data/cn4/evan/ROIs/Craddock/Craddock_350_R.func.gii';...
% '/data/cn4/evan/RestingState/FC_Mapping_120/Yeo_17_L_parcels.func.gii','/data/cn4/evan/RestingState/FC_Mapping_120/Yeo_17_R_parcels.func.gii';...
% '/data/cn4/evan/ROIs/mode.L.aparc.a2009s.32k_fs_LR.func.gii','/data/cn4/evan/ROIs/mode.R.aparc.a2009s.32k_fs_LR.func.gii'};

allparcelfilenames = {'/data/cn4/evan/RestingState/FC_Mapping_120/120_L_wateredgethresh_watershedmerge_0.45.func.gii','/data/cn4/evan/RestingState/FC_Mapping_120/120_R_wateredgethresh_watershedmerge_0.45.func.gii';...
    '/data/cn4/evan/RestingState/FC_Mapping_120/120_L_wateredgethresh_watershedmerge_0.25.func.gii','/data/cn4/evan/RestingState/FC_Mapping_120/120_R_wateredgethresh_watershedmerge_0.25.func.gii';...
    '/data/cn4/evan/RestingState/FC_Mapping_120/120_L_wateredgethresh_watershedmerge_0.65.func.gii','/data/cn4/evan/RestingState/FC_Mapping_120/120_R_wateredgethresh_watershedmerge_0.65.func.gii'};

iterations = 1000;



    
    for i = 1:size(allparcelfilenames,1)
        
        slashloc = strfind(allparcelfilenames{i,1},'/');
        outfileprefix = [allparcelfilenames{i,1}(slashloc(end)+1 : end-9)];
        
        if ~exist([outfileprefix '.mat'])
            load('/data/cn4/evan/RestingState/FC_Mapping_120/cov_corr_normalwall_L.mat');
        
            %load('/data/cn4/evan/RestingState/FC_Mapping_120/Homogeneity_testing_in_60/C2_cov_corr_L_normalwall.mat')
            
            generate_rotated_parcels_andPCA7(allparcelfilenames{i,1},iterations,cov_corr_L,1,'L',outfileprefix);
            clear cov_corr_L
        end
            
        slashloc = strfind(allparcelfilenames{i,2},'/');
        outfileprefix = [allparcelfilenames{i,2}(slashloc(end)+1 : end-9)];
        
        if ~exist([outfileprefix '.mat'])
            load('/data/cn4/evan/RestingState/FC_Mapping_120/cov_corr_normalwall_R.mat');
        
            %load('/data/cn4/evan/RestingState/FC_Mapping_120/Homogeneity_testing_in_60/C2_cov_corr_R_normalwall.mat')
            
            generate_rotated_parcels_andPCA7(allparcelfilenames{i,2},iterations,cov_corr_R,1,'R',outfileprefix);
            clear cov_corr_R
        end
    end
    

    
