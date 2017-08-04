allparcelfilenames = {'/data/cn4/evan/RestingState/FC_Mapping_120/120_L_wateredgethresh_watershedmerge_0.45.func.gii','/data/cn4/evan/RestingState/FC_Mapping_120/120_R_wateredgethresh_watershedmerge_0.45.func.gii';...
'/data/cn4/evan/ROIs/264_surfvert_ROIs_L.func.gii','/data/cn4/evan/ROIs/264_surfvert_ROIs_R.func.gii';...
'/data/cn4/evan/ROIs/Craddock/Craddock_350_L.func.gii','/data/cn4/evan/ROIs/Craddock/Craddock_350_R.func.gii';...
'/data/cn4/evan/ROIs/Shen/Shen_300_L.func.gii','/data/cn4/evan/ROIs/Shen/Shen_300_R.func.gii';...
'/data/cn4/evan/RestingState/FC_Mapping_120/Yeo_17_L_parcels.func.gii','/data/cn4/evan/RestingState/FC_Mapping_120/Yeo_17_R_parcels.func.gii';...
'/data/cn4/evan/ROIs/mode.L.aparc.a2009s.32k_fs_LR.func.gii','/data/cn4/evan/ROIs/mode.R.aparc.a2009s.32k_fs_LR.func.gii';...
%'/data/cn4/evan/RestingState/FC_Mapping_120/ConsensusParcels_L.func.gii','/data/cn4/evan/RestingState/FC_Mapping_120/ConsensusParcels_R.func.gii';...
};

% allparcelfilenames = {'/data/cn4/evan/ROIs/264_surfvert_ROIs_L_removeedge.func.gii','/data/cn4/evan/ROIs/264_surfvert_ROIs_R_removeedge.func.gii';...
%     '/data/cn4/evan/RestingState/FC_Mapping_120/ConsensusParcels_L_removeedge.func.gii','/data/cn4/evan/RestingState/FC_Mapping_120/ConsensusParcels_R_removeedge.func.gii';...
% '/data/cn4/evan/ROIs/Craddock/Craddock_350_L_removeedge.func.gii','/data/cn4/evan/ROIs/Craddock/Craddock_350_R_removeedge.func.gii';...
% '/data/cn4/evan/RestingState/FC_Mapping_120/Yeo_17_L_parcels_removeedge.func.gii','/data/cn4/evan/RestingState/FC_Mapping_120/Yeo_17_R_parcels_removeedge.func.gii';...
% '/data/cn4/evan/ROIs/mode.L.aparc.a2009s.32k_fs_LR_removeedge.func.gii','/data/cn4/evan/ROIs/mode.R.aparc.a2009s.32k_fs_LR_removeedge.func.gii'};

%allparcelfilenames = {'/data/cn4/evan/RestingState/FC_Mapping_120/Cohorts/C1_Consensus_L_parcels.func.gii','/data/cn4/evan/RestingState/FC_Mapping_120/Cohorts/C1_Consensus_R_parcels.func.gii'};

% for i = 4:42
%     
%     allparcelfilenames(i-3,:) = {['/data/cn4/evan/ROIs/Craddock/Craddock' sprintf('%04i',i) '_L.func.gii'],['/data/cn4/evan/ROIs/Craddock/Craddock' sprintf('%04i',i) '_R.func.gii']};
% end
% allparcelfilenames(end+1,:) = {'/data/cn4/evan/ROIs/Shen/Shen_200_L.func.gii','/data/cn4/evan/ROIs/Shen/Shen_200_R.func.gii'};
% allparcelfilenames(end+1,:) = {'/data/cn4/evan/ROIs/Shen/Shen_100_L.func.gii','/data/cn4/evan/ROIs/Shen/Shen_100_R.func.gii'};

mkdir ReplaceNan

for i = 1:size(allparcelfilenames,1)
    for j = 1:2
        slashloc = strfind(allparcelfilenames{i,j},'/');
        outfileprefix = [allparcelfilenames{i,j}(slashloc(end)+1 : end-9)];
        
        load([outfileprefix '.mat'])
        
        rot1 = rotated_eigvals;
        
        load([outfileprefix '2.mat'])
        
        rot1 = [rot1 rotated_eigvals];
        
        load([outfileprefix '3.mat'])
        
        rot1 = [rot1 rotated_eigvals];
        
        load([outfileprefix '4.mat'])
        
        rotated_eigvals = [rot1 rotated_eigvals];
        
%         save([outfileprefix '.mat'],'rotated_eigvals','realsizes','realeigvals_per_first')
%         
%     end
%     
% end
% 
% 
% 
% for i = 1:size(allparcelfilenames,1)
%     for j = 1:2
%         slashloc = strfind(allparcelfilenames{i,j},'/');
%         outfileprefix = [allparcelfilenames{i,j}(slashloc(end)+1 : end-9)];
%        
%        load([outfileprefix '.mat'])
        
        for k = 1:size(rotated_eigvals,1)
            naninds = logical(isnan(rotated_eigvals(k,:)));
            rotated_eigvals(k,naninds) = nanmean(rotated_eigvals(k,:));
        end
        
        
        save(['ReplaceNan/' outfileprefix '.mat'],'rotated_eigvals','realsizes','realeigvals_per_first')
        
    end
    
end
