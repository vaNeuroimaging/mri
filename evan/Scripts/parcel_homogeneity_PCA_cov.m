function multioutputmetric = parcel_homogeneity_PCA_cov(multiwatershed,cov_corr,iscifti,hem)
%outputmetric = parcel_homogeneity_PCA_cov(watershed,cov_corr,iscifti,hem)

%corrdir = '/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/modified_cifti_network/BOTH';
%corrdir = '/data/cn4/evan/Temp/';
%'/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/cifti_timeseries/';
%waterdir = '/data/cn4/evan/RestingState/Ind_variability/Poldrome/';
%outputname = [waterdir '/PCA_eigval_per_first_watermerge_Poldrome_LR'];




%%
%HEMS = {'L';'R'};
%hemname = {'LEFT';'RIGHT'};
%waterdir = '/data/cn4/laumannt/watershed_network/FCPROCESS_NEW/';

%for hem = 1%:2
%watershed = gifti([waterdir '/avgcorrofcorr_smooth255_allgrad_' HEMS{hem} '_smooth255_edge_avg_uc_smooth_minima5_iter200_frach1_watershed.func.gii']);
%watershed = gifti([waterdir '/avgcorrofcorr_smooth2.55_allgrad_' HEMS{hem} '_smooth2.55_edge_avg_minima3_watershedmerged.func.gii']);
%watershedname{hem} = ['Poldrome_' HEMS{hem} '_smoothed_testingwatershedmerge.func.gii'];
%watershed = gifti([waterdir watershedname{hem}]);



maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii'];
mask = gifti(maskname);
mask = mask.cdata;

    if iscifti == 1
        multiwatershed = multiwatershed(mask==0,:);
    elseif iscifti == 2
        maskname = ['/data/cn4/laumannt/subcortical_mask/' hem '.atlasroi_erode3.32k_fs_LR.shape.gii'];
        mask = gifti(maskname);
        mask = ~mask.cdata;
        multiwatershed = multiwatershed(mask==0,:);
    end



%corrmat_use = corrmat.avgcrosscorr(1:29696,:);
%corrmat_use = corrmat.avgcrosscorr_reshape;

index = 1;
for map = 1:size(multiwatershed,2)
    watershed = multiwatershed(:,map);
    waternum = unique(watershed);
    waternum(waternum==0) = [];
    for w = 1:length(waternum)
        brainindices{index} = find(watershed==waternum(w));
        mapnum(index) = map;
        eigval_per_first{index} = 0;
        index = index+1;
    end
end

if iscifti
        multioutputmetric = zeros(nnz(mask==0),size(multiwatershed,2));
else
        multioutputmetric = zeros(length(mask),size(multiwatershed,2));
end

%fprintf(repmat(' ',1,30))

%% Calculate PCA, determine explained variance first component
%matlabpool open 2
%parfor parcel = 1:size(brainindices,2)
for parcel = 1:length(brainindices)
    
%     if rem(parcel,50)==0 || parcel == length(brainindices)
%         fprintf(repmat('\b',1,30))
%         fprintf('%30s',['Parcel ' num2str(parcel) ' of ' num2str(length(brainindices))])
%     end
    
        watercovcorr = cov_corr(brainindices{parcel},brainindices{parcel});
        
        
        if size(watercovcorr,1) > 2
            
            eigvals_per = PCA_reduction_cov_onecomp(watercovcorr);
            eigval_per_first{parcel} = eigvals_per(1);
            
%         else
%             eigval_per_first{parcel} = 0;
        end
    
end
%disp(' ')

%matlabpool close

for parcel = 1:length(eigval_per_first)
    multioutputmetric(brainindices{parcel},mapnum(parcel)) = eigval_per_first{parcel};
end
    
if iscifti
    temp = multioutputmetric;
    multioutputmetric = zeros(length(mask),size(multioutputmetric,2));
    multioutputmetric(mask==0,:) = temp;
end
