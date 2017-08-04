function multioutputmetric = parcel_homogeneity_PCA_multiple(multiwatershed,avgcrosscorrname,iscifti,hem)
%outputmetric = parcel_homogeneity_PCA(watershedname,avgcrosscorrname,iscifti,outputname,hem)

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

    if iscifti
        multiwatershed = multiwatershed(mask==0,:);
    end

%corrmat = load([corrdir '/cifti_avgcrosscorr.mat']);
if ~ischar(avgcrosscorrname)
    avgcrosscorr = avgcrosscorrname;
    clear avgcrosscorrname
elseif iscifti
    avgcrosscorr = cifti_read(avgcrosscorrname);
else
    
    try
        
        avgcrosscorr = gifti(avgcrosscorrname);
        avgcrosscorr = avgcrosscorr.cdata;
        
    catch
        avgcrosscorr = cifti_read(avgcrosscorrname);
    end
    
end

if iscifti==2
    erodedmask = gifti('/data/cn4/laumannt/subcortical_mask/L.atlasroi_erode3.32k_fs_LR.shape.gii'); erodedmask = erodedmask.cdata;
    [ign indices_within_eroded ign2] = intersect(find(erodedmask),find(mask==0));
    avgcrosscorr = avgcrosscorr(indices_within_eroded,:);
end

%avgcrosscorr(isnan(avgcrosscorr)) = 0;

%corrmat_use = corrmat.avgcrosscorr(1:29696,:);
%corrmat_use = corrmat.avgcrosscorr_reshape;





%% Calculate PCA, determine explained variance first component

for map = 1:size(multiwatershed,2)
    watershed = multiwatershed(:,map);
    

    %watershed(logical(mask)) = 0;
    waternum = unique(watershed);
    waternum(waternum==0) = [];
    
    clear Y eigval_per_first
    for w = 1:length(waternum)
        
        %waterval = watershed(logical(~mask))==waternum(w);
        waterval = watershed==waternum(w);
        %watercorr = corrmat_use(waterval,:);
        watercorr = avgcrosscorr(waterval,:);
        
        nanmat = isnan(watercorr);
        nancolsum = sum(nanmat,2);
        watercorr((nancolsum>0),:)= [];
        
        if size(watercorr,1) > 2
            
            [Y(:,:,w) V_use eigvals_per{w}] = PCA_reduction(watercorr','comps',2);
            eigval_per_first(w) = eigvals_per{w}(1);
            
        else
            eigval_per_first(w) = 0;
        end
        
    end
    
    if iscifti
        metric = zeros(nnz(mask==0),1);
    else
        metric = zeros(size(mask));
    end
    
    
    for w = 1:length(waternum)
        
        waterval = watershed==waternum(w);
        
        metric(waterval) = eigval_per_first(w);
    end
    
    outputmetric = zeros(length(mask),1);
    if iscifti
        outputmetric(mask==0) = metric;
    else
        outputmetric = metric;
    end
    
    multioutputmetric(:,map) = outputmetric;
end





