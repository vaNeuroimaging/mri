function outputmetric = parcel_homogeneity_PCA_vector(watershed,covmat,hem)

iscifti = 0;

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
   
%     if ~ischar(watershedname)
%     
%     watershed = gifti([watershedname]);
%     watershed = watershed.cdata;
%     else
%         watershed = watershedname;
%     end
    
    maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii'];
    mask = gifti(maskname);
    mask = mask.cdata;

    %watershed(logical(mask)) = 0;
    waternum = unique(watershed);
    waternum(waternum==0) = [];
    
    
    
    
    
    
    %% Calculate PCA, determine explained variance first component
    
    clear Y eigval_per_first
    for w = 1:length(waternum)
        
        %waterval = watershed(logical(~mask))==waternum(w);
        waterval = watershed==waternum(w);
        %watercorr = corrmat_use(waterval,:);
        watercovmat = covmat(waterval,waterval);
                
        nanmat = isnan(watercovmat);
        nancolsum = sum(nanmat,2);
        watercovmat((nancolsum>0),:)= [];
        
        if size(watercovmat,1) > 2
        
        [ign ign2 eigvals_per{w}] = PCA_reduction_cov(watercovmat,'comps',2);
        eigval_per_first(w) = eigvals_per{w}(1);
        
        else
            eigval_per_first(w) = 0;
        end
        
    end
    
    %%
    %metric = zeros(length(find(mask==0)),1);
    metric = zeros(size(mask));
    
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
    
    %save(gifti(single(outputmetric)),outputname)
    
    %%
    
%     for w = 1:length(waternum)
%         
%         waterval = watershed==waternum(w);
%         
%         waterval_num(w) = nnz(waterval);
%     end
%     
%     
%     corr(waterval_num',eigval_per_first')
    
%end

%gifti_to_cifti([waterdir '/PCA_eigval_per_first_water_' watershedname{1}],[waterdir '/PCA_eigval_per_first_water_' watershedname{2}],outputname);





