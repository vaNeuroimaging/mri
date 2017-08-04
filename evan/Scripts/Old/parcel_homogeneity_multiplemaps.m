threshs = [1.01 1.05 1.1 1.15 1.2 1.25 1.3 1.35 1.4 1.45 1.5 1.6 1.7 1.8 1.9 2];

hem = 'L';

maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii'];
mask = gifti(maskname);
mask = mask.cdata;
corticalindices = find(mask==0);

gooddataname = ['/data/cn4/evan/ROIs/all_meanimage_' hem '.func.gii'];
gooddata = gifti(gooddataname); gooddata = gooddata.cdata;
gooddata = gooddata>750;

for threshnum = 1:length(threshs)
    thresh = threshs(threshnum);
    watermap = gifti(['AdjustedWatershedEdges_' num2str(thresh) '_watershedmerge.func.gii']); watermap = watermap.cdata;
    PCAmap = gifti(['PCA_eigval_per_first_AdjustedWatershedEdges_' num2str(thresh) '_watershedmerge.func.gii']); PCAmap = PCAmap.cdata;
    
    watermap = watermap .* (mask==0);
    
    watersheds = unique(watermap); watersheds(watersheds==0) = [];
    
    for waternum = 1:length(watersheds)
        if any((~gooddata) .* (watermap==watersheds(waternum)))
            watermap(watermap==watersheds(waternum)) = 0;
        end
    end
    
    watersheds = unique(watermap); watersheds(watersheds==0) = [];
    
    for waternum = 1:length(watersheds)
        watershed = watersheds(waternum);
        
        allPCAvals{threshnum}(waternum) = mean(PCAmap(watermap==watershed));
        allsizes{threshnum}(waternum) = nnz(watermap==watershed);
        
    end
    
    meanPCAvals(threshnum) = mean(allPCAvals{threshnum});
end
