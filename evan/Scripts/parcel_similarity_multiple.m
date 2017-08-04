function multioutputmetric = parcel_similarity_multiple(multiwatershed,avgcrosscorrname,iscifti,hem)
%outputmetric = parcel_similarity(watershedname,avgcrosscorrname,iscifti,outputname,hem)



bufsize=16384;
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;



maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii'];
mask = gifti(maskname);
mask = mask.cdata;
corticalindices = find(mask==0);



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

avgcrosscorr(isnan(avgcrosscorr)) = 0;

multioutputmetric = zeros(size(mask,1),size(multiwatershed,2));


for map = 1:size(multiwatershed,2)
    
    watershed = multiwatershed(:,map);
    
    if iscifti
        ciftispace_watershed = watershed(corticalindices);
    else
        ciftispace_watershed = watershed;
    end
    
    waternum = unique(watershed);
    waternum(waternum==0) = [];
    
    for waterindex = 1:length(waternum)
        thiswaternum = waternum(waterindex);
        ciftiindices{waterindex} = find(ciftispace_watershed==thiswaternum);
        watershed_correlpattern(waterindex,:) = mean(avgcrosscorr(ciftiindices{waterindex},:),1);
        adjacentwatersheds{waterindex} = [];
    end
    
    
    borderindices = intersect(find(watershed==0),corticalindices)';
    
    for bordervertex = borderindices
        
        borderneighs = neighbors(bordervertex,2:7);
        borderneighs(isnan(borderneighs)) = [];
        borderneighvals = watershed(borderneighs);
        borderneighvals(borderneighvals==0) = [];
        
        watershedneighbors = unique(borderneighvals);
        
        for waterneighbornum = 1:length(watershedneighbors)
            
            thiswatershedindex = find(waternum==watershedneighbors(waterneighbornum));
            
            otherwatershedneighbors = watershedneighbors;
            otherwatershedneighbors(otherwatershedneighbors==watershedneighbors(waterneighbornum)) = [];
            
            adjacentwatersheds{thiswatershedindex} = unique([adjacentwatersheds{thiswatershedindex} otherwatershedneighbors']);
        end
    end
    
    outputmetric = zeros(size(mask));
    
    for waterindex = 1:length(waternum)
        
        thiswaternum = waternum(waterindex);
        
        if ~isempty(adjacentwatersheds{waterindex}) > 0
            
            for neighborparcelindex = 1:length(adjacentwatersheds{waterindex})
                
                adjacentwatershedindex = find(waternum==adjacentwatersheds{waterindex}(neighborparcelindex));
                
                neighborpattern = watershed_correlpattern(adjacentwatershedindex,:);
                
                nonwatershedindices = setdiff(corticalindices,ciftiindices{waterindex});
                nonwatershedindices = setdiff(nonwatershedindices,ciftiindices{adjacentwatershedindex});
                
                similarities(neighborparcelindex) = paircorr_mod(watershed_correlpattern(waterindex,nonwatershedindices)',neighborpattern(nonwatershedindices)');
                
            end
            maxsimilarity = max(similarities);
            
            outputmetric(watershed==thiswaternum) = maxsimilarity;
            
        end
        
        clear similarities
        
    end
    
    multioutputmetric(:,map) = outputmetric .* 100;
    
end

%save(gifti(single(outputmetric)),outputname);


