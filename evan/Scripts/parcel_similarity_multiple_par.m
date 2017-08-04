function multioutputmetric = parcel_similarity_multiple_par(multiwatershed,avgcrosscorrname,iscifti,hem)
%outputmetric = parcel_similarity(watershedname,avgcrosscorrname,iscifti,outputname,hem)



bufsize=16384;
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;
neighbors(isnan(neighbors)) = 617;



maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii'];
mask = gifti(maskname);
mask = mask.cdata;
corticalindices = (mask==0);



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


index = 1;
for map = 1:size(multiwatershed,2)
    disp(num2str(map))
    watershed = multiwatershed(:,map);
    waternum = unique(watershed);
    waternum(waternum==0) = [];
    
    if iscifti
        ciftispace_watershed = watershed(corticalindices);
    else
        ciftispace_watershed = watershed;
    end
    
    for w = 1:length(waternum)
        brainindices{index} = find(watershed==waternum(w));
        mapnum(index) = map;
        thismapwatershedindex(w) = index;
                
        ciftiindices = ciftispace_watershed==waternum(w);
        watershed_correlpattern(index,:) = mean(avgcrosscorr(ciftiindices,:),1);
                
        watershedvertexneighbors = unique(neighbors(brainindices{index},2:7));
        watershedborders = intersect(watershedvertexneighbors,find(watershed==0));
        watershedborderneighbors = unique(neighbors(watershedborders,2:7));
        adjacentwatershednums = unique(watershed(watershedborderneighbors));
        adjacentwatershednums(adjacentwatershednums==0) = [];
        adjacentwatershednums(adjacentwatershednums==waternum(w)) = [];
        for i = 1:length(adjacentwatershednums)
            adjacentwatershedindices{w}(i) = find(waternum==adjacentwatershednums(i));
        end
        index = index+1;
    end
    
    for w = 1:length(waternum)
        adjacentwatersheds{thismapwatershedindex(w)} = thismapwatershedindex(adjacentwatershedindices{w});
    end
    
    clear thismapwatershedindex adjacentwatershedindices
        
     
end
 
multioutputmetric = zeros(size(mask,1),size(multiwatershed,2));


%matlabpool open 2
for parcel = 1:length(brainindices)
    disp(num2str(parcel))
        
    if ~isempty(adjacentwatersheds{parcel})
    
        similarities = paircorr_mod(watershed_correlpattern(parcel,:)', watershed_correlpattern(adjacentwatersheds{parcel},:)');
        maxsimilarity(parcel) = max(similarities) .* 100;
        
                
    end
        
end

%matlabpool close
for parcel = 1:length(brainindices)
    multioutputmetric(brainindices{parcel},mapnum(parcel)) = maxsimilarity(parcel);
end




%save(gifti(single(outputmetric)),outputname);


