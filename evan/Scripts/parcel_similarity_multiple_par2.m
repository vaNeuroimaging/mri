function multioutputmetric = parcel_similarity_multiple_par2(multiwatershed,avgcrosscorrname,iscifti,hem)
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

watershed = zeros(numel(multiwatershed),1);


for map = 1:size(multiwatershed,2)
    watershed((map-1)*size(multiwatershed,1) +1 : map*size(multiwatershed,1),1) = multiwatershed(:,map) + max(watershed);
    if iscifti
        ciftispace_watershed((map-1)*length(corticalindices) +1 : map*length(corticalindices),1) = watershed(corticalindices);
    end
    allneighbors((map-1)*size(multiwatershed,1) +1 : map*size(multiwatershed,1),:) = neighbors + (map-1)*size(multiwatershed,1);
end

waternum = unique(watershed);
waternum(waternum==0) = [];

maxsimilarity = zeros(length(waternum),1);

for w = 1:length(waternum)
    disp([num2str(w) ' of ' num2str(length(waternum))])
    waterindices = find(watershed==waternum(w));
    
    if iscifti
        ciftiindices = find(ciftispace_watershed==waternum(w));
        brainindices = mod(ciftiindices,length(corticalindices));
        brainindices(brainindices==0) = length(corticalindices);
    else
        brainindices = mod(waterindices,size(multiwatershed,1));
        brainindices(brainindices==0) = size(multiwatershed,1);
    end
    
    
    watershed_correlpattern = mean(avgcrosscorr(brainindices,:),1);
    
    watershedvertexneighbors = unique(allneighbors(waterindices,2:7));
    watershedborders = intersect(watershedvertexneighbors,find(watershed==0));
    watershedborderneighbors = unique(allneighbors(watershedborders,2:7));
    adjacentwatershednums = unique(watershed(watershedborderneighbors));
    adjacentwatershednums(adjacentwatershednums==0) = [];
    adjacentwatershednums(adjacentwatershednums==waternum(w)) = [];
    
    for i = 1:length(adjacentwatershednums)
        if iscifti
            ciftiindices = find(ciftispace_watershed==adjacentwatershednums(i));
            brainindices = mod(ciftiindices,length(corticalindices));
            brainindices(brainindices==0) = length(corticalindices);
        else
            theseindices = find(watershed==adjacentwatershednums(i));
            brainindices = mod(theseindices,size(multiwatershed,1));
            brainindices(brainindices==0) = size(multiwatershed,1);
        end
        adjacentwatershed_correlpattern = mean(avgcrosscorr(brainindices,:),1);
 
        similarity = corr(watershed_correlpattern',adjacentwatershed_correlpattern');
        
        if similarity > maxsimilarity(w)
            maxsimilarity(w) = similarity;
        end
    end
end
%matlabpool close

multioutputmetric = zeros(size(multiwatershed));
for parcel = 1:length(waternum)
    waterindices = find(watershed==waternum(parcel));
    brainindices = mod(waterindices,size(multiwatershed,1));
    brainindices(brainindices==0) = size(multiwatershed,1);
    map = ceil(waterindices(1)/size(multiwatershed,1));
    multioutputmetric(brainindices{parcel},map) = maxsimilarity(parcel);
end




%save(gifti(single(outputmetric)),outputname);


