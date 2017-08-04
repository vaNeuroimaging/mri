function [maxoutputmetric meanoutputmetric]= parcel_similarity(watershed,avgcrosscorrname,iscifti,hem,neighdist)
%[maxoutputmetric meanoutputmetric]= parcel_similarity(watershed,avgcrosscorrname,iscifti,hem,neighdist)

%  watershed = data.cdata;
%  avgcrosscorrname = avgcorr.cdata;
%  iscifti = 2;
%  hem = 'L';
%  neighdist = 3;

if ischar(watershed)
    watershed = gifti(watershed); watershed = watershed.cdata;
end

bufsize=16384;
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;
%neighbors(isnan(neighbors)) = 617;

origwatershed = watershed;

maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii'];
mask = gifti(maskname);
mask = mask.cdata;
corticalindices = (mask==0);


    if iscifti == 1
        watershed = watershed(mask==0,:);
    elseif iscifti == 2
        maskname = ['/data/cn4/laumannt/subcortical_mask/' hem '.atlasroi_erode3.32k_fs_LR.shape.gii'];
        mask = gifti(maskname);
        mask = ~mask.cdata;
        watershed = watershed(mask==0,:);
    end

%corrmat = load([corrdir '/cifti_avgcrosscorr.mat']);
if ~ischar(avgcrosscorrname)
    avgcrosscorr = avgcrosscorrname;
    clear avgcrosscorrname

else
    
    try
        
        avgcrosscorr = gifti(avgcrosscorrname);
        avgcrosscorr = avgcrosscorr.cdata;
        
    catch
        avgcrosscorr = cifti_read(avgcrosscorrname);
    end
    
end

avgcrosscorr(isnan(avgcrosscorr)) = 0;

parcelIDs = unique(watershed); parcelIDs(parcelIDs==0) = [];

parcelcorrelpatterns = zeros(length(parcelIDs),size(avgcrosscorr,2));

for parcelnum = 1:length(parcelIDs)
    parcelneighbors{parcelnum} = [];
    parcelindices{parcelnum} = find(watershed==parcelIDs(parcelnum));
    parcelcorrelpatterns(parcelnum,:) = mean(avgcrosscorr(parcelindices{parcelnum},:),1);
end

zeroindices = setdiff(find(origwatershed==0),find(mask));
for index = zeroindices'
    
    nodeneigh = index;
    newneigh = nodeneigh;
    curneigh = newneigh;
    
    for n = 1:neighdist
        
        for t = 1:length(curneigh)
            newneigh = [newneigh neighbors(curneigh(t),2:7)];
            newneigh(isnan(newneigh)) = [];
        end
        curneigh = setdiff(newneigh,nodeneigh);
        nodeneigh = union(nodeneigh,newneigh);
        
    end
    
    parcelsnearby = unique(origwatershed(nodeneigh)); parcelsnearby(parcelsnearby==0) = [];
    
    for parcelnum = 1:length(parcelsnearby)
        parcelindicesnearby(parcelnum) = find(parcelIDs==parcelsnearby(parcelnum));
    end
    for parcelnum = 1:length(parcelsnearby)
        parcelneighbors{parcelindicesnearby(parcelnum)} = union(parcelneighbors{parcelindicesnearby(parcelnum)},parcelindicesnearby);
    end
    clear parcelindicesnearby
end

maxoutputmetric = zeros(size(origwatershed));
meanoutputmetric = zeros(size(origwatershed));

for parcelnum = 1:length(parcelIDs)
    
    parcelneighbors{parcelnum}(parcelneighbors{parcelnum}==parcelnum) = [];
    
    correlations = paircorr_mod(parcelcorrelpatterns(parcelnum,:)',parcelcorrelpatterns(parcelneighbors{parcelnum},:)');
    
    maxoutputmetric(origwatershed==parcelIDs(parcelnum)) = max(correlations);
    meanoutputmetric(origwatershed==parcelIDs(parcelnum)) = mean(correlations);
    
    
end
