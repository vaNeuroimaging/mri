%distances = smartload('/data/cn4/evan/fsaverage_LR32k/Cifti_geo_distances.mat');

parcels = cifti_read('/data/cn4/evan/Published_parcels/Parcels_LR.dtseries.nii');

coords = cifti_read('/data/cn4/evan/fsaverage_LR32k/Conte69.LR.midthickness.32k_fs_LR_ASCIIformat_coords.dtseries.nii');

IDs = unique(parcels); IDs(IDs==0) = [];

centroidverts = zeros(length(IDs),1);
centroidcoords = zeros(length(IDs),3);
centroidmap = zeros(size(parcels));

for IDnum = 1:length(IDs)
    ID = IDs(IDnum);
    
    parcelinds = find(parcels==ID);
    
    parcel_totaldistances = sum(distances(parcelinds,parcelinds),2);
    
    [ign mini] = min(parcel_totaldistances);
    
    centroidverts(IDnum) = parcelinds(mini);
    centroidcoords(IDnum,:) = coords(parcelinds(mini),:);
    centroidmap(parcelinds(mini)) = 1;
end

centroidcoords = round(centroidcoords*10) / 10;

%cifti_write_wHDR(centroidmap,[],'Parcel_centroids');
dlmwrite('Parcel_centroids.txt',centroidcoords,' ')
