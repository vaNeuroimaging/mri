function outmatrix = parcel_physical_adjacency(parcels,distances,distthresh)

if ischar(parcels)
    parcels = ft_read_cifti_mod(parcels); parcels = parcels.data;
elseif isstruct(parcels)
    parcels = parcels.data;
end

if ischar(distances)
    distances = smartload(distances);
end

if ~exist('distthresh')
    distthresh = 1000000;
end

parcelIDs = unique(parcels); parcelIDs(parcelIDs<1) = [];

outmatrix = zeros(length(parcelIDs));

inds_withparcels = find(parcels>0);

for parcelnum = 1:length(parcelIDs)
    parcelID = parcelIDs(parcelnum);
    
    thisparcel_inds = find(parcels==parcelID);
    otherparcel_inds = setdiff(inds_withparcels,thisparcel_inds);
    
    [closest_distances,closest_temp] = min(distances(thisparcel_inds,otherparcel_inds),[],2);
    closest_inds = otherparcel_inds(closest_temp);
    closest_inds(closest_distances>distthresh) = [];
    closest_parcels = unique(parcels(closest_inds));
    for i = 1:length(closest_parcels)   
        outmatrix(parcelID,find(parcelIDs==closest_parcels(i))) = 1;
    end
end

