function create_parcel_center_rois(parcelsLfile,parcelsRfile,radius)

% distances = smartload('/data/cn4/evan/fsaverage_LR32k/Surface_distances_L.mat');
% 
% parcels = gifti(parcelsLfile); parcels = parcels.cdata;
% 
% IDs = unique(parcels); IDs(IDs==0) = [];
% 
% centroids = zeros(length(IDs),1);
% for IDnum = 1:length(IDs)
%     parcelinds = find(parcels==IDs(IDnum));
%     [ign mini] = min(sum(distances(parcelinds,parcelinds),2));
%     centroids(IDnum) = parcelinds(mini) - 1;
% end
% 
% dlmwrite('ROICentersL.txt',centroids)
% 
% system(['wb_command -surface-geodesic-rois /data/cn4/evan/fsaverage_LR32k/Conte69.L.midthickness.32k_fs_LR.surf.gii ' num2str(radius) ' ROICentersL.txt ParcelCenterROIs_' num2str(radius) 'mm_L_separated.func.gii -overlap-logic CLOSEST'])
% 
% ROIs_all = gifti(['ParcelCenterROIs_' num2str(radius) 'mm_L_separated.func.gii']); ROIs_all = ROIs_all.cdata;
% ROIs = zeros(length(ROIs_all),1);
% for i = 1:size(ROIs_all,2)
%     ROIs(logical(ROIs_all(:,i))) = IDs(i);
% end
% mask = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii']);
% ROIs(logical(mask.cdata)) = 0;
% save(gifti(single(ROIs)),['ParcelCenterROIs_' num2str(radius) 'mm_L.func.gii'])






distances = smartload('/data/cn4/evan/fsaverage_LR32k/Surface_distances_R.mat');

parcels = gifti(parcelsRfile); parcels = parcels.cdata;

IDs = unique(parcels); IDs(IDs==0) = [];

centroids = zeros(length(IDs),1);
for IDnum = 1:length(IDs)
    parcelinds = find(parcels==IDs(IDnum));
    [ign mini] = min(sum(distances(parcelinds,parcelinds),2));
    centroids(IDnum) = parcelinds(mini) - 1;
end

dlmwrite('ROICentersR.txt',centroids)

system(['wb_command -surface-geodesic-rois /data/cn4/evan/fsaverage_LR32k/Conte69.R.midthickness.32k_fs_LR.surf.gii ' num2str(radius) ' ROICentersR.txt ParcelCenterROIs_' num2str(radius) 'mm_R_separated.func.gii -overlap-logic CLOSEST'])

ROIs_all = gifti(['ParcelCenterROIs_' num2str(radius) 'mm_R_separated.func.gii']); ROIs_all = ROIs_all.cdata;
ROIs = zeros(length(ROIs_all),1);
for i = 1:size(ROIs_all,2)
    ROIs(logical(ROIs_all(:,i))) = IDs(i);
end
mask = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.R.32k_fs_LR.func.gii']);
ROIs(logical(mask.cdata)) = 0;
save(gifti(single(ROIs)),['ParcelCenterROIs_' num2str(radius) 'mm_R.func.gii'])
