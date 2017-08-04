function Hub_detection_summedPC(parcelsfile,corrmatfile,dmatfile,communitiesfile,consensus_modulesfile,kdenthresholds,xdist)
%Hub_detection_summedPC(parcelsfile,corrmatfile,dmatfile,communitiesfile,consensus_modulesfile,kdenthresh,xdist)

PC_cutoffs = [.5 : .05 : .95 .98];
pct_kdenthreshes_tocount_for_connection = .75;
pctconnections_needed_tocount_foracommunity = 0;%.1;

parcels = ft_read_cifti_mod(parcelsfile);
cifti_template = parcels;
parcels = parcels.data;
parcelIDs = unique(parcels); parcelIDs(parcelIDs<1) = [];

corrmat = smartload(corrmatfile);
dmat = smartload(dmatfile);
communities = load(communitiesfile);

corrmat(dmat<xdist) = 0;
triumask = triu(true(length(corrmat)));
sorted = sort(corrmat(triumask),'descend');

corrmat_bin_all = zeros(size(corrmat,1),size(corrmat,1),length(kdenthresholds));

PCvals = zeros(size(communities));

for threshnum = 1:length(kdenthresholds)
    kdenthresh = kdenthresholds(threshnum);
    rthresh = sorted(ceil(length(sorted) .* kdenthresh));
    corrmat_bin = corrmat > rthresh;
    corrmat_bin_all(:,:,threshnum) = corrmat_bin;
    
    thresh_communities = unique(communities(:,threshnum));
    thresh_communities(thresh_communities<1) = [];
    
    parcel_degree = sum(corrmat_bin,2);
    moduleconnectionssum = zeros(size(corrmat_bin,1),1);
    
    for communitynum = 1:length(thresh_communities)
        communityID = thresh_communities(communitynum);
        communityindices = communities(:,threshnum)==communityID;
        parcel_community_degree = sum(corrmat_bin(:,communityindices),2);
        parcel_community_ratio = (parcel_community_degree ./ parcel_degree);
        moduleconnectionssum = moduleconnectionssum + (parcel_community_ratio.^2);
    end
    
    PCvals(:,threshnum) = 1-moduleconnectionssum;
end
summed_PCs = nanmean(PCvals,2);
summed_PCs(isnan(summed_PCs)) = 0;

modules = ft_read_cifti_mod(consensus_modulesfile); modules = modules.data;
parcelmodules = zeros(size(corrmat,1),1);
for parcelIDnum = 1:length(parcelIDs)
    parcelmodules(parcelIDnum) = mode(modules(parcels==parcelIDs(parcelIDnum)));
end
moduleIDs = unique(parcelmodules); moduleIDs(moduleIDs<1) = [];

corrmat_bin_sum = sum(corrmat_bin_all,3);
corrmat_bin_sum_thresh = corrmat_bin_sum >= (length(kdenthresholds).*pct_kdenthreshes_tocount_for_connection);
parcel_connections = zeros(size(corrmat,1),length(moduleIDs));
for parcelIDnum = 1:length(parcelIDs)
    this_parcel_allconnections = unique(parcelmodules(corrmat_bin_sum_thresh(parcelIDnum,:)));
    this_parcel_allconnections(this_parcel_allconnections<1) = [];
    this_parcel_connections = [];
    for ID = this_parcel_allconnections(:)'
        if (nnz(parcelmodules(corrmat_bin_sum_thresh(parcelIDnum,:))==ID) / nnz(corrmat_bin_sum_thresh(parcelIDnum,:))) > pctconnections_needed_tocount_foracommunity
            this_parcel_connections(end+1) = ID;
        end
    end
    
    
    %this_parcel_connections = unique(parcelmodules(corrmat_bin_sum_thresh(parcelIDnum,:)));
    %this_parcel_connections(this_parcel_connections<1) = [];
    this_parcel_connections((end+1):length(moduleIDs)) = 0;
    parcel_connections(parcelIDnum,:) = this_parcel_connections;
end



save('summed_PCs.mat','summed_PCs')
save('parcel_connections.mat','parcel_connections')

[~,sortedindex] = sort(summed_PCs);
ordered(sortedindex) = 1:length(summed_PCs);
summed_PCs_percentile = ordered ./ length(summed_PCs);

summed_PCs_map = zeros(size(parcels));
parcel_connections_map = zeros(size(parcels,1),size(parcel_connections,2));
for parcelnum = 1:length(parcelIDs)
    summed_PCs_map(parcels==parcelIDs(parcelnum)) = summed_PCs_percentile(parcelnum);
    parcel_connections_map(parcels==parcelIDs(parcelnum),:) = repmat(parcel_connections(parcelnum,:),nnz(parcels==parcelIDs(parcelnum)),1);
end

cifti_template.data = summed_PCs_map;
ft_write_cifti_mod('parcel_PCs_percentile',cifti_template)  

cifti_template.data = parcel_connections_map;
ft_write_cifti_mod('parcel_connections',cifti_template)  

parcel_connectorhubs_map = zeros(length(parcels),length(PC_cutoffs));

for i = 1:length(PC_cutoffs)
    
    parcel_connectorhubs_map(:,i) = (summed_PCs_map >= PC_cutoffs(i))  .* (parcels > 0);
    cifti_template.mapname{i} = ['PC percentile thresh=' num2str(PC_cutoffs(i))];
    
end

cifti_template.data = parcel_connectorhubs_map;
cifti_template.dimord = 'scalar_pos';
ft_write_cifti_mod('parcel_ConnectorHubs',cifti_template)


evalc(['!wb_command -cifti-resample ' parcelsfile ' COLUMN /home/data/scripts/Resources/Conte69_atlas.LR.164k_fs_LR/Conte69.LR.164k_fs_LR.normalwall_surfaceonly_template.dtseries.nii COLUMN BARYCENTRIC ENCLOSING_VOXEL Temp_164.dtseries.nii -surface-largest -left-spheres /home/data/scripts/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.L.sphere.32k_fs_LR.surf.gii /home/data/subjects/MAV006/fs_LR/MNI/MAV006.L.sphere.164k_fs_LR.surf.gii -right-spheres /home/data/scripts/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.R.sphere.32k_fs_LR.surf.gii /home/data/subjects/MAV006/fs_LR/MNI/MAV006.R.sphere.164k_fs_LR.surf.gii']);
parcels_164 = ft_read_cifti_mod('Temp_164.dtseries.nii');
delete('Temp_164.dtseries.nii')

parcel_connectorhubs_connections_map = zeros(length(parcels),length(moduleIDs));
for parcel = 1:length(parcelIDs)
    if summed_PCs_percentile(parcel) >= PC_cutoffs(1)
        thisparcel_connections = unique(parcel_connections(parcel,:));
        parcel_connectorhubs_connections_map(parcels==parcelIDs(parcel),1:length(thisparcel_connections)) = repmat(thisparcel_connections(:)',[nnz(parcels==parcelIDs(parcel)),1]);
    end
end
make_striped_cifti(parcel_connectorhubs_connections_map,0,'temp_thresh',1/25);
temp = ft_read_cifti_mod('temp_thresh.dtseries.nii');
delete('temp_thresh.dtseries.nii')
stripedconnectors = temp;
stripedconnectors.dimord = 'scalar_pos';
stripedconnectors.data = zeros(size(temp.data,1),length(PC_cutoffs));
stripedconnectors.data(:,1) = temp.data;

for i = 2:(length(PC_cutoffs))
    for parcel = 1:length(parcelIDs)
        if summed_PCs_percentile(parcel) >= PC_cutoffs(i)
            stripedconnectors.data(parcels_164.data==parcelIDs(parcel),i) = temp.data(parcels_164.data==parcelIDs(parcel),1);
        end
    end
end
stripedconnectors.dimord = 'scalar_pos';   
    

stripedconnectors.mapname = cifti_template.mapname;
ft_write_cifti_mod('parcel_ConnectorHubs_connections',stripedconnectors)





    