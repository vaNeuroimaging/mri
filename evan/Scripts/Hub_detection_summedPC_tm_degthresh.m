function Hub_detection_summedPC_tm_degthresh(parcelsfile,corrmatfile,dmatfile,consensus_modulesfile,kdenthresholds,xdist)
%Hub_detection_summedPC(parcelsfile,corrmatfile,dmatfile,consensus_modulesfile,kdenthresh,xdist)

PC_cutoffs = [.5 : .05 : .95 .98];
pct_kdenthreshes_tocount_for_connection = .75;
pctconnections_needed_tocount_foracommunity = .1;
degreethreshold = .25;

parcels = ft_read_cifti_mod(parcelsfile);
cifti_template = parcels;
parcels = parcels.data;
parcelIDs = unique(parcels); parcelIDs(parcelIDs<1) = [];

corrmat = smartload(corrmatfile);
dmat = smartload(dmatfile);
%communities = load(communitiesfile);

corrmat(dmat<xdist) = 0;
triumask = triu(true(length(corrmat)));
sorted = sort(corrmat(triumask),'descend');

modules = ft_read_cifti_mod(consensus_modulesfile); modules = modules.data;
parcelmodules = zeros(size(corrmat,1),1);
for parcelIDnum = 1:length(parcelIDs)
    parcelmodules(parcelIDnum) = mode(modules(parcels==parcelIDs(parcelIDnum)));
end
moduleIDs = unique(parcelmodules); moduleIDs(moduleIDs<1) = [];

corrmat_bin_all = zeros(size(corrmat,1),size(corrmat,1),length(kdenthresholds));



parcel_degree_pctiles = zeros(length(parcelIDs),length(kdenthresholds));
for threshnum = 1:length(kdenthresholds);
    sortedvec = sort(corrmat(triu(true(size(corrmat)),1)),'descend');
    rthresh = sortedvec(floor(length(sortedvec) * kdenthresholds(threshnum)));
     bin_corrmat = corrmat > rthresh; %bin_corrmat(diag(true(size(corrmat,1)),0)) = 0;
     corrmat_bin_all(:,:,threshnum) = bin_corrmat;
%     degree = sum(bin_corrmat,2);
    wei_corrmat = corrmat .* (corrmat > rthresh); 
    %wei_corrmat(diag(true(size(corrmat,1)),0)) = 0;
    degree = sum(wei_corrmat,2);
    [~,sorti] = sort(degree);
    parcel_degree_pctiles(sorti,threshnum) = [1:length(degree)] / length(degree);
end
mean_degree_pctile = mean(parcel_degree_pctiles,2);


parcel_correlation_vec = corrmat(triu(true(size(corrmat)),1));
sorted_vals = sort(parcel_correlation_vec,'descend');
all_PCs = zeros(length(parcelIDs),length(kdenthresholds));
for threshnum = 1:length(kdenthresholds)
    kdenthresh = kdenthresholds(threshnum);
    rthresh = sorted_vals(ceil(length(sorted_vals) .* kdenthresh));
    corrmat_weighted = corrmat .* (corrmat > rthresh);
    
    parcel_degree = sum(corrmat_weighted,2);
    moduleconnectionssum = zeros(size(corrmat_weighted,1),1);
    
    for communitynum = 1:length(moduleIDs)
        communityID = moduleIDs(communitynum);
        communityindices = parcelmodules==communityID;
        parcel_community_degree = sum(corrmat_weighted(:,communityindices),2);
        parcel_community_ratio = (parcel_community_degree ./ parcel_degree);
        moduleconnectionssum = moduleconnectionssum + (parcel_community_ratio.^2);
    end
    
    all_PCs(:,threshnum) = 1-moduleconnectionssum;
end

% all_PC_pcts = zeros(size(all_PCs));
% for i = 1:size(all_PCs,2);
%     nonanvec = all_PCs(:,i); nonanvec(isnan(nonanvec)) = 0;
%     [~,sorti] = sort(nonanvec);
%     all_PC_pcts(sorti,i) = [1:length(nonanvec)] ./ length(nonanvec);
% end
% all_PC_pcts(isnan(all_PCs)) = NaN;
%
% mean_PCs = nanmean(all_PC_pcts,2);
mean_PCs = nanmean(all_PCs,2);

mean_PCs(isnan(mean_PCs)) = 0;



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



save('summed_PCs.mat','mean_PCs')
save('parcel_connections.mat','parcel_connections')
save('mean_degree_pctile.mat','mean_degree_pctile');

mean_PCs_orig = mean_PCs;
mean_PCs(mean_degree_pctile<degreethreshold) = 0;

[~,sortedindex] = sort(mean_PCs);
ordered(sortedindex) = 1:length(mean_PCs);
summed_PCs_percentile = ordered ./ length(mean_PCs);

[~,sortedindex] = sort(mean_PCs_orig);
ordered(sortedindex) = 1:length(mean_PCs_orig);
summed_PCs_percentile_orig = ordered ./ length(mean_PCs_orig);

summed_PCs_map = zeros(size(parcels));
summed_PCs_orig_map = zeros(size(parcels));
degree_percentile_map = zeros(size(parcels));
parcel_connections_map = zeros(size(parcels,1),size(parcel_connections,2));
for parcelnum = 1:length(parcelIDs)
    summed_PCs_map(parcels==parcelIDs(parcelnum)) = summed_PCs_percentile(parcelnum);
    summed_PCs_orig_map(parcels==parcelIDs(parcelnum)) = summed_PCs_percentile_orig(parcelnum);
    degree_percentile_map(parcels==parcelIDs(parcelnum)) = mean_degree_pctile(parcelnum);
    parcel_connections_map(parcels==parcelIDs(parcelnum),:) = repmat(parcel_connections(parcelnum,:),nnz(parcels==parcelIDs(parcelnum)),1);
end

cifti_template.data = summed_PCs_map;
ft_write_cifti_mod('parcel_PCs_percentile',cifti_template)  

cifti_template.data = summed_PCs_orig_map;
ft_write_cifti_mod('parcel_PCs_percentile_orig',cifti_template)  

cifti_template.data = parcel_connections_map;
ft_write_cifti_mod('parcel_connections',cifti_template)  

cifti_template.data = degree_percentile_map;
ft_write_cifti_mod('degree_percentile',cifti_template)  


parcel_connectorhubs_map = zeros(length(parcels),length(PC_cutoffs));

for i = 1:length(PC_cutoffs)
    
    parcel_connectorhubs_map(:,i) = (summed_PCs_map >= PC_cutoffs(i))  .* (parcels > 0);
    cifti_template.mapname{i} = ['PC percentile thresh=' num2str(PC_cutoffs(i))];
    
end

cifti_template.data = parcel_connectorhubs_map;
cifti_template.dimord = 'scalar_pos';
ft_write_cifti_mod('parcel_ConnectorHubs',cifti_template)


% evalc(['!wb_command -cifti-resample ' parcelsfile ' COLUMN /data/cn/data1/scripts/CIFTI_RELATED/Resources/Conte69_atlas.LR.164k_fs_LR/Conte69.LR.164k_fs_LR.normalwall_surfaceonly_template.dtseries.nii COLUMN BARYCENTRIC ENCLOSING_VOXEL Temp_164.dtseries.nii -surface-largest -left-spheres /data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k/Conte69.L.sphere.32k_fs_LR.surf.gii /data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/vc25125/7112b_fs_LR/vc25125.L.sphere.164k_fs_LR.surf.gii -right-spheres /data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k/Conte69.R.sphere.32k_fs_LR.surf.gii /data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/vc25125/7112b_fs_LR/vc25125.R.sphere.164k_fs_LR.surf.gii']);
% parcels_164 = ft_read_cifti_mod('Temp_164.dtseries.nii');
% delete('Temp_164.dtseries.nii')
% 
% parcel_connectorhubs_connections_map = zeros(length(parcels),length(moduleIDs));
% for parcel = 1:length(parcelIDs)
%     if summed_PCs_percentile(parcel) >= PC_cutoffs(1)
%         thisparcel_connections = unique(parcel_connections(parcel,:));
%         parcel_connectorhubs_connections_map(parcels==parcelIDs(parcel),1:length(thisparcel_connections)) = repmat(thisparcel_connections(:)',[nnz(parcels==parcelIDs(parcel)),1]);
%     end
% end
% make_striped_cifti(parcel_connectorhubs_connections_map,0,'temp_thresh',1/25);
% temp = ft_read_cifti_mod('temp_thresh.dtseries.nii');
% delete('temp_thresh.dtseries.nii')
% stripedconnectors = temp;
% stripedconnectors.dimord = 'scalar_pos';
% stripedconnectors.data = zeros(size(temp.data,1),length(PC_cutoffs));
% stripedconnectors.data(:,1) = temp.data;
% 
% for i = 2:(length(PC_cutoffs))
%     for parcel = 1:length(parcelIDs)
%         if summed_PCs_percentile(parcel) >= PC_cutoffs(i)
%             stripedconnectors.data(parcels_164.data==parcelIDs(parcel),i) = temp.data(parcels_164.data==parcelIDs(parcel),1);
%         end
%     end
% end
% stripedconnectors.dimord = 'scalar_pos';   
%     
% 
% stripedconnectors.mapname = cifti_template.mapname;
% ft_write_cifti_mod('parcel_ConnectorHubs_connections',stripedconnectors)





    