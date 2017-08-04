function calc_PC_wmDeg(parcelsfile,corrmatfile,dmatfile,modulesfile)

cutoff_percentiles = [.5 : .05 : .9]; 


parcels = ft_read_cifti_mod(parcelsfile);
cifti_template = parcels;
parcels = parcels.data;

corrmat = smartload(corrmatfile);
dmat = smartload(dmatfile);

modules = ft_read_cifti_mod(modulesfile);
thisthresh_info = modules.mapname{1};
modules = modules.data;
modules(isnan(modules)) = 0;
moduleIDs = unique(modules); moduleIDs(moduleIDs<=0) = [];

infotokens = tokenize(thisthresh_info,',');
kdentokens = tokenize(infotokens{1},'=');
kdenthresh = str2num(kdentokens{2});
xdisttokens = tokenize(infotokens{1},'=');
xdist = str2num(xdisttokens{2});

corrmat(dmat<xdist) = 0;
triumask = triu(true(length(corrmat)));
sorted = sort(corrmat(triumask),'descend');
kdennum = ceil(nnz(triumask .* (dmat>xdist)) .* kdenthresh);
rthresh = sorted(kdennum);
corrmat(corrmat < rthresh) = 0;




parcelIDs = unique(parcels); parcelIDs(parcelIDs<1) = [];

parcelmodules = zeros(length(parcelIDs),1);
for parcel = 1:length(parcelIDs);
    parcelmodules(parcel) = mode(modules(parcels==parcelIDs(parcel)));
end

parcel_PCs = zeros(length(parcelIDs),1);
parcel_PCs_map = zeros(length(parcels),1);
parcel_wmDegs = zeros(length(parcelIDs),1);
parcel_degs_bin = zeros(length(parcelIDs),1);
parcel_wmDegsZ = zeros(length(parcelIDs),1);
parcel_wmDegsZ_map = zeros(length(parcels),1);

for parcel = 1:length(parcelIDs)
    thismodule = parcelmodules(parcel);
    parcel_wmDegs(parcel) = sum(corrmat(parcel,parcelmodules==thismodule));
    parcel_degs_bin(parcel) = sum(logical(corrmat(parcel,:)));
end

for parcel = 1:length(parcelIDs)
    thismodule = parcelmodules(parcel);
    thismodule_Degs = parcel_wmDegs(parcelmodules==thismodule);
    parcel_wmDegsZ(parcel) = (parcel_wmDegs(parcel) - mean(thismodule_Degs)) ./ std(thismodule_Degs);
    parcel_wmDegsZ_map(parcels==parcelIDs(parcel)) = parcel_wmDegsZ(parcel);
    
    moduleconnectionssum = 0;
    for m = 1:length(moduleIDs)
        linkstothismodule = nnz(corrmat(parcel,parcelmodules==moduleIDs(m)));
        moduleconnectionssum = moduleconnectionssum + ((linkstothismodule / parcel_degs_bin(parcel)).^2);
    end
    parcel_PCs(parcel) = 1-moduleconnectionssum;
    parcel_PCs_map(parcels==parcelIDs(parcel)) = parcel_PCs(parcel);
end

cifti_template.data = parcel_wmDegsZ_map;
ft_write_cifti_mod('parcel_wmDegZs',cifti_template)

cifti_template.data = parcel_PCs_map;
ft_write_cifti_mod('parcel_PCs',cifti_template)


parcel_connectorhubs_map = zeros(length(parcels),length(cutoff_percentiles));
PCthresh = zeros(1,length(cutoff_percentiles));
Degthresh = zeros(1,length(cutoff_percentiles));

for i = 1:length(cutoff_percentiles)
    percentile = cutoff_percentiles(i);
    sortedPCs = sort(parcel_PCs,'ascend');
    PCthresh(i) = sortedPCs(ceil(numel(sortedPCs) .* percentile));
    sortedDegs = sort(parcel_wmDegsZ,'ascend');
    Degthresh(i) = sortedDegs(ceil(numel(sortedDegs) .* percentile));
    parcel_connectorhubs_map(:,i) = ((parcel_PCs_map >= PCthresh(i)) + (parcel_wmDegsZ_map >= Degthresh(i)).*2) .* (parcels > 0);
    cifti_template.mapname{i} = ['%tile thresh = ' num2str(percentile) ' for PC and w/in module Degree'];
    
end

cifti_template.data = parcel_connectorhubs_map;
cifti_template.dimord = 'scalar_pos';
ft_write_cifti_mod('parcel_ConnectorHubs',cifti_template)


evalc(['!wb_command -cifti-resample ' parcelsfile ' COLUMN /home/data/scripts/Resources/Conte69_atlas.LR.164k_fs_LR/Conte69.LR.164k_fs_LR.normalwall_surfaceonly_template.dtseries.nii COLUMN BARYCENTRIC ENCLOSING_VOXEL Temp_164.dtseries.nii -surface-largest -left-spheres /home/data/scripts/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.L.sphere.32k_fs_LR.surf.gii /home/data/subjects/MAV006/fs_LR/MNI/MAV006.L.sphere.164k_fs_LR.surf.gii -right-spheres /home/data/scripts/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.R.sphere.32k_fs_LR.surf.gii /home/data/subjects/MAV006/fs_LR/MNI/MAV006.R.sphere.164k_fs_LR.surf.gii']);
parcels_164 = ft_read_cifti_mod('Temp_164.dtseries.nii');
delete('Temp_164.dtseries.nii')

parcel_connectorhubs_connections_map = zeros(length(parcels),length(moduleIDs));
for parcel = 1:length(parcelIDs)
    if (parcel_PCs(parcel) >= PCthresh(1)) && (parcel_wmDegsZ(parcel) >= Degthresh(1))
        thisparcel_connections = unique(parcelmodules(logical(corrmat(parcel,:))));
        parcel_connectorhubs_connections_map(parcels==parcelIDs(parcel),1:length(thisparcel_connections)) = repmat(thisparcel_connections(:)',[nnz(parcels==parcelIDs(parcel)),1]);
    end
end
make_striped_cifti(parcel_connectorhubs_connections_map,0,'temp_thresh');
temp = ft_read_cifti_mod('temp_thresh.dtseries.nii');
delete('temp_thresh.dtseries.nii')
stripedconnectors = temp;
stripedconnectors.dimord = 'scalar_pos';
stripedconnectors.data = zeros(size(temp.data,1),length(cutoff_percentiles));
stripedconnectors.data(:,1) = temp.data;

for i = 2:(length(cutoff_percentiles))
    for parcel = 1:length(parcelIDs)
        if (parcel_PCs(parcel) >= PCthresh(i)) && (parcel_wmDegsZ(parcel) >= Degthresh(i))
            stripedconnectors.data(parcels_164.data==parcelIDs(parcel),i) = temp.data(parcels_164.data==parcelIDs(parcel),1);
        end
    end
end
stripedconnectors.dimord = 'scalar_pos';   
    

stripedconnectors.mapname = cifti_template.mapname;
ft_write_cifti_mod('parcel_ConnectorHubs_connections',stripedconnectors)
    