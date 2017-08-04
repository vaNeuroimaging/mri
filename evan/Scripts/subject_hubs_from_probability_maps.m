function subject_hubs_from_probability_maps(probmapfile,mapnums,directory)
%subject_patches_from_probability_maps(probmapfile,mapnums,patch_probability_map_params)

if exist('directory')
cd(directory)
end

subjectlist = '/home/data/subjects/DART.txt';

cifti_out_164_templatefile = '/home/data/scripts/Resources/Conte69_atlas.LR.164k_fs_LR/Conte69.LR.164k_fs_LR.normalwall_surfaceonly_template.dtseries.nii';

homogeneity_thresh = .7;

Hubdefinitioncolumn = 3;


hubcount = 0;
subjects = textread(subjectlist,'%s');

%subnetworks_bycluster = zeros(size(subnetworks));

for s = 1:length(subjects)
    sub_parcels = ft_read_cifti_mod(['/home/data/subjects/' subjects{s} '/parcellation/RSFC_parcels_edgethresh_0.5.dtseries.nii']); sub_parcels = sub_parcels.data;
    sub_hubs = ft_read_cifti_mod(['/home/data/subjects/' subjects{s} '/parcellation/parcel_infomap/parcel_ConnectorHubs.dscalar.nii']); sub_hubs = sub_hubs.data;
    sub_connections = ft_read_cifti_mod(['/home/data/subjects/' subjects{s} '/parcellation/parcel_infomap/parcel_connections.dtseries.nii']); sub_connections = sub_connections.data;
    sub_homogeneity = ft_read_cifti_mod(['/home/data/subjects/' subjects{s} '/parcellation/homogeneity_testing/PCA_eigval_per_first_RSFC_parcels_edgethresh_0.5.dtseries.nii']); sub_homogeneity = sub_homogeneity.data;
    
    if s==1
        subnetworks_bycluster = zeros(size(sub_parcels,1),length(subjects));
    end
    
    sub_hubs = sub_hubs .* repmat((sub_homogeneity > homogeneity_thresh),1,size(sub_hubs,2));
    
    connector_IDs = unique(sub_parcels .* ((sub_hubs(:,Hubdefinitioncolumn)==3) | (sub_hubs(:,Hubdefinitioncolumn)==1))); 
    connector_IDs(connector_IDs<1) = [];
    
    for ID = connector_IDs(:)'
        hubcount = hubcount + 1;
        
        cluster_connections = unique(sub_connections(sub_parcels==ID,:)); cluster_connections(cluster_connections<1) = [];
        
        clusters_subs(hubcount,1) = s;
        
        clusters_IDs{hubcount,1} = cluster_connections;
        
%        clusters_SAs(hubcount,1) = sum(surfacearea(sub_parcels==ID));
        
        subnetworks_bycluster(sub_parcels==ID,s) = hubcount;
        
    end
    
end





probmaps = ft_read_cifti_mod(probmapfile);
cifti_template = probmaps; %cifti_template.data = [];

simple_assigns = load('rawassn_minsize2.txt');

communities = unique(simple_assigns); communities(communities<1) = [];


for mapnum = mapnums

    rs = zeros(length(communities),1);

    probmap = probmaps.data(1:59412,mapnum);
    
    
    
    for communitynum = 1:length(communities)
        community = communities(communitynum);
        
        communitymap = zeros(size(subnetworks_bycluster,1),1);
        
        patchnums_incommunity = find(simple_assigns==community);
        
        for p = 1:length(patchnums_incommunity)
            patchnum = patchnums_incommunity(p);
            sub = clusters_subs(patchnum,1); 
            communitymap = communitymap + (subnetworks_bycluster(:,sub)==patchnum);
        end
        
        rs(communitynum) = paircorr_mod(communitymap,probmap);
        
        if rs(communitynum)>.999
            submaps = false(size(subnetworks_bycluster));
            for p = 1:length(patchnums_incommunity)
                patchnum = patchnums_incommunity(p);
                sub = clusters_subs(patchnum,1);
                submaps(:,sub) = submaps(:,sub) | (subnetworks_bycluster(:,sub)==patchnum);
            end
%             tokens = tokenize(probmaps.mapname{mapnum},':');
%             columnnetworkname = tokens{1};
%             color = find(strcmp(columnnetworkname,networklabels));
%             submaps = submaps .* color;
            break
        else
            clear submaps
        end
    end
    
    if ~exist('submaps')
        error('No matching probability map')
    end
    
    cifti_template.mapname = [];
    for s = 1:size(submaps,1);
        cifti_template.mapname{s} = ['Subject ' num2str(s)];
    end
    cifti_template.data = zeros(size(cifti_template.data,1),size(submaps,2)); cifti_template.data(1:59412,:) = submaps;
    
    clear submaps
    
    ft_write_cifti_mod('TEMP_Submaps',cifti_template);
    
    system(['wb_command -cifti-resample TEMP_Submaps.dscalar.nii COLUMN ' cifti_out_164_templatefile ' COLUMN BARYCENTRIC ENCLOSING_VOXEL Temp_164.dscalar.nii -surface-largest -left-spheres /home/data/scripts/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.L.sphere.32k_fs_LR.surf.gii /home/data/subjects/MAV006/fs_LR/MNI/MAV006.L.sphere.164k_fs_LR.surf.gii -right-spheres /home/data/scripts/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.R.sphere.32k_fs_LR.surf.gii /home/data/subjects/MAV006/fs_LR/MNI/MAV006.R.sphere.164k_fs_LR.surf.gii']);
    
    masks = ft_read_cifti_mod('Temp_164.dscalar.nii');
    
    out = masks;
    out.data = zeros(size(out.data));
    
    for s = 1:length(subjects)
        
        striped = ft_read_cifti_mod(['/home/data/subjects/' subjects{s} '/parcellation/parcel_infomap/parcel_ConnectorHubs_connections.dscalar.nii']);
        out.data(:,s) = masks.data(:,s) .* striped.data(:,1);
        
    end
    
    ft_write_cifti_mod(['Submaps_probmap' num2str(mapnum) '_' probmapfile],out);
    
end



        
        
        
    