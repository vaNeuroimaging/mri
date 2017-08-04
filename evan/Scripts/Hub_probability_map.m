subjectlist = '/home/data/subjects/DART.txt';

Hubdefinitioncolumn = 8;


%Location of a file describing the surface area of each cortical vertex
surfaceareasfile = '/home/data/scripts/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.LR.midthickness.32k_fs_LR_surfaceareas.dtseries.nii';


%Maximum distance metric (mean closest point-to-point distance) beyond
%which two subjects' patches are not allowed to be matched.
distthresh = 10;

homogeneity_thresh = .7;

common_connections_thresh = 2;

%Minimum subjects who must have a cluster of patches for it to survive in
%the final thresholded output.
pctsubthresh = 1/2;

%Name of file containing data from subjects' discrete contiguous system
%pieces. Will be created or loaded.
system_clustering_names = 'Hub_clustering.mat';

%Labels of systems. The index of each name in this variable should match
%the numeric value of the system in the subjects' maps.
networklabels = {'Default','Visual','FrontoPar','FrontoPar 2','DorsalAttn','DorsalAttn2','VentAttn','Salience','CingOperc','MotorHand','MotorMouth','Auditory','MTL1','MTL2','MedPar','ParOccip'}; %Standard Power colors
%networklabels = {'Default','Visual','FrontoPar','FrontoPar 2','DorsalAttn','DorsalAttn2','VentAttn','Salience','CingOperc','MortorHand','MotorMouth','Auditory','MTL1','STS','MedPar','ParOccip','Orbitofrontal','LatVisual'};

%A file describing vertex-to-vertex geodesic distances. Use a version
%with large cross-hemispheric distances (this prevents an unlikely scenario
%of matching clusters across the medial wall). Using a uint8 version is
%recommended as it will be faster to load, and the script reduces it to
%uint8 anyway to save time.
distances_file = '/home/data/scripts/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.LR.32k_fs_LR.distances_surfgeo_voleuc_normalwall_standardsubcort_uint8.mat';

%Folder where results will be written to.
outputfolder = ['/home/data/Analysis/Hub_clustering/Col' num2str(Hubdefinitioncolumn) '_' num2str(common_connections_thresh) 'commonconnections_' num2str(distthresh) 'mm_homogenous_PC/'];
warning off;
mkdir(outputfolder)
disp(['Writing to ' outputfolder])

ncortverts = 59412;

%% Get info about hubs


neighbors = cifti_neighbors(surfaceareasfile);

hubcount = 0;
subjects = textread(subjectlist,'%s');

surfacearea = ft_read_cifti_mod(surfaceareasfile); 
out_template = surfacearea;
surfacearea = surfacearea.data;

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
        
        clusters_SAs(hubcount,1) = sum(surfacearea(sub_parcels==ID));
        
        subnetworks_bycluster(sub_parcels==ID,s) = hubcount;
        
    end
    
end


    
    %Save all information about all subjects' discontiguous system pieces
    save([outputfolder '/' system_clustering_names],'clusters_subs','clusters_IDs','clusters_SAs','subnetworks_bycluster','-v7.3')

    %% Conduct subject to subject matching

%Load point-to-point geodesic distances
distances = smartload(distances_file);
distances = uint8(distances(1:size(subnetworks_bycluster,1),1:size(subnetworks_bycluster,1)));
distances(1:29696,29697:end) = 1000; distances(29697:end,1:29696) = 1000;

%Set up variables
indices = zeros(0,2);

prevstring = [];
for s1 = 1:length(subjects)
    s1clus_dmat = [];
    for s2 = (s1+1):length(subjects)
        
        string = ['Matching subjects pairwise: Subject ' num2str(s1) ' vs Subject ' num2str(s2)];
        fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
        prevstring = string;
        
        %Run script which matches all clusters in subject 1 with all
        %clusters in subject 2, for each system separately. This spits
        %out columns of matched patches, such that each column in
        %s1matches corresponds with the same numbered column in
        %s2matches. It also spits out patch-to-point distances from
        %subject 1's patches, which will be useful in later loop
        %iterations comparing subject 1 to subjects 3,4,5,etc.
        
        %[matches,s1clus_dmat] = Hubs_match_preclustered(subnetworks_bycluster(:,s1),clusters_IDs(clusters_subs==s1),distances,s1clus_dmat,subnetworks_bycluster(:,s2),clusters_IDs(clusters_subs==s2),distthresh,common_connections_thresh);
        %[matches,s1clus_dmat] = Hubs_match_preclustered2(subnetworks_bycluster(:,s1),clusters_IDs(clusters_subs==s1),distances,s1clus_dmat,subnetworks_bycluster(:,s2),clusters_IDs(clusters_subs==s2),distthresh);
        [matches,s1clus_dmat] = Hubs_match_preclustered3(subnetworks_bycluster(:,s1),clusters_IDs(clusters_subs==s1),distances,s1clus_dmat,subnetworks_bycluster(:,s2),clusters_IDs(clusters_subs==s2),distthresh,common_connections_thresh);
        indices = [indices ; matches];
        
    end
    
end

disp(' ')

%Make a huge [subject patches X subject matches] matrix
match_matrix = false(hubcount);
match_matrix(sub2ind(size(match_matrix),indices(:,1),indices(:,2))) = true;


%Save the huge match matrix
save([outputfolder '/match_matrix.mat'],'match_matrix','-v7.3')

%% Infomap to group matched patches together

prevdir = pwd;
cd(outputfolder)

%Run infomap on the match matrix
mat2pajek_byindex(match_matrix,find(triu(match_matrix,1)),'pajek.net')
rawclrs = run_infomap_on_pajekfile('pajek.net',100);
dlmwrite(['rawassn.txt'],rawclrs,'\t')
delete('pajek.clu')

%Eliminate groups of patches that have only one patch
simple_assigns = modify_clrfile('simplify','rawassn.txt',2);

cd(prevdir)









%% Make outputs

load([outputfolder '/' system_clustering_names])
simple_assigns = load([outputfolder '/rawassn_minsize2.txt']);

%Get a list of all groups of matched patches ID'd by infomap
communities = unique(simple_assigns); communities(communities<1) = [];

%Set up output variables
outdata = zeros(size(subnetworks_bycluster,1),length(communities));
median_size_outlines = zeros(size(subnetworks_bycluster,1),length(communities));

community_IDs = zeros(length(communities),length(networklabels));


%Set up a variable keeping track of how many subjects have a patch in each
%group of patches
numsubs_all = zeros(length(communities),1);

%Set up a variable keeping track of relevant information for each group of
%patches; this will be written into the dscalar output file
columnnames = cell(1,length(communities));

robustIDs = cell(length(communities),1);

%Get the connection distribution for each group of hubs
for communitynum = 1:length(communities)
    community = communities(communitynum);
    clusternums_in_community = find(simple_assigns==community);
    
    for clusternum = clusternums_in_community(:)'
        
        vector_toadd = zeros(1,length(networklabels));
        IDindices = clusters_IDs{clusternum};
        IDindices(logical(mod(IDindices,1))) = []; IDindices(IDindices > length(networklabels)) = [];
        vector_toadd(IDindices) = 1;
        community_IDs(communitynum,:) = community_IDs(communitynum,:) + vector_toadd;
        
    end
    
    communityIDs_string = [];
    IDspresent = find(community_IDs(communitynum,:));
    for thisID = IDspresent(:)'
        
        if (community_IDs(communitynum,thisID) / length(clusternums_in_community)) >= .25
        
            IDstring = [num2str(round(community_IDs(communitynum,thisID) / length(clusternums_in_community) * 100)) '% ' networklabels{thisID} '; '];
            communityIDs_string = [communityIDs_string IDstring];
            
            robustIDs{communitynum}(end+1) = thisID;
            
        end
    end
    communityIDs_string = communityIDs_string(1:end-2);
        
    
    %Get a list of subjects who have patches in this group
    these_subjects = unique(clusters_subs(clusternums_in_community));
    numsubs = length(these_subjects);
    numsubs_all(communitynum) = numsubs;
    
%     %Sometimes a subject has a really big patch that overlaps both this
%     %group of patches and a different group (it's like two patches stuck
%     %together). In this case it will have matched with many patches in the
%     %other group and many in this group, but it only gets put into one of
%     %the two groups. If we're counting how many subjects have a patch in
%     %each group, this should probably count in both groups. This code
%     %checks for that, and add the subject to a separate count of "total"
%     %subjects represented in each group.
%     clusternums_not_in_community = find(simple_assigns~=community);
%     additionalsubs = [];
%     for clusternum = clusternums_not_in_community(:)'
%         if ~any(these_subjects==clusters_subs_IDs_SAs(clusternum,1))
%             matches_to_thiscommunity = clusternums_in_community(logical(match_matrix(clusternum,clusternums_in_community)));
%             if length(unique(clusters_subs_IDs_SAs(matches_to_thiscommunity,1))) > (numsubs/2)
%                 additionalsubs = [additionalsubs clusters_subs_IDs_SAs(clusternum,1)];
%             end
%         end
%     end
%     numsubs_all(communitynum) = numsubs + length(unique(additionalsubs));
%     
%     %Calculate the average dice metrics across patches in the group
%     avg_maxdice = mean(clusters_subs_IDs_SAs(clusternums_in_community,3));
%     avg_maxdicestr = sprintf('%0.3f',avg_maxdice); avg_maxdicestr = avg_maxdicestr(2:end);
%     avg_dicediff = mean(clusters_subs_IDs_SAs(clusternums_in_community,4));
%     avg_dicediffstr = sprintf('%0.3f',avg_dicediff); avg_dicediffstr = avg_dicediffstr(2:end);
    
    %For each subject, calculate total size of all patches they have in the group
    combined_cluster_sizes = zeros(numsubs,1);
    for subnum = 1:numsubs
        thissub = these_subjects(subnum);
        thissub_clusternums = find(clusters_subs==thissub);
        sub_sizes = clusters_SAs(intersect(thissub_clusternums,clusternums_in_community),1);
        combined_cluster_sizes(subnum) = sum(sub_sizes);
    end
    
    %Calculate median patch size across subjects
    median_size = median(combined_cluster_sizes);
    
    %Save information about each patch group, to be put into the output file
    columnnames{1,communitynum} = [num2str(numsubs) ' of ' num2str(length(subjects)) ' subs; median ' num2str(round(median_size)) 'mm2; ' communityIDs_string];
    
    %Add patches to probability map for this group
    for this_clusternum = clusternums_in_community(:)'
        outdata(:,communitynum) = outdata(:,communitynum) + ((subnetworks_bycluster(:,clusters_subs(this_clusternum))==this_clusternum) ./ numsubs);
    end
    
    
    
    %This section of code builds an outline around an object the same size
    %as the median patch size of this group, following the contours of
    %the probability map. This is the "central tendency" object.
    
    %Define this probability map
    thismap = outdata(:,communitynum);
    
    %Get the max value in the probability map
    [maxval, maxi] = max(thismap);
    thissizemap = false(size(thismap));
    
    %Seed the max value
    thissizemap(maxi) = 1;
    
    %Get all unique values in the probability map
    sortedvals = sort(unique(thismap),'descend');
    
    %For each ascending value, grow from the seed if vertices next to the growing region 
    %are below that map value
    for i = 1:length(sortedvals)
        mask = (thismap > sortedvals(i)) & (~thissizemap);
        newverts = 1;
        while ~isempty(newverts)
            maskpos = find(mask);
            thissizemap_neighs = unique(neighbors(logical(thissizemap),2:7));
            thissizemap_neighs(isnan(thissizemap_neighs)) = [];
        
            newverts = intersect(thissizemap_neighs,maskpos);
        
            thissizemap(newverts) = 1;
        
            mask = mask & (~thissizemap);
            
        end
        
        %Stop if the region gets to be the median size
        if sum(surfacearea(logical(thissizemap))) >= median_size;
           break
        end
    end
    
    %Get the outline of that grown region
    thissizemap_outline = zeros(size(thissizemap));
    inds_in_sizemap = find(thissizemap);
    for ind = inds_in_sizemap(:)'
        indneigh = neighbors(ind,2:end); 
        indneigh_withnan = indneigh;
        indneigh(isnan(indneigh)) = [];
        if any(thissizemap(indneigh)==0)
            thissizemap_outline(indneigh(thissizemap(indneigh)==0)) = 1;
        end
        if any(isnan(indneigh_withnan))
            thissizemap_outline(ind) = 1;
        end
    end
    
    %The outline can have funny extrusions sometimes. Get only the biggest
    %contiguous piece of outline
    clustereddata = cifti_cluster_surfacearea_volume(thissizemap_outline,.5,1.5,0,0,ncortverts,surfacearea,neighbors);
    
    %Save that outline to the output file
    median_size_outlines(:,communitynum) = (clustereddata==1);
        
    
end


%% Save out results

%Save probability maps
out_template.data = zeros(size(out_template.data,1),size(outdata,2));
out_template.data(1:ncortverts,:) = outdata(1:ncortverts,:);
out_template.dimord = 'scalar_pos';
out_template.mapname = columnnames;
ft_write_cifti_mod([outputfolder '/Hub_probability_maps_' num2str(distthresh) 'mm'],out_template)

%Save outlines
out_template.data(1:ncortverts,:) = median_size_outlines(1:ncortverts,:);
ft_write_cifti_mod([outputfolder '/Hub_probability_maps_' num2str(distthresh) 'mm_mediansizeobjects'],out_template)



%Get the probability maps where more than the defined pct of subjects have
%clusters in the map
nsubthresh = ceil(length(subjects) * pctsubthresh);
communityindex = find(numsubs_all >= nsubthresh);

%Save out those thresholded probability maps
out_template.data = zeros(size(out_template.data,1),nnz(communityindex));
out_template.data(1:ncortverts,:) = outdata(1:ncortverts,communityindex);
out_template.mapname = columnnames(communityindex);
ft_write_cifti_mod([outputfolder '/Hub_probability_maps_' num2str(distthresh) 'mm_' num2str(nsubthresh) 'sub'],out_template)

%And outlines
out_template.data(1:ncortverts,:) = median_size_outlines(1:ncortverts,communityindex);
ft_write_cifti_mod([outputfolder '/Hub_probability_maps_' num2str(distthresh) 'mm_mediansizeobjects_' num2str(nsubthresh) 'sub'],out_template)



outcollapsed = zeros(size(outdata,1),length(communityindex));
for communitynum = 1:length(communityindex)
    
    for IDnum = 1:length(robustIDs{communityindex(communitynum)})
    
        outcollapsed(logical(median_size_outlines(:,communityindex(communitynum))),IDnum) = robustIDs{communityindex(communitynum)}(IDnum);
        
    end
end

singlecol_out_template = out_template;
singlecol_out_template.data = zeros(size(out_template.data,1),size(outcollapsed,2));
singlecol_out_template.data(1:ncortverts,:) = outcollapsed(1:ncortverts,:);
singlecol_out_template.dimord = 'pos_time';
singlecol_out_template = rmfield(singlecol_out_template,'mapname');
ft_write_cifti_mod([outputfolder '/Hubs_' num2str(distthresh) 'mm_connections_' num2str(nsubthresh) 'sub'],singlecol_out_template)
make_striped_cifti([outputfolder '/Hubs_' num2str(distthresh) 'mm_connections_' num2str(nsubthresh) 'sub.dtseries.nii'],0,[outputfolder '/Hubs_' num2str(distthresh) 'mm_' num2str(nsubthresh) 'sub'])
