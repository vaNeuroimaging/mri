%function Patch_probability_map_grouptoinds
%Patch_probability_map(patch_probability_map_params)
%
% Given a previously-calculated set of subject-specific system maps, this
% function breaks these maps into discrete contiguous pieces ("patches")
% and then attempts to match these patches across individuals using a
% distance-based metric (geodesic distance needed to travel from
% each point in one subject's patch to the closest point in the other
% subject's patch). Matched patches are then grouped across subjects
% using Infomap to identify clusters of nearby patches of the same system.
% These clusters are represented as probabilistic maps. 
%
% The script also creates "central tendency"-esque objects that are the
% median size of the grouped subject-level patches and follow the contours
% of the proability map (these are saved as outlines for easy superposition
% on top of the probability maps).
%
% This can take a looong time to run. The time needed increases with the
% number of subjects being matched (because all pairwise subject
% comparisons are made).
%
%
% Requires a parameters file (a .m file) which will be executed to load
% needed parameters, including:
%
% the subject-level systems maps
% the surface area of each cortical vertex
% the distance threshold beyond which patches do not match
% the number of subjects required to have a patch
% whether to define contiguous clusters (or use definitions from a previous
%   run of the script) 
% the name of a file containing contigiuous cluster information, to be
%   saved or loaded
% names for systems
% a file containing all point-to-point geodesic distances. Use a version
%   with large cross-hemispheric distances (this prevents an unlikely scenario
%   of matching clusters across the medial wall). 
% the location where outputs will be written to
%
%
% Requires Infomap scripts to be in your path (e.g.,
% /data/cn/data1/scripts/CIFTI_RELATED/Infomap/ and subfolders)
%
% Requires the Cifti Resources scripts to be in your path (e.g.,
% /data/cn/data1/scripts/CIFTI_RELATED/Resources/ and subfolders)
%
%EMG 09/01/15


%% Load data

% %Find params file
% [paramspath,paramsname,paramsextension] = fileparts(patch_probability_map_params);
% origpath = pwd;
% if ~isempty(paramspath)
%     cd(paramspath)
% end
% 
% %Load parameters
% params = feval(paramsname);
% varnames = fieldnames(params);
% for i = 1:length(varnames)
%     evalc([varnames{i} ' = params.' varnames{i}]);
% end
% clear varnames params
% 
% 
% %Load dice metrics, if they exist
% [path, temp, ~] = fileparts(subnetworksfile);
% [~,filestem,~] = fileparts(temp);
% maxdicesfile = [path '/' filestem '_dicemaxes.dscalar.nii'];
% dicediffsfile = [path '/' filestem '_dicediffs.dscalar.nii'];
% 
% if exist(maxdicesfile) && exist(maxdicesfile)
%     calc_dices = true;
%     maxdices = ft_read_cifti_mod(maxdicesfile); maxdices = maxdices.data;
%     dicediffs = ft_read_cifti_mod(dicediffsfile); dicediffs = dicediffs.data;
% else
%     calc_dices = false;
% end

system_clustering_names = 'System_clustering.mat';
distthresh = 20;

%Load subject systems
subnetworks = ft_read_cifti_mod('Allsubsplusgroup_maps_15verts.dtseries.nii'); 
out_template = subnetworks; out_template.data = [];
subnetworks = subnetworks.data;

%Get the identities of all systems
IDs = unique(subnetworks(:,1)); IDs(IDs<1) = [];

distances_file = '/home/data/scripts/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.LR.32k_fs_LR.distances_surfgeo_voleuc_normalwall_xhemlarge_standardsubcort_uint8.mat';
%Load point-to-point geodesic distances
distances = smartload(distances_file);
distances = distances(1:59412,1:59412);
distances = uint8(distances);

new_cluster = 0;

%Load neighbor adjacency information
%neighbors = cifti_neighbors(subnetworks);

%Load surface areas
%surfacearea = ft_read_cifti_mod(surfaceareasfile);
%surfacearea = surfacearea.data;

%% Find discrete contiguous clusters of systems

if new_cluster
    
    %ID the number of cortical vertices
    ncortverts = nnz(out_template.brainstructure == 1) + nnz(out_template.brainstructure == 2);
    
    %Set up clustering outputs
    clusters_subs_IDs_dices_SAs = zeros(0,5);
    subnetworks_bycluster = zeros(size(subnetworks));
    clustercount = 0;
    
    prevstring = [];
    
    %For each subject
    for s = 1:size(subnetworks,2)
        string = ['Clustering systems: Subject ' num2str(s)];
        fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
        prevstring = string;
        
        %For each system
        for ID = IDs(:)'
            %Find discontiguous pieces of that system
            clustereddata = cifti_cluster(subnetworks(:,s),ID-.05,ID+.05,0);
            clustereddata = sum((clustereddata .* repmat([1:size(clustereddata,2)],size(clustereddata,1),1)) , 2);
            
            %Get IDs of those pieces within the clustered map
            unique_clusters = unique(clustereddata); unique_clusters(unique_clusters==0) = [];
            
            %For each piece
            for clusternum = 1:length(unique_clusters)
                
                %Add to the running count of pieces
                clustercount = clustercount + 1;
                
                %Get the identity fo this piece
                cluster_ID = unique_clusters(clusternum);
                
                %Get dice quality information of this piece
%                 if calc_dices
%                     this_maxdice = mean(maxdices((clustereddata==cluster_ID),s));
%                     this_dicediff = mean(dicediffs((clustereddata==cluster_ID),s));
%                 else
                    this_maxdice = 0;
                    this_dicediff = 0;
                %end
                
                %Get total surface area of this piece
                cluster_surfacearea = nnz(clustereddata==cluster_ID);
                
                %Save information about this piece
                clusters_subs_IDs_dices_SAs(clustercount,:) = [s ID this_maxdice this_dicediff cluster_surfacearea];
                
                %Save the piece itself
                subnetworks_bycluster((clustereddata==cluster_ID),s) = clustercount;
            end
        end
    end
    disp(' ')
    
    %Save all information about all subjects' discontiguous system pieces
    save(system_clustering_names,'clusters_subs_IDs_dices_SAs','subnetworks_bycluster','clustercount','-v7.3')

else
    %Load all information about all subjects' discontiguous system pieces
    load(system_clustering_names)
end

%% Match patches subject-to-subject 


%Set up variables
indices = {[]};
count = 0;

prevstring = [];
for s1 = 1%:size(subnetworks,2)
    s1clus_dmat = [];
    for s2 = (s1+1):size(subnetworks,2)
        
        
        string = ['Matching subjects pairwise: Subject ' num2str(s1) ' vs Subject ' num2str(s2) ', '];
        fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
        prevstring = [string 'Network #' num2str(max(IDs)) '  '];
        
        %Run script which matches all clusters in subject 1 with all
        %clusters in subject 2, for each system separately. This spits
        %out columns of matched patches, such that each column in
        %s1matches corresponds with the same numbered column in
        %s2matches. It also spits out patch-to-point distances from
        %subject 1's patches, which will be useful in later loop
        %iterations comparing subject 1 to subjects 3,4,5,etc.
        [s1matches,s2matches,s1clus_dmat] = networkClusters_match_preclustered(subnetworks_bycluster(:,s1),[find(clusters_subs_IDs_dices_SAs(:,1)==s1) clusters_subs_IDs_dices_SAs(clusters_subs_IDs_dices_SAs(:,1)==s1,2)],distances,s1clus_dmat,subnetworks_bycluster(:,s2),[find(clusters_subs_IDs_dices_SAs(:,1)==s2) clusters_subs_IDs_dices_SAs(clusters_subs_IDs_dices_SAs(:,1)==s2,2)],distthresh);
        
        %For each matched patch
        for matchnum = 1:size(s1matches,2)
            
            %Get the overall indices of the matched patches (within the
            %huge subnetworks_bycluster variable)
            s1_clusternums = unique(subnetworks_bycluster(logical(s1matches(:,matchnum)),s1));
            s2_clusternums = unique(subnetworks_bycluster(logical(s2matches(:,matchnum)),s2));
            
            %Save those indices
            count = count+1;
            indices{count,1} = s1_clusternums;
            indices{count,2} = s2_clusternums;
            
        end
    end
end
disp(' ')
clear distances

%Make a huge [subject patches X subject matches] matrix
match_matrix = false(clustercount);

%For each subject-to-subject match
for i = 1:size(indices,1)
    
    %Make the indices of the huge match matrix true if those patches were
    %matched between the two subjects
    match_matrix(indices{i,1},indices{i,2}) = true;  match_matrix(indices{i,2},indices{i,1}) = true;
end


nfirstsub_patches = nnz(clusters_subs_IDs_dices_SAs(:,1) == 1);
match_matrix = match_matrix(:,1:nfirstsub_patches);

%Save the huge match matrix
save(['match_matrix.mat'],'match_matrix','-v7.3')

for s = 2:size(subnetworks,2)
clus_inds = find(clusters_subs_IDs_dices_SAs(:,1)==s);
pairwise_mat = match_matrix(clus_inds,:);
numgroupclusters_matched_bysubject(s) = nnz(any(pairwise_mat,1));
end

%% Infomap to group matched patches together
% 
% prevdir = pwd;
% cd(outputfolder)
% 
% %Run infomap on the match matrix
% mat2pajek_byindex(match_matrix,find(triu(match_matrix,1)),'pajek.net')
% rawclrs = run_infomap_on_pajekfile('pajek.net',100);
% dlmwrite(['rawassn.txt'],rawclrs,'\t')
% delete('pajek.clu')
% 
% %Eliminate groups of patches that have only one patch
% simple_assigns = modify_clrfile('simplify','rawassn.txt',2);
% 
% cd(prevdir)
% 
% 
% %% Make outputs
% 
% %Get a list of all groups of matched patches ID'd by infomap
% communities = unique(simple_assigns); communities(communities<1) = [];
% 
% %Set up output variables
% outdata = zeros(size(subnetworks,1),length(communities));
% median_size_outlines = zeros(size(subnetworks,1),length(communities));
% 
% %Get the system identity for each group of patches
% for communitynum = 1:length(communities)
%     community = communities(communitynum);
%     clusternums_in_community = find(simple_assigns==community);
%     community_IDs(communitynum) = mode(clusters_subs_IDs_dices_SAs(clusternums_in_community,2));
% end
% 
% %Sort the patches by system (to group e.g. all Default patches together)
% [~,sorti] = sort(community_IDs,'ascend');
% sorted_communities = communities(sorti);
% 
% %Set up a variable keeping track of how many subjects have a patch in each
% %group of patches
% numsubs_all = zeros(length(sorted_communities),1);
% 
% %Set up a variable keeping track of relevant information for each group of
% %patches; this will be written into the dscalar output file
% columnnames = cell(1,length(sorted_communities));
% 
% %For each group of patches
% for communitynum = 1:length(sorted_communities)
%     community = sorted_communities(communitynum);
%     
%     %Get a list of patches in the group
%     clusternums_in_community = find(simple_assigns==community);
%     
%     %Get the system ID of this group
%     this_ID = mode(clusters_subs_IDs_dices_SAs(clusternums_in_community,2));
%     
%     %Get a list of subjects who have patches in this group
%     these_subjects = unique(clusters_subs_IDs_dices_SAs(clusternums_in_community,1));
%     numsubs = length(these_subjects);
%     
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
%         if ~any(these_subjects==clusters_subs_IDs_dices_SAs(clusternum,1))
%             matches_to_thiscommunity = clusternums_in_community(logical(match_matrix(clusternum,clusternums_in_community)));
%             if length(unique(clusters_subs_IDs_dices_SAs(matches_to_thiscommunity,1))) > (numsubs/2)
%                 additionalsubs = [additionalsubs clusters_subs_IDs_dices_SAs(clusternum,1)];
%             end
%         end
%     end
%     numsubs_all(communitynum) = numsubs + length(unique(additionalsubs));
%     
%     %Calculate the average dice metrics across patches in the group
%     avg_maxdice = mean(clusters_subs_IDs_dices_SAs(clusternums_in_community,3));
%     avg_maxdicestr = sprintf('%0.3f',avg_maxdice); avg_maxdicestr = avg_maxdicestr(2:end);
%     avg_dicediff = mean(clusters_subs_IDs_dices_SAs(clusternums_in_community,4));
%     avg_dicediffstr = sprintf('%0.3f',avg_dicediff); avg_dicediffstr = avg_dicediffstr(2:end);
%     
%     %For each subject, calculate total size of all patches they have in the group
%     combined_cluster_sizes = zeros(numsubs,1);
%     for subnum = 1:numsubs
%         thissub = these_subjects(subnum);
%         thissub_clusternums = find(clusters_subs_IDs_dices_SAs(:,1)==thissub);
%         sub_sizes = clusters_subs_IDs_dices_SAs(intersect(thissub_clusternums,clusternums_in_community),5);
%         combined_cluster_sizes(subnum) = sum(sub_sizes);
%     end
%     
%     %Calculate median patch size across subjects
%     median_size = median(combined_cluster_sizes);
%     
%     %Save information about each patch group, to be put into the output file
%     columnnames{1,communitynum} = [networklabels{this_ID} ': ' num2str(numsubs) ' of ' num2str(size(subnetworks,2)) ' subs; ' num2str(numsubs_all(communitynum)) ' subs total; median size=' num2str(round(median_size)) 'mm2; Dice=' avg_maxdicestr ', Dice diff=' avg_dicediffstr];
%     
%     %Add patches to probability map for this group
%     for this_clusternum = clusternums_in_community(:)'
%         outdata(:,communitynum) = outdata(:,communitynum) + ((subnetworks_bycluster(:,clusters_subs_IDs_dices_SAs(this_clusternum,1))==this_clusternum) ./ numsubs);
%     end
%     
%     
%     
%     %This section of code builds an outline around an object the same size
%     %as the median patch size of this group, following the contours of
%     %the probability map. This is the "central tendency" object.
%     
%     %Define this probability map
%     thismap = outdata(:,communitynum);
%     
%     %Get the max value in the probability map
%     [maxval, maxi] = max(thismap);
%     thissizemap = false(size(thismap));
%     
%     %Seed the max value
%     thissizemap(maxi) = 1;
%     
%     %Get all unique values in the probability map
%     sortedvals = sort(unique(thismap),'descend');
%     
%     %For each ascending value, grow from the seed if vertices next to the growing region 
%     %are below that map value
%     for i = 1:length(sortedvals)
%         mask = (thismap > sortedvals(i)) & (~thissizemap);
%         newverts = 1;
%         while ~isempty(newverts)
%             maskpos = find(mask);
%             thissizemap_neighs = unique(neighbors(logical(thissizemap),2:7));
%             thissizemap_neighs(isnan(thissizemap_neighs)) = [];
%         
%             newverts = intersect(thissizemap_neighs,maskpos);
%         
%             thissizemap(newverts) = 1;
%         
%             mask = mask & (~thissizemap);
%             
%         end
%         
%         %Stop if the region gets to be the median size
%         if sum(surfacearea(logical(thissizemap))) >= median_size;
%            break
%         end
%     end
%     
%     %Get the outline of that grown region
%     thissizemap_outline = zeros(size(thissizemap));
%     inds_in_sizemap = find(thissizemap);
%     for ind = inds_in_sizemap(:)'
%         indneigh = neighbors(ind,2:end); 
%         indneigh_withnan = indneigh;
%         indneigh(isnan(indneigh)) = [];
%         if any(thissizemap(indneigh)==0)
%             thissizemap_outline(indneigh(thissizemap(indneigh)==0)) = 1;
%         end
%         if any(isnan(indneigh_withnan))
%             thissizemap_outline(ind) = 1;
%         end
%     end
%     
%     %The outline can have funny extrusions sometimes. Get only the biggest
%     %contiguous piece of outline
%     clustereddata = cifti_cluster_surfacearea_volume(thissizemap_outline,.5,1.5,0,0,ncortverts,surfacearea,neighbors);
%     
%     %Save that outline to the output file
%     median_size_outlines(:,communitynum) = (clustereddata==1);
%         
%     
% end
% 
% 
% %% Save out results
% 
% %Save probability maps
% out_template.data = outdata;
% out_template.dimord = 'scalar_pos';
% out_template.mapname = columnnames;
% ft_write_cifti_mod([outputfolder '/Cluster_probability_maps_sorted_' num2str(distthresh) 'mm'],out_template)
% 
% %Save outlines
% out_template.data = median_size_outlines;
% ft_write_cifti_mod([outputfolder '/Cluster_probability_maps_sorted_' num2str(distthresh) 'mm_mediansizeoutlines'],out_template)
% 
% %Get the probability maps where more than the defined pct of subjects have
% %clusters in the map
% nsubthresh = ceil(size(subnetworks,2) * pctsubthresh);
% communityindex = logical(numsubs_all >= nsubthresh);
% 
% %Save out those thresholded probability maps
% out_template.data = outdata(:,communityindex);
% out_template.mapname = columnnames(communityindex);
% ft_write_cifti_mod([outputfolder '/Cluster_probability_maps_sorted_' num2str(distthresh) 'mm_' num2str(nsubthresh) 'sub'],out_template)
% 
% %And outlines
% out_template.data = median_size_outlines(:,communityindex);
% ft_write_cifti_mod([outputfolder '/Cluster_probability_maps_sorted_' num2str(distthresh) 'mm_mediansizeoutlines_' num2str(nsubthresh) 'sub'],out_template)
