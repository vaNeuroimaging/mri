
subnetworksfile = ['Templatematch_dice_bysubject_kden0.05.dscalar.nii'];
maxdicesfile = ['Templatematch_dice_bysubject_kden0.05_dicemaxes.dscalar.nii'];
dicediffsfile = ['Templatematch_dice_bysubject_kden0.05_dicediffs.dscalar.nii'];


distthresh = 10;

pctsubthresh = 1/3;

new_cluster = 0;

calc_dices = 0;


subnetworks = ft_read_cifti_mod(subnetworksfile); 
out_template = subnetworks; out_template.data = [];
subnetworks = subnetworks.data;

if calc_dices
    maxdices = ft_read_cifti_mod(maxdicesfile); maxdices = maxdices.data;
    dicediffs = ft_read_cifti_mod(dicediffsfile); dicediffs = dicediffs.data;
end

IDs = unique(subnetworks); IDs(IDs<1) = [];

networklabels = {'Default','Visual','FrontoPar','FrontoPar 2','DorsalAttn','DorsalAttn2','VentAttn','Salience','CingOperc','MortorHand','MotorMouth','Auditory','MTL1','MTL2','MedPar','ParOccip'};
%networklabels = {'Default','Visual','FrontoPar','FrontoPar 2','DorsalAttn','DorsalAttn2','VentAttn','Salience','CingOperc','MortorHand','MotorMouth','Auditory','MTL1','STS','MedPar','ParOccip','Orbitofrontal','LatVisual'};



if new_cluster
    
    surfacearea = ft_read_cifti_mod('/data/cn4/evan/fsaverage_LR32k/Conte69.LR.midthickness.32k_fs_LR_surfaceareas_normalwall.dtseries.nii');
    surfacearea = surfacearea.data;
    
    clusters_subs_IDs_dices_SAs = zeros(0,5);
    subnetworks_bycluster = zeros(size(subnetworks));
    
    prevstring = [];
    clustercount = 0;
    for s = 1:size(subnetworks,2)
        string = ['Clustering systems: Subject ' num2str(s)];
        fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
        prevstring = string;
        
        for ID = IDs(:)'
            outputcifti = metric_cluster_cifti(subnetworks(:,s),ID-.05,ID+.05,0);
            for clusternum = 1:size(outputcifti,2)
                clustercount = clustercount + 1;
                
                if calc_dices
                    this_maxdice = mean(maxdices(logical(outputcifti(:,clusternum)),s));
                    this_dicediff = mean(dicediffs(logical(outputcifti(:,clusternum)),s));
                else
                    this_maxdice = 0;
                    this_dicediff = 0;
                end
                
                cluster_surfacearea = sum(surfacearea(logical(outputcifti(:,clusternum))));
                
                clusters_subs_IDs_dices_SAs(clustercount,:) = [s ID this_maxdice this_dicediff cluster_surfacearea];
                subnetworks_bycluster(logical(outputcifti(:,clusternum)),s) = clustercount;
            end
        end
    end
    disp(' ')
    
    save('System_clustering.mat','clusters_subs_IDs_dices_SAs','subnetworks_bycluster','clustercount','-v7.3')
else
    load('System_clustering.mat')
end

%%

distances = smartload('/data/cn4/evan/fsaverage_LR32k/Cifti_geo_distances_uint8.mat');
%distances = uint8(distances);

indices = {[]};
count = 0;

prevstring = [];
for s1 = 1:size(subnetworks,2)
    s1clus_dmat = [];
    for s2 = (s1+1):size(subnetworks,2)
        
       
            string = ['Matching subjects pairwise: Subject ' num2str(s1) ' vs Subject ' num2str(s2) ', '];
            fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
            prevstring = [string 'Network #' num2str(max(IDs)) '  '];
            
            [s1matches,s2matches,s1clus_dmat] = networkClusters_match_preclustered(subnetworks_bycluster(:,s1),[find(clusters_subs_IDs_dices_SAs(:,1)==s1) clusters_subs_IDs_dices_SAs(clusters_subs_IDs_dices_SAs(:,1)==s1,2)],distances,s1clus_dmat,subnetworks_bycluster(:,s2),[find(clusters_subs_IDs_dices_SAs(:,1)==s2) clusters_subs_IDs_dices_SAs(clusters_subs_IDs_dices_SAs(:,1)==s2,2)],distthresh);
            
            for matchnum = 1:size(s1matches,2)
                
                s1_clusternums = unique(subnetworks_bycluster(logical(s1matches(:,matchnum)),s1));
                s2_clusternums = unique(subnetworks_bycluster(logical(s2matches(:,matchnum)),s2));
                
                count = count+1;
                indices{count,1} = s1_clusternums;
                indices{count,2} = s2_clusternums;
                
                %match_matrix(s1_clusternums,s2_clusternums) = 1;  match_matrix(s2_clusternums,s1_clusternums) = 1;
                
            end
        
    end
end
disp(' ')
clear distances

match_matrix = false(clustercount);
for i = 1:size(indices,1)
    match_matrix(indices{i,1},indices{i,2}) = true;  match_matrix(indices{i,2},indices{i,1}) = true;
end

save('match_matrix.mat','match_matrix','-v7.3')

%%

surfacearea = ft_read_cifti_mod('/data/cn4/evan/fsaverage_LR32k/Conte69.LR.midthickness.32k_fs_LR_surfaceareas_normalwall.dtseries.nii');
surfacearea = surfacearea.data;
neighbors = cifti_neighbors('/data/cn4/evan/fsaverage_LR32k/Conte69.LR.midthickness.32k_fs_LR_surfaceareas_normalwall.dtseries.nii');

mat2pajek_byindex(match_matrix,find(triu(match_matrix,1)),'pajek.net')
rawclrs = run_infomap_on_pajekfile('pajek.net',100);
dlmwrite(['rawassn.txt'],rawclrs,'\t')
delete('pajek.clu')
simple_assigns = modify_clrfile('simplify','rawassn.txt',2);
%%
communities = unique(simple_assigns); communities(communities<1) = [];

outdata = zeros(size(subnetworks,1),length(communities));
median_size_outlines = zeros(size(subnetworks,1),length(communities));

for communitynum = 1:length(communities)
    community = communities(communitynum);
    
    clusternums_in_community = find(simple_assigns==community);
    
    community_IDs(communitynum) = mode(clusters_subs_IDs_dices_SAs(clusternums_in_community,2));
end

[sortedIDs sorti] = sort(community_IDs,'ascend');

sorted_communities = communities(sorti);

numsubs_all = zeros(length(sorted_communities),1);

for communitynum = 1:length(sorted_communities)
    community = sorted_communities(communitynum);
    
    clusternums_in_community = find(simple_assigns==community);
    
    this_ID = mode(clusters_subs_IDs_dices_SAs(clusternums_in_community,2));
    
    these_subjects = unique(clusters_subs_IDs_dices_SAs(clusternums_in_community,1));
    
    numsubs = length(these_subjects);
    
    clusternums_not_in_community = find(simple_assigns~=community);
    additionalsubs = [];
    for clusternum = clusternums_not_in_community(:)'
        if ~any(these_subjects==clusters_subs_IDs_dices_SAs(clusternum,1))
            matches_to_thiscommunity = clusternums_in_community(logical(match_matrix(clusternum,clusternums_in_community)));
            if length(unique(clusters_subs_IDs_dices_SAs(matches_to_thiscommunity,1))) > (numsubs/2)
                additionalsubs = [additionalsubs clusters_subs_IDs_dices_SAs(clusternum,1)];
            end
        end
    end
    numsubs_all(communitynum) = numsubs + length(unique(additionalsubs));
    
    avg_maxdice = mean(clusters_subs_IDs_dices_SAs(clusternums_in_community,3));
    avg_maxdicestr = sprintf('%0.3f',avg_maxdice); avg_maxdicestr = avg_maxdicestr(2:end);
    avg_dicediff = mean(clusters_subs_IDs_dices_SAs(clusternums_in_community,4));
    avg_dicediffstr = sprintf('%0.3f',avg_dicediff); avg_dicediffstr = avg_dicediffstr(2:end);
    
    combined_cluster_sizes = zeros(numsubs,1);
    for subnum = 1:numsubs
        thissub = these_subjects(subnum);
        thissub_clusternums = find(clusters_subs_IDs_dices_SAs(:,1)==thissub);
        sub_sizes = clusters_subs_IDs_dices_SAs(intersect(thissub_clusternums,clusternums_in_community),5);
        combined_cluster_sizes(subnum) = sum(sub_sizes);
    end
    
    median_size = median(combined_cluster_sizes);
    
    
    columnnames{1,communitynum} = [networklabels{this_ID} ': ' num2str(numsubs) ' of ' num2str(size(subnetworks,2)) ' subs; ' num2str(numsubs_all(communitynum)) ' subs total; median size=' num2str(round(median_size)) 'mm2; Dice=' avg_maxdicestr ', Dice diff=' avg_dicediffstr];
    
    for this_clusternum = clusternums_in_community(:)'
        
        outdata(:,communitynum) = outdata(:,communitynum) + ((subnetworks_bycluster(:,clusters_subs_IDs_dices_SAs(this_clusternum,1))==this_clusternum) ./ numsubs);
    
    end
    
    
    
    thismap = outdata(:,communitynum);
    
    [maxval maxi] = max(thismap);
    thissizemap = false(size(thismap));
    thissizemap(maxi) = 1;
    sortedvals = sort(unique(thismap),'descend');
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
        if sum(surfacearea(logical(thissizemap))) >= median_size;
           break
        end
    end
    
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
    
    outputclusters = cifti_cluster_surfacearea(thissizemap_outline,.5,1.5,0,neighbors);
    
    median_size_outlines(:,communitynum) = outputclusters(:,1);
    
    
    
end


out_template.data = outdata;
out_template.dimord = 'scalar_pos';
out_template.mapname = columnnames;
ft_write_cifti_mod(['Cluster_probability_maps_sorted_' num2str(distthresh) 'mm'],out_template)


out_template.data = median_size_outlines;
ft_write_cifti_mod(['Cluster_probability_maps_sorted_' num2str(distthresh) 'mm_mediansizeoutlines'],out_template)


nsubthresh = ceil(size(subnetworks,2) * pctsubthresh);
communityindex = logical(numsubs_all >= nsubthresh);

out_template.data = outdata(:,communityindex);
out_template.mapname = columnnames(communityindex);
ft_write_cifti_mod(['Cluster_probability_maps_sorted_' num2str(distthresh) 'mm_' num2str(nsubthresh) 'sub'],out_template)

out_template.data = median_size_outlines(:,communityindex);
ft_write_cifti_mod(['Cluster_probability_maps_sorted_' num2str(distthresh) 'mm_mediansizeoutlines_' num2str(nsubthresh) 'sub'],out_template)
