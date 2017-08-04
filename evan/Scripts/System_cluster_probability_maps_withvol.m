subnetworksfile = ['Templatematch_dice_bysubject_kden0.05.dtseries.nii'];
maxdicesfile = ['Templatematch_dice_bysubject_kden0.05_dicemaxes.dtseries.nii'];
dicediffsfile = ['Templatematch_dice_bysubject_kden0.05_dicediffs.dtseries.nii'];

calc_dices = 1;


subnetworks = ft_read_cifti_mod(subnetworksfile); 
out_template = subnetworks; out_template.data = [];
subnetworks = subnetworks.data;

if calc_dices
    maxdices = ft_read_cifti_mod(maxdicesfile); maxdices = maxdices.data;
    dicediffs = ft_read_cifti_mod(dicediffsfile); dicediffs = dicediffs.data;
end

IDs = unique(subnetworks); IDs(IDs<1) = [];

networklabels = {'Default','Visual','FrontoPar','FrontoPar 2','DorsalAttn','DorsalAttn2','VentAttn','Salience','CingOperc','MortorHand','MotorMouth','Auditory','MTL1','MTL2','MedPar','ParOccip'};

clusters_subs_IDs_dices = zeros(0,4);
subnetworks_bycluster = zeros(size(subnetworks));

neighbors = cifti_neighbors(out_template);

prevstring = [];
clustercount = 0;
for s = 1:size(subnetworks,2)
    string = ['Clustering systems: Subject ' num2str(s)];
    fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
    prevstring = string;
    for ID = IDs(:)'
        outputcifti = cifti_cluster_surfacearea_volume(subnetworks(:,s),ID-.05,ID+.05,0,0,59412,ones(length(subnetworks(:,s)),1),neighbors);
        clusterIDs = unique(outputcifti);
        for clusternum = 1:length(clusterIDs)
            clustercount = clustercount + 1;
            
            if calc_dices
                this_maxdice = mean(maxdices(outputcifti==clusterIDs(clusternum),s));
                this_dicediff = mean(dicediffs(outputcifti==clusterIDs(clusternum),s));
            else
                this_maxdice = 0;
                this_dicediff = 0;
            end
            
            clusters_subs_IDs_dices(clustercount,:) = [s ID this_maxdice this_dicediff];
            subnetworks_bycluster(outputcifti==clusterIDs(clusternum),s) = clustercount;
        end
    end
end
disp(' ')

%%

distances = smartload('/data/cn4/evan/fsaverage_LR32k/Cifti_geo_distances.mat');

indices = {[]};
count = 0;

prevstring = [];
for s1 = 1:size(subnetworks,2)
    for s2 = (s1+1):size(subnetworks,2)
        
        
        string = ['Matching subjects pairwise: Subject ' num2str(s1) ' vs Subject ' num2str(s2) ', '];
        fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
        prevstring = [string 'Network #' num2str(max(IDs)) '  '];
        
        [s1matches,s2matches] = networkClusters_match_preclustered_withvol(subnetworks_bycluster(:,s1),[find(clusters_subs_IDs_dices(:,1)==s1) clusters_subs_IDs_dices(clusters_subs_IDs_dices(:,1)==s1,2)],distances,subnetworks_bycluster(:,s2),[find(clusters_subs_IDs_dices(:,1)==s2) clusters_subs_IDs_dices(clusters_subs_IDs_dices(:,1)==s2,2)],15);
        
        for matchnum = 1:size(s1matches,2)
            
            s1_clusternums = unique(subnetworks_bycluster(logical(s1matches(:,matchnum)),s1));
            s2_clusternums = unique(subnetworks_bycluster(logical(s2matches(:,matchnum)),s2));
            
            count = count+1;
            indices{count,1} = s1_clusternums;
            indices{count,2} = s2_clusternums;
            
            
        end
        
    end
end
disp(' ')
clear distances

match_matrix = zeros(clustercount);
for i = 1:size(indices,1)
    match_matrix(indices{i,1},indices{i,2}) = 1;  match_matrix(indices{i,2},indices{i,1}) = 1;
end

save('match_matrix.mat','match_matrix','-v7.3')

%%


mat2pajek_mod_EG(match_matrix,find(triu(match_matrix,1)),'pajek.net')
rawclrs = infomap_wrapper_mod_EG('pajek.net',100);
dlmwrite(['rawassn.txt'],rawclrs,'\t')
delete('pajek.clu')
simple_assigns = modify_clrfile('simplify','rawassn.txt',2);

communities = unique(simple_assigns); communities(communities<1) = [];

outdata = zeros(size(subnetworks,1),length(communities));

for communitynum = 1:length(communities)
    community = communities(communitynum);
    
    clusternums_in_community = find(simple_assigns==community);
    
    community_IDs(communitynum) = mode(clusters_subs_IDs_dices(clusternums_in_community,2));
end

[sortedIDs sorti] = sort(community_IDs,'ascend');

sorted_communities = communities(sorti);

for communitynum = 1:length(sorted_communities)
    community = sorted_communities(communitynum);
    
    clusternums_in_community = find(simple_assigns==community);
    
    this_ID = mode(clusters_subs_IDs_dices(clusternums_in_community,2));
    
    these_subjects = unique(clusters_subs_IDs_dices(clusternums_in_community,1));
    
    numsubs = length(these_subjects);
    
    clusternums_not_in_community = find(simple_assigns~=community);
    additionalsubs = [];
    for clusternum = clusternums_not_in_community(:)'
        if ~any(these_subjects==clusters_subs_IDs_dices(clusternum,1))
            matches_to_thiscommunity = clusternums_in_community(logical(match_matrix(clusternum,clusternums_in_community)));
            if length(unique(clusters_subs_IDs_dices(matches_to_thiscommunity,1))) > (numsubs/2)
                additionalsubs = [additionalsubs clusters_subs_IDs_dices(clusternum,1)];
            end
        end
    end
    numsubs_additional = numsubs + length(unique(additionalsubs));
    
    avg_maxdice = mean(clusters_subs_IDs_dices(clusternums_in_community,3));
    avg_maxdicestr = sprintf('%0.3f',avg_maxdice); avg_maxdicestr = avg_maxdicestr(2:end);
    avg_dicediff = mean(clusters_subs_IDs_dices(clusternums_in_community,4));
    avg_dicediffstr = sprintf('%0.3f',avg_dicediff); avg_dicediffstr = avg_dicediffstr(2:end);
    
    columnnames{1,communitynum} = [networklabels{this_ID} ': ' num2str(numsubs) ' of ' num2str(size(subnetworks,2)) ' subs; ' num2str(numsubs_additional) ' subs total; Dice=' avg_maxdicestr ', Dice diff=' avg_dicediffstr];
    
    for this_clusternum = clusternums_in_community(:)'
        
        outdata(:,communitynum) = outdata(:,communitynum) + ((subnetworks_bycluster(:,clusters_subs_IDs_dices(this_clusternum,1))==this_clusternum) ./ numsubs);
    
    end
end


out_template.data = outdata;
out_template.dimord = 'scalar_pos';
out_template.mapname = columnnames;