function [groupmatches, submatches] = networkClusters_matchtogroup(subnetworksfilename,distances,groupnetworksfile)
%networkIDs_matchtogroup(subnetworksfilename,[distances],[groupnetworksfile])

thresh = 10;

if ischar(subnetworksfilename)
    slashlocs = strfind(subnetworksfilename,'/');
    if isempty(slashlocs)
        subnetworksfolder = pwd;
        subnetworksfile = subnetworksfilename;
    else
        subnetworksfolder = subnetworksfilename(1:slashlocs(end));
        subnetworksfile = subnetworksfilename(slashlocs(end)+1 : end);
    end
    subnetworks = cifti_read([subnetworksfolder '/' subnetworksfile]);
else
    subnetworks = subnetworksfilename;
    clear subnetworksfilename
end

if ~exist('groupnetworksfile')
    groupnetworksfile = '/data/cn4/evan/RestingState/Consensus/120_LR_minsize400_recolored_manualconsensus_LR.dtseries.nii';
elseif ischar(groupnetworksfile)
    groupnetworks = cifti_read(groupnetworksfile);
else
    groupnetworks = groupnetworksfile;
    clear groupnetworksfile
end

nsurfverts = 29696 + 29716;


if ~exist('distances') || isempty(distances)
    load /data/cn4/evan/fsaverage_LR32k/Cifti_geo_distances.mat
    distances = distances(1:nsurfverts,1:nsurfverts);
end



subnetworks = subnetworks(1:nsurfverts);

groupnetworks = groupnetworks(1:nsurfverts);
groupIDs = unique(groupnetworks); groupIDs(groupIDs<1) = [];

groupmatches = zeros(length(groupnetworks),0);
submatches = zeros(length(groupnetworks),0);


for groupIDnum = 1:length(groupIDs)
    groupID = groupIDs(groupIDnum);
    
    string{groupIDnum} = ['Network #' num2str(groupID)];
    if groupIDnum==1; fprintf('%s',string{groupIDnum}); else fprintf([repmat('\b',1,length(string{groupIDnum-1})) '%s'],string{groupIDnum}); end
    if groupIDnum==length(groupIDs); disp(' '); end
    %lowest_mindist = Inf;
    %mindist_ID = 0;
    
    groupclusters = metric_cluster_cifti(groupnetworks,groupID-.01,groupID+.05,0); groupclusters = groupclusters(1:nsurfverts,:);
    subclusters = metric_cluster_cifti(subnetworks,groupID-.01,groupID+.05,0); subclusters = subclusters(1:nsurfverts,:);
    
    if size(subclusters,2) > 0
    
    mean_mindists = zeros(size(groupclusters,2),size(subclusters,2));
    
    for groupclusnum = 1:size(groupclusters,2)
        for subclusnum = 1:size(subclusters,2)
            
            submindists = min(distances(logical(groupclusters(:,groupclusnum)),logical(subclusters(:,subclusnum))),[],1);
            groupmindists = min(distances(logical(groupclusters(:,groupclusnum)),logical(subclusters(:,subclusnum))),[],2);
            mean_mindists(groupclusnum,subclusnum) = min(mean(submindists),mean(groupmindists));%mean([submindists(:); groupmindists(:)]));

    
        end
    
    end
    
    mean_mindists(mean_mindists>thresh) = Inf;
    
    [subtogroup_assign,grouptosub_assign,subtogroup_cost, grouptosub_cost] = munkres_mult(mean_mindists,5);
    
    
    unique_subtogroup = unique(subtogroup_assign); unique_subtogroup(unique_subtogroup==0) = [];
    thisnetwork_groupmatchclusters = zeros(size(groupclusters,1),length(unique_subtogroup));
    thisnetwork_submatchclusters = zeros(size(groupclusters,1),length(unique_subtogroup));
    for i = 1:length(unique_subtogroup)
        thisnetwork_groupmatchclusters(:,i) = sum(groupclusters(:,logical(subtogroup_assign==unique_subtogroup(i))),2);
        groupclusters_inthiscluster{i} = find(subtogroup_assign==unique_subtogroup(i));
    end
    unique_grouptosub = unique(grouptosub_assign); unique_grouptosub(unique_grouptosub==0) = [];
    for j = 1:length(unique_grouptosub)
        for i = 1:length(groupclusters_inthiscluster)
            if any(groupclusters_inthiscluster{i}==unique_grouptosub(j))
                targetcluster = i;
                break
            end
        end
        thisnetwork_submatchclusters(:,targetcluster) = sum(subclusters(:,logical(grouptosub_assign==unique_grouptosub(j))),2);
    end
    clear groupclusters_inthiscluster
    
    inds_toremove = logical(all(thisnetwork_submatchclusters==0,1) + all(thisnetwork_groupmatchclusters==0,1));
    thisnetwork_groupmatchclusters(:,inds_toremove) = [];
    thisnetwork_submatchclusters(:,inds_toremove) = [];
    
    groupmatches(:,(end+1):(end+size(thisnetwork_groupmatchclusters,2))) = thisnetwork_groupmatchclusters;
    submatches(:,(end+1):(end+size(thisnetwork_submatchclusters,2))) = thisnetwork_submatchclusters;
    end
    
end
    

%out = zeros(66697,size(matched,2)); out(1:nsurfverts,:) = matched;
%cifti_write_wHDR(out,[],[subnetworksfolder '/Clustermatch'])

