function [submatches,groupmatches] = networkClusters_match_preclustered_par(subclustermap,subclusternetworkIDs,distances,groupclustermap,groupclusternetworkIDs)
%networkIDs_matchtogroup(subnetworksfilename,[distances],[groupnetworksfile])

distthresh = 20;

%subnetworksfilename = networks(:,1);
%groupnetworksfile = networks(:,2);





nsurfverts = 29696 + 29716;


orig_groupsize = size(groupclustermap,1);
orig_subsize = size(subclustermap,1);

%distances = distances(1:nsurfverts,1:nsurfverts);

subclustermap = subclustermap(1:nsurfverts);

groupclustermap = groupclustermap(1:nsurfverts);
groupIDs = unique(intersect(subclusternetworkIDs(:,2),groupclusternetworkIDs(:,2))); groupIDs(groupIDs<1) = [];

groupmatches = zeros(nsurfverts,0);
submatches = zeros(nsurfverts,0);


for groupIDnum = 1:length(groupIDs)
    groupID = groupIDs(groupIDnum);
    
    %string{groupIDnum} = ['Network #' num2str(groupID)];
    %if groupIDnum==1; fprintf('%s',string{groupIDnum}); else fprintf([repmat('\b',1,length(string{groupIDnum-1})) '%s'],string{groupIDnum}); end
    %if groupIDnum==length(groupIDs); disp(' '); end
    %lowest_mindist = Inf;
    %mindist_ID = 0;
    
    
    subclusternums_thisID = subclusternetworkIDs(subclusternetworkIDs(:,2)==groupID,1);
    subclusters = zeros(size(subclustermap,1),length(subclusternums_thisID));
    for i = 1:size(subclusters,2)
        subclusters(subclustermap==subclusternums_thisID(i),i) = 1;
    end
    
    groupclusternums_thisID = groupclusternetworkIDs(groupclusternetworkIDs(:,2)==groupID,1);
    groupclusters = zeros(size(groupclustermap,1),length(groupclusternums_thisID));
    for i = 1:size(groupclusters,2)
        groupclusters(groupclustermap==groupclusternums_thisID(i),i) = 1;
    end
    
    
    %groupclusters = metric_cluster_cifti(groupnetworks,groupID-.01,groupID+.05,0); groupclusters = groupclusters(1:nsurfverts,:);
    %subclusters = metric_cluster_cifti(subnetworks,groupID-.01,groupID+.05,0); subclusters = subclusters(1:nsurfverts,:);
    
    if (size(subclusters,2) > 0)
    
    mean_mindists = zeros(size(groupclusters,2),size(subclusters,2));
    
    for groupclusnum = 1:size(groupclusters,2)
        for subclusnum = 1:size(subclusters,2)
            
            submindists = min(distances(logical(groupclusters(:,groupclusnum)),logical(subclusters(:,subclusnum))),[],1);
            groupmindists = min(distances(logical(groupclusters(:,groupclusnum)),logical(subclusters(:,subclusnum))),[],2);
            mean_mindists(groupclusnum,subclusnum) = mean([submindists(:); groupmindists(:)]);%mean([submindists(:); groupmindists(:)]));

    
        end
    
    end
    
    best_matches = zeros(size(mean_mindists));
    [ign mini] = min(mean_mindists,[],2);
    groupmatchinds = sub2ind(size(mean_mindists),[1:size(mean_mindists,1)]',mini);
    [ign mini] = min(mean_mindists,[],1);
    submatchinds = sub2ind(size(mean_mindists),mini',[1:size(mean_mindists,2)]');
    best_matchinds = intersect(groupmatchinds,submatchinds);
    best_matches(best_matchinds) = 1;
    
    groupmatchedclusters = find(any(best_matches,2));
    no_groupmatch = find(any(best_matches,2)==0);
    
    for groupclusnum_nomatch = no_groupmatch(:)'
        improvement = zeros(length(groupmatchedclusters),1);
        for groupclusnum_matchi = 1:length(groupmatchedclusters)
            groupclusnum_match = groupmatchedclusters(groupclusnum_matchi);
            matched_subclusnum{groupclusnum_matchi} = find(best_matches(groupclusnum_match,:));
            
            submindists = min(distances(logical(sum(groupclusters(:,[groupclusnum_match groupclusnum_nomatch]),2)),logical(subclusters(:,matched_subclusnum{groupclusnum_matchi}))),[],1);
            groupmindists = min(distances(logical(sum(groupclusters(:,[groupclusnum_match groupclusnum_nomatch]),2)),logical(subclusters(:,matched_subclusnum{groupclusnum_matchi}))),[],2);
            mean_mindist_test = mean([submindists(:); groupmindists(:)]);%mean([submindists(:); groupmindists(:)]));
            
            improvement(groupclusnum_matchi) = mean_mindists(groupclusnum_match,matched_subclusnum{groupclusnum_matchi}) - mean_mindist_test;
            
        end
        
        %if any(improvement > 0)
            
         [max_improvement maxi] = max(improvement);
            
         if (max_improvement > 0) && (nnz(groupclusters(:,groupclusnum_nomatch) .* sum(subclusters(:,matched_subclusnum{maxi}),2)) > 0)
            
            best_matches(groupclusnum_nomatch,matched_subclusnum{maxi}) = 1;
            
            %if mean_mindist_test < mean_mindists(groupclusnum_match,matched_subclusnum)
                
                
                
            %    best_matches(groupclusnum_nomatch,matched_subclusnum) = 1;
                %mean_mindists([groupclusnum_nomatch groupclusnum_match],matched_subclusnum) = mean_mindist_test;
            %end
        end
        clear matched_subclusnum
    end
    
    submatchedclusters = find(any(best_matches,1));
    no_submatch = find(any(best_matches,1)==0);
    
    for subclusnum_nomatch = no_submatch(:)'
        improvement = zeros(length(submatchedclusters),1);
        for subclusnum_matchi = 1:length(submatchedclusters)
            subclusnum_match = submatchedclusters(subclusnum_matchi);
            matched_groupclusnum{subclusnum_matchi} = find(best_matches(:,subclusnum_match));
            
            submindists = min(distances(logical(sum(groupclusters(:,matched_groupclusnum{subclusnum_matchi}),2)),logical(sum(subclusters(:,[subclusnum_match]),2))),[],1);
            groupmindists = min(distances(logical(sum(groupclusters(:,matched_groupclusnum{subclusnum_matchi}),2)),logical(sum(subclusters(:,[subclusnum_match]),2))),[],2);
            mean_mindist_current = mean([submindists(:); groupmindists(:)]);
            
            submindists = min(distances(logical(sum(groupclusters(:,matched_groupclusnum{subclusnum_matchi}),2)),logical(sum(subclusters(:,[subclusnum_match subclusnum_nomatch]),2))),[],1);
            groupmindists = min(distances(logical(sum(groupclusters(:,matched_groupclusnum{subclusnum_matchi}),2)),logical(sum(subclusters(:,[subclusnum_match subclusnum_nomatch]),2))),[],2);
            mean_mindist_test = mean([submindists(:); groupmindists(:)]);%mean([submindists(:); groupmindists(:)]));
            
            improvement(subclusnum_matchi) = mean_mindist_current - mean_mindist_test;
            
        end
        
        %if any(improvement > 0)
            
        [max_improvement maxi] = max(improvement);
        
        if (max_improvement > 0) && (nnz(subclusters(:,subclusnum_nomatch) .* sum(groupclusters(:,matched_groupclusnum{maxi}),2)) > 0)
            
            best_matches(matched_groupclusnum{maxi},subclusnum_nomatch) = 1;
            
%             if mean_mindist_test < mean_mindist_current
%                 
%                 best_matches(matched_groupclusnum,subclusnum_nomatch) = 1;
%                 %mean_mindists(matched_groupclusnum,[subclusnum_match subclusnum_nomatch]) = mean_mindist_test;
%             end
        end
        clear matched_groupclusnum
    end
    
    clusterstowrite = any(best_matches,2);
    for groupclusnum = 1:size(best_matches,1)
        if logical(clusterstowrite(groupclusnum))
            
            matchedsubclusnums = logical(best_matches(groupclusnum,:));
            matchedgroupclusnums = logical(any(best_matches(:,matchedsubclusnums),2));
            clusterstowrite(matchedgroupclusnums) = 0;
            
            submindists = min(distances(logical(any(groupclusters(:,matchedgroupclusnums),2)),logical(any(subclusters(:,matchedsubclusnums),2))),[],1);
            groupmindists = min(distances(logical(any(groupclusters(:,matchedgroupclusnums),2)),logical(any(subclusters(:,matchedsubclusnums),2))),[],2);
            this_mean_mindist = mean([submindists(:); groupmindists(:)]);
            
            if this_mean_mindist <= distthresh
                groupmatches = [groupmatches logical(any(groupclusters(:,matchedgroupclusnums),2))];
                submatches = [submatches logical(any(subclusters(:,matchedsubclusnums),2))];
            end
                        
        end
    end
    end
end

groupmatches(end+1:orig_groupsize,:) = 0;
submatches(end+1:orig_subsize,:) = 0;
            
                
% 
%     
%     
%     
%     
%     
%     mean_mindists(mean_mindists>thresh) = Inf;
%     
%     [subtogroup_assign,grouptosub_assign,subtogroup_cost, grouptosub_cost] = munkres_mult(mean_mindists,5);
%     
%     
%     unique_subtogroup = unique(subtogroup_assign); unique_subtogroup(unique_subtogroup==0) = [];
%     thisnetwork_groupmatchclusters = zeros(size(groupclusters,1),length(unique_subtogroup));
%     thisnetwork_submatchclusters = zeros(size(groupclusters,1),length(unique_subtogroup));
%     for i = 1:length(unique_subtogroup)
%         thisnetwork_groupmatchclusters(:,i) = sum(groupclusters(:,logical(subtogroup_assign==unique_subtogroup(i))),2);
%         groupclusters_inthiscluster{i} = find(subtogroup_assign==unique_subtogroup(i));
%     end
%     unique_grouptosub = unique(grouptosub_assign); unique_grouptosub(unique_grouptosub==0) = [];
%     for j = 1:length(unique_grouptosub)
%         for i = 1:length(groupclusters_inthiscluster)
%             if any(groupclusters_inthiscluster{i}==unique_grouptosub(j))
%                 targetcluster = i;
%                 break
%             end
%         end
%         thisnetwork_submatchclusters(:,targetcluster) = sum(subclusters(:,logical(grouptosub_assign==unique_grouptosub(j))),2);
%     end
%     clear groupclusters_inthiscluster
%     
%     inds_toremove = logical(all(thisnetwork_submatchclusters==0,1) + all(thisnetwork_groupmatchclusters==0,1));
%     thisnetwork_groupmatchclusters(:,inds_toremove) = [];
%     thisnetwork_submatchclusters(:,inds_toremove) = [];
%     
%     groupmatches(:,(end+1):(end+size(thisnetwork_groupmatchclusters,2))) = thisnetwork_groupmatchclusters;
%     submatches(:,(end+1):(end+size(thisnetwork_submatchclusters,2))) = thisnetwork_submatchclusters;
%     end
%     
% end
%     
% 
% %out = zeros(66697,size(matched,2)); out(1:nsurfverts,:) = matched;
% %cifti_write_wHDR(out,[],[subnetworksfolder '/Clustermatch'])
% 
