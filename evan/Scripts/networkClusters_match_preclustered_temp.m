function [groupmatches,submatches,groupclus_dmat] = networkClusters_match_preclustered_temp(groupclustermap,groupclusternetworkIDs, distances, groupclus_dmat,subclustermap,subclusternetworkIDs,distthresh)
%[groupmatches,submatches,groupclus_dmat] = networkClusters_match_preclustered_temp(groupclustermap,groupclusternetworkIDs, distances, groupclus_dmat,subclustermap,subclusternetworkIDs,distthresh)

if ~exist('distthresh')
    distthresh = 20;
end


makenew_groupclus_dmat = 0;
if isempty(groupclus_dmat)
    makenew_groupclus_dmat = 1;
end


nsurfverts = 29696 + 29716;


orig_groupsize = size(groupclustermap,1);
orig_subsize = size(subclustermap,1);

%distances = distances(1:nsurfverts,1:nsurfverts);

subclustermap = subclustermap(1:nsurfverts);

groupclustermap = groupclustermap(1:nsurfverts);
groupIDs = unique(intersect(subclusternetworkIDs(:,2),groupclusternetworkIDs(:,2))); groupIDs(groupIDs<1) = [];

groupmatches = zeros(nsurfverts,0);
submatches = zeros(nsurfverts,0);


if makenew_groupclus_dmat; groupclus_dmat = cell(length(groupIDs),1); end

for groupIDnum = 1:length(groupIDs)
    groupID = groupIDs(groupIDnum);
    
    string{groupIDnum} = ['Network #' num2str(groupID)];
    if groupIDnum==1; fprintf('%s',string{groupIDnum}); else fprintf([repmat('\b',1,length(string{groupIDnum-1})) '%s'],string{groupIDnum}); end
    if groupIDnum==length(groupIDs); disp(' '); end
    
    
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
    
    
    if (size(subclusters,2) > 0)
    
    mean_mindists = zeros(size(groupclusters,2),size(subclusters,2));
    
    if makenew_groupclus_dmat; groupclus_dmat{groupIDnum} = cell(size(groupclusters,2),1); end
    for groupclusnum = 1:size(groupclusters,2)
        if makenew_groupclus_dmat; groupclus_dmat{groupIDnum}{groupclusnum} = distances(logical(groupclusters(:,groupclusnum)),:); end
        for subclusnum = 1:size(subclusters,2)
            
            intersect_dmat = groupclus_dmat{groupIDnum}{groupclusnum}(:,logical(subclusters(:,subclusnum)));
            
            submindists = min(intersect_dmat,[],1);
            groupmindists = min(intersect_dmat,[],2);
            
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
            
            intersect_dmat = [groupclus_dmat{groupIDnum}{groupclusnum_match}(:,logical(subclusters(:,matched_subclusnum{groupclusnum_matchi}))),; groupclus_dmat{groupIDnum}{groupclusnum_nomatch}(:,logical(subclusters(:,matched_subclusnum{groupclusnum_matchi})))];
            submindists = min(intersect_dmat,[],1);
            groupmindists = min(intersect_dmat,[],2);
            %submindists = min(distances(logical(sum(groupclusters(:,[groupclusnum_match groupclusnum_nomatch]),2)),logical(subclusters(:,matched_subclusnum{groupclusnum_matchi}))),[],1);
            %groupmindists = min(distances(logical(sum(groupclusters(:,[groupclusnum_match groupclusnum_nomatch]),2)),logical(subclusters(:,matched_subclusnum{groupclusnum_matchi}))),[],2);
            mean_mindist_test = mean([submindists(:); groupmindists(:)]);%mean([submindists(:); groupmindists(:)]));
            
            improvement(groupclusnum_matchi) = mean_mindists(groupclusnum_match,matched_subclusnum{groupclusnum_matchi}) - mean_mindist_test;
            
        end
        
            
         [max_improvement maxi] = max(improvement);
            
         if (max_improvement > 0) && (nnz(groupclusters(:,groupclusnum_nomatch) .* sum(subclusters(:,matched_subclusnum{maxi}),2)) > 0)
            
            best_matches(groupclusnum_nomatch,matched_subclusnum{maxi}) = 1;
            
           
        end
        clear matched_subclusnum
    end
    
    %clear groupclus_dmat
    submatchedclusters = find(any(best_matches,1));
    no_submatch = find(any(best_matches,1)==0);
    
    for subclusnum_nomatch = no_submatch(:)'
        improvement = zeros(length(submatchedclusters),1);
        for subclusnum_matchi = 1:length(submatchedclusters)
            subclusnum_match = submatchedclusters(subclusnum_matchi);
            matched_groupclusnum{subclusnum_matchi} = find(best_matches(:,subclusnum_match));
            
            current_intersect_dmat = distances(logical(sum(groupclusters(:,matched_groupclusnum{subclusnum_matchi}),2)),logical(sum(subclusters(:,[subclusnum_match]),2)));
            
            submindists = min(current_intersect_dmat,[],1);
            groupmindists = min(current_intersect_dmat,[],2);
            mean_mindist_current = mean([submindists(:); groupmindists(:)]);
            
            test_intersect_dmat = [current_intersect_dmat distances(logical(sum(groupclusters(:,matched_groupclusnum{subclusnum_matchi}),2)),logical(sum(subclusters(:,[subclusnum_nomatch]),2)))];
            
            submindists = min(test_intersect_dmat,[],1);
            groupmindists = min(test_intersect_dmat,[],2);
            mean_mindist_test = mean([submindists(:); groupmindists(:)]);%mean([submindists(:); groupmindists(:)]));
            
            improvement(subclusnum_matchi) = mean_mindist_current - mean_mindist_test;
            
        end
        
            
        [max_improvement maxi] = max(improvement);
        
        if (max_improvement > 0) && (nnz(subclusters(:,subclusnum_nomatch) .* sum(groupclusters(:,matched_groupclusnum{maxi}),2)) > 0)
            
            best_matches(matched_groupclusnum{maxi},subclusnum_nomatch) = 1;
            
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
            
                

