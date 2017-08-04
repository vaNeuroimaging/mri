function [matches, groupclus_dmat] = Hubs_match_preclustered3(groupclustermap,groupclusternetworkIDs, distances, groupclus_dmat,subclustermap,subclusternetworkIDs,distthresh,common_connections_thresh)
%[groupmatches,submatches,groupclus_dmat] = networkClusters_match_preclustered(groupclustermap,groupclusternetworkIDs, distances, groupclus_dmat,subclustermap,subclusternetworkIDs,distthresh,common_connections_thresh)

if ~exist('distthresh') || isempty(distthresh)
    distthresh = 10;
end

common_networks = intersect(vertcat(groupclusternetworkIDs{:}),vertcat(subclusternetworkIDs{:}));


nsurfverts = 29696 + 29716;


%orig_groupsize = size(groupclustermap,1);
%orig_subsize = size(subclustermap,1);

subclustermap = subclustermap(1:nsurfverts);

groupclustermap = groupclustermap(1:nsurfverts);
%groupIDs = unique(intersect(subclusternetworkIDs(:,2),groupclusternetworkIDs(:,2))); groupIDs(groupIDs<1) = [];
%groupIDs = unique(groupclusternetworkIDs(:,2)); groupIDs(groupIDs<1) = [];

%groupmatches = zeros(nsurfverts,0);
%submatches = zeros(nsurfverts,0);

groupclusterIDs = unique(groupclustermap); groupclusterIDs(groupclusterIDs<1) = [];
if isempty(groupclus_dmat)
    groupclus_dmat = cell(length(groupclusterIDs),1);
    for i = 1:length(groupclusterIDs)
        groupclus_dmat{i} = distances(logical(groupclustermap==groupclusterIDs(i)),:);
    end
end

subclusterIDs = unique(subclustermap); subclusterIDs(subclusterIDs<1) = [];
% subclus_dmat = cell(length(subclusterIDs),1);
% for i = 1:length(groupclusterIDs)
%     subclus_dmat{i} = distances(logical(subclustermap==subclusterIDs(i)),:);
% end


mean_mindists = zeros(length(groupclusterIDs),length(subclusterIDs));

matches = zeros(0,2);

for groupclusnum = 1:length(groupclusterIDs)
    for subclusnum = 1:length(subclusterIDs)
        
        commonconnections = intersect(intersect(groupclusternetworkIDs{groupclusnum},subclusternetworkIDs{subclusnum}),common_networks);
        if (length(commonconnections) >= common_connections_thresh) || ...
            ((length(commonconnections) == length(intersect(groupclusternetworkIDs{groupclusnum},common_networks))) && ...
           (length(commonconnections) == length(intersect(subclusternetworkIDs{subclusnum},common_networks))))
            
            
            intersect_dmat = groupclus_dmat{groupclusnum}(:,subclustermap==subclusterIDs(subclusnum));
            
            submindists = min(intersect_dmat,[],1);
            groupmindists = min(intersect_dmat,[],2);
            
            mean_mindists = mean([submindists(:); groupmindists(:)]);%mean([submindists(:); groupmindists(:)]));
            
            
            if mean_mindists <= distthresh
                
                matches(end+1,:) = [groupclusterIDs(groupclusnum) subclusterIDs(subclusnum)];
                
                
            end
            
            
            
            
        end
        
        
    end
    
end


            
                

