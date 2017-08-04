function Distance_from_infomap_borders_mm(giftispace_infomap,edges,thresh,hem,outputname)

load(['/data/cn4/evan/fsaverage_LR32k/Surface_distances_' hem '.mat'])

bufsize=16384;
%caretdir = ['/data/cn4/laumannt/assignment_problem_v2/caret2/PALS_B12.BOTH-HEMS.CLEAN.73730/LEFT/18K'];
caretdir = '/data/cn4/evan/fsaverage_LR32k/';
%cd(caretdir)
% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir 'node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;

edgesthresh = logical(edges >= thresh);

netborders = zeros(size(giftispace_infomap));

for i = 1:length(netborders)
    indneighs = neighbors(i,2:end); indneighs(isnan(indneighs)) = [];
    neighvals = giftispace_infomap(indneighs); neighvals(neighvals<1) = [];
    if any(neighvals~=giftispace_infomap(i))
        
        dist = min(geo_distances(i,edgesthresh));
        
        netborders(i) = dist;
        
    else
        netborders(i) = -1;
    end
end

save(gifti(single(netborders)),outputname);