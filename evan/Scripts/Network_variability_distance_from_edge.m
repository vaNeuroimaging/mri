variabilitymapname = 'Subject_clustering_vertexwise_allnetworks_L_minQ_0.1.func.gii';
variabilitymap = gifti(variabilitymapname); variabilitymap = variabilitymap.cdata(:,1:4);

minsize = 4;

bufsize=16384;
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';

% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;

% for vertex = 1:size(variabilitymap,1)
%     vertneighs = neighbors(vertex,2:7);
%     vertneighs(isnan(vertneighs)) = [];
%     if ~any(variabilitymap(vertneighs,1)==variabilitymap(vertex,1))
%         variabilitymap(vertex,1) = mode(variabilitymap(vertneighs,1));
%     end
% end

networkIDs = unique(variabilitymap(variabilitymap>0));
giftimap = variabilitymap(:,1);

for networkID = networkIDs'
    clusteredmetric = zeros(size(giftimap));
    thiscolorverts = find(giftimap==networkID);
    for vertex = thiscolorverts'
        
        %find the neighbors of this vertex
        vertexneighbors = neighbors(vertex,:);
        
        %find which of those neighbors also pass the thresholds
        vertexneighbors_thiscolor = intersect(thiscolorverts,vertexneighbors);
        
        %find if those neighbors have already been assigned different cluster values
        uniqueneighborvals = unique(clusteredmetric(vertexneighbors_thiscolor));
        uniqueneighborvals(uniqueneighborvals==0) = [];
        
        %if no neighbors have cluster identifiers, assign them the number of this vertex as a unique cluster identifier
        if isempty(uniqueneighborvals)
            clusteredmetric(vertexneighbors_thiscolor) = vertex;
            %if there is only one previous cluster identifier present, make all the neighbors that value
        elseif length(uniqueneighborvals)==1
            clusteredmetric(vertexneighbors_thiscolor) = uniqueneighborvals;
            %if there are multiple cluster identifier values in the neighborhood, merge them into one
        else
            for valuenum = 2:length(uniqueneighborvals)
                clusteredmetric(clusteredmetric==uniqueneighborvals(valuenum)) = uniqueneighborvals(1);
            end
        end
    end
    uniqueclustervals = unique(clusteredmetric);
    uniqueclustervals(uniqueclustervals==0) = [];
    
    for clusternum = uniqueclustervals'
        if nnz(clusteredmetric==clusternum) < minsize
            neighborverts = unique(neighbors((clusteredmetric==clusternum),2:7));
            neighborverts(isnan(neighborverts)) = [];
            borderverts = setdiff(neighborverts,find(clusteredmetric==clusternum));
            borderverts(giftimap(borderverts)<1) = [];
            mode_neighborval = mode(giftimap(borderverts));
            giftimap(clusteredmetric==clusternum) = mode_neighborval;
        end
    end
end
variabilitymap(:,1) = giftimap;



distancemap = zeros(size(variabilitymap,1),length(networkIDs));
for IDnum = 1:length(networkIDs)
    networkID = networkIDs(IDnum);
    disp(['Network ' num2str(networkID)])
    networkindices = find(variabilitymap(:,1)==networkID);
    distancemap(networkindices,IDnum) = -1;
    
    outside_network_indices = find(sum((variabilitymap(:,2:end)==networkID),2));
    
    for vertex = outside_network_indices'
        
        dist = 0;
        found = 0;
        checked = vertex;
        while found==0
            dist = dist+1;
            neighs = unique(neighbors(checked,2:7));
            neighs(isnan(neighs)) = [];
            checked = [checked; neighs(:)];
            if any(variabilitymap(checked,1)==networkID)
                found = 1;
                distancemap(vertex,IDnum) = dist;
            end
        end
    end
end

save(gifti(single(distancemap)),[variabilitymapname(1:end-9) '_distancefromedge.func.gii'])
        
        