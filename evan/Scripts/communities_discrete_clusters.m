function outputmetric = communities_discrete_clusters(metricname)
%outputmetric = communities_discrete_clusters(metricname)

bufsize=16384;
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';

% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;

%load the gifti
if isnumeric(metricname)
    metric = metricname;
else
    metric = gifti(metricname);
    metric = metric.cdata;
end

outputmetric = [];

values = unique(metric); values(values<1) = [];
clustercount = 1;
for value = values(:)'
    
    
    %find which verticies meet the threshold criteria
    data_inthresh = find(metric==value);
    
    %initialize the metric keeping track of unique cluster identifiers
    clusteredmetric = zeros(size(metric));
    
    for vertex = data_inthresh'
        
        %find the neighbors of this vertex
        vertexneighbors = neighbors(vertex,:);
        
        %find which of those neighbors also pass the thresholds
        vertexneighbors_inthresh = intersect(data_inthresh,vertexneighbors);
        
        %find if those neighbors have already been assigned different cluster values
        uniqueneighborvals = unique(clusteredmetric(vertexneighbors_inthresh));
        uniqueneighborvals(uniqueneighborvals==0) = [];
        
        %if no neighbors have cluster identifiers, assign them the number of this vertex as a unique cluster identifier
        if isempty(uniqueneighborvals)
            clusteredmetric(vertexneighbors_inthresh) = vertex;
            %if there is only one previous cluster identifier present, make all the neighbors that value
        elseif length(uniqueneighborvals)==1
            clusteredmetric(vertexneighbors_inthresh) = uniqueneighborvals;
            %if there are multiple cluster identifier values in the neighborhood, merge them into one
        else
            for valuenum = 2:length(uniqueneighborvals)
                clusteredmetric(clusteredmetric==uniqueneighborvals(valuenum)) = uniqueneighborvals(1);
            end
        end
        
    end
    
    %find out what the unique cluster identifier values are
    uniqueclustervals = unique(clusteredmetric);
    uniqueclustervals(uniqueclustervals==0) = [];
    
    %save each unique cluster that passes the cluster size minimum into a column of the output metric
    
    for clusternum = 1:length(uniqueclustervals)
        outputmetric(logical(clusteredmetric == uniqueclustervals(clusternum)),clustercount) = value;
        clustercount = clustercount + 1;
    end
    
end