function clustereddata = cifti_cluster_surfacearea_volume(data,minthresh,maxthresh,minsize_mm2,minsize_mm3,ncortverts,surfacearea_voxvol,neighbors)
%outputclusters = cifti_cluster_surfacearea_volume(data,minthresh,maxthresh,minsize_mm2,minsize_mm3,ncortverts,surfacearea_voxvol,neighbors)



%find which verticies meet the threshold criteria
data_inthresh = find((data >= minthresh) .* (data <= maxthresh));

%initialize the metric keeping track of unique cluster identifiers
clustereddata = zeros(size(data));

for vertex = data_inthresh'
    
    %find the neighbors of this vertex
    vertexneighbors = neighbors(vertex,:);
    
    %find which of those neighbors also pass the thresholds
    vertexneighbors_inthresh = intersect(data_inthresh,vertexneighbors);
    
    %find if those neighbors have already been assigned different cluster values
    uniqueneighborvals = unique(clustereddata(vertexneighbors_inthresh));
    uniqueneighborvals(uniqueneighborvals==0) = [];
    
    %if no neighbors have cluster identifiers, assign them the number of this vertex as a unique cluster identifier 
    if isempty(uniqueneighborvals)
        clustereddata(vertexneighbors_inthresh) = vertex;
    %if there is only one previous cluster identifier present, make all the neighbors that value 
    elseif length(uniqueneighborvals)==1
        clustereddata(vertexneighbors_inthresh) = uniqueneighborvals;
    %if there are multiple cluster identifier values in the neighborhood, merge them into one 
    else
        for valuenum = 2:length(uniqueneighborvals)
            clustereddata(clustereddata==uniqueneighborvals(valuenum)) = uniqueneighborvals(1);
        end
    end
    
end

%find out what the unique cluster identifier values are
uniqueclustervals = unique(clustereddata);
uniqueclustervals(uniqueclustervals==0) = [];

clustervalorder = [];
clustersizes = [];
for clusternum = uniqueclustervals(:)'
    
    if any(find(clustereddata==clusternum) < ncortverts)
        if sum(surfacearea_voxvol(clustereddata==clusternum)) < minsize_mm2
            clustereddata(clustereddata==clusternum) = 0;
        else
            clustervalorder(end+1) = clusternum;
            clustersizes(end+1) = sum(surfacearea_voxvol(clustereddata==clusternum));
        end
    end
end
for clusternum = uniqueclustervals(:)'
    if ~any(find(clustereddata==clusternum) < ncortverts)
        if sum(surfacearea_voxvol(clustereddata==clusternum)) < minsize_mm3
            clustereddata(clustereddata==clusternum) = 0;
            else
            clustervalorder(end+1) = clusternum;
            clustersizes(end+1) = sum(surfacearea_voxvol(clustereddata==clusternum));
        end
    end
end

clustereddata_temp = zeros(size(clustereddata));
[~,sortedinds] = sort(clustersizes,'descend');
for i = 1:length(sortedinds)
    clustereddata_temp(clustereddata==clustervalorder(sortedinds(i))) = i;
end
clustereddata = clustereddata_temp;


    
