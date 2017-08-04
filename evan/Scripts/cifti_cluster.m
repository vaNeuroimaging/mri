function outputclusters = cifti_cluster(ciftiname,minthresh,maxthresh,minsize,cifti_template)
%outputclusters = cifti_cluster(ciftiname,minthresh,maxthresh,minsize,[cifti_template])

if exist('cifti_template')
    if ischar(cifti_template)
        neighbors = cifti_neighbors(cifti_template);
    else
        neighbors = cifti_template;
    end
else
%     if ischar(ciftiname)
%         neighbors = cifti_neighbors(ciftiname);
%     else
        neighbors = smartload('/home/data/scripts/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Cifti_surf_neighbors_LR_normalwall.mat');
%    end
end

if ischar(ciftiname)
    data = ft_read_cifti_mod(ciftiname); data = data.data;
else
    data = ciftiname;
end


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

%If there isn't a minimum cluster size defined, set it to zero. Otherwise, add the value to the output suffix 
if ~exist('minsize') || isempty(minsize)
    minsize=0;
end

outputclusters = zeros(size(clustereddata,1),0);

%save each unique cluster that passes the cluster size minimum into a column of the output metric 
clustercount = 1;

for clusternum = 1:length(uniqueclustervals)
    
    if nnz(clustereddata==uniqueclustervals(clusternum)) > minsize
        outputclusters(:,clustercount) = (clustereddata == uniqueclustervals(clusternum));
        clustercount = clustercount + 1;
    end
    
end
