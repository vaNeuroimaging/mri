function outputclusters = cifti_cluster_surface(ciftiname,minthresh,maxthresh,minsize,surfaceareas,cifti_template)
%outputclusters = cifti_cluster_surface(ciftiname,minthresh,maxthresh,minsize,surfaceareas,cifti_template)

if exist('cifti_template')
    if ischar(cifti_template)
        neighbors = cifti_neighbors(cifti_template);
    else
        neighbors = cifti_template;
    end
else
    if ischar(ciftiname)
        neighbors = cifti_neighbors(ciftiname);
    else
        neighbors = smartload('/data/cn4/evan/fsaverage_LR32k/Cifti_surf_neighbors_LR_normalwall.mat');
    end
end

if ~exist('surfaceareas')
    surfaceareas = cifti_read('/data/cn4/evan/fsaverage_LR32k/Conte69.LR.midthickness.32k_fs_LR_surfaceareas_normalwall.dtseries.nii');
end

if ischar(ciftiname)
    data = cifti_read(ciftiname);
else
    data = ciftiname;
end

clustercount = 0;

%find which verticies meet the threshold criteria
data_inthresh = find((data >= minthresh) .* (data <= maxthresh));

if ~isempty(data_inthresh)
    
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
    
    %save each unique cluster that passes the cluster size minimum into a column of the output metric
    for clusternum = 1:length(uniqueclustervals)
        
        if sum(surfaceareas(clustereddata==uniqueclustervals(clusternum))) > minsize
            clustercount = clustercount + 1;
            outputclusters(:,clustercount) = (clustereddata == uniqueclustervals(clusternum));
        end
        
    end
end
if clustercount == 0
    outputclusters = zeros(size(data));
end
