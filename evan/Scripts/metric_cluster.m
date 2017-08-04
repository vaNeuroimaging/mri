function outputmetric = metric_cluster(metricname,minthresh,maxthresh,minsize,varargin)
%outputmetric = metric_cluster(metricname,minthresh,maxthresh,minsize,[outputstem])

bufsize=16384;
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';

% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;

%load the gifti
if islogical(metricname)
    metric = single(metricname);
elseif isnumeric(metricname)
    metric = metricname;
else
    metric = gifti(metricname);
    metric = metric.cdata;
end


%start building the output file suffix
suffix = 'clusters_';

%if min and max thresholds aren't specified, make them outside the range of values of the image
%Add to the suffix depending on the presence of these values
if isempty(minthresh)
    minthresh = min(metric)-1;
    suffix = [suffix 'below' num2str(maxthresh)];
elseif isempty(maxthresh)
    maxthresh = max(metric)+1;
    suffix = [suffix 'above' num2str(minthresh)];
else
    suffix = [suffix num2str(minthresh) '_to_' num2str(maxthresh)];
end

%find which verticies meet the threshold criteria
data_inthresh = find((metric >= minthresh) .* (metric <= maxthresh));

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

%If there isn't a minimum cluster size defined, set it to zero. Otherwise, add the value to the output suffix 
if isempty(minsize)
    minsize=0;
else
    suffix = [suffix '_minsize' num2str(minsize)];
end

%save each unique cluster that passes the cluster size minimum into a column of the output metric 
outputmetric = zeros(size(clusteredmetric,1),0);

clustercount = 1;
for clusternum = 1:length(uniqueclustervals)
    
    if length(find(clusteredmetric==uniqueclustervals(clusternum))) > minsize
        outputmetric(:,clustercount) = (clusteredmetric == uniqueclustervals(clusternum));
        clustercount = clustercount + 1;
    end
    
end

if ~isempty(varargin)
%save the output metric
save(gifti(single(outputmetric)),[varargin{1} '_' suffix '.func.gii'])
end