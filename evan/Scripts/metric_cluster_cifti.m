function outputcifti = metric_cluster_cifti(cifti,minthresh,maxthresh,minsize)
%outputcifti = metric_cluster_cifti(cifti,minthresh,maxthresh,minsize)

bufsize=16384;
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';

% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;

medial_wall{1} = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii']);
medial_wall{1} = medial_wall{1}.cdata;
medial_wall{2} = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.R.32k_fs_LR.func.gii']);
medial_wall{2} = medial_wall{2}.cdata;



        ciftiinds{1} = 1:nnz(medial_wall{1}==0);
        ciftiinds{2} = (nnz(medial_wall{1}==0)+1) : (nnz(medial_wall{1}==0) + nnz(medial_wall{2}==0));


metric_hem{2} = zeros(size(medial_wall{2}));
metric_hem{2}(medial_wall{2}==0) = cifti(ciftiinds{2});

outputcifti = zeros(size(cifti,1),0);

for hem = 1:2
  
metric = zeros(size(medial_wall{hem}));
metric(medial_wall{hem}==0) = cifti(ciftiinds{hem});

%load the gifti



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
end

%save each unique cluster that passes the cluster size minimum into a column of the output metric 
outputmetric = zeros(size(clusteredmetric,1),1);
clustercount = 1;
for clusternum = 1:length(uniqueclustervals)
    
    if length(find(clusteredmetric==uniqueclustervals(clusternum))) > minsize
        outputmetric(:,clustercount) = (clusteredmetric == uniqueclustervals(clusternum));
        clustercount = clustercount + 1;
    end
    
end
    if any(outputmetric)
        outputcifti(ciftiinds{hem}, (end+1) : (end+size(outputmetric,2))) = outputmetric(medial_wall{hem}==0,:);
    end
    
    clear outputmetric

end