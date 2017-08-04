function Fill_Craddock_parcel_holes(parcelname,outname,hem)

% parcelname = 'tcorr05_2level_all0030_R.func.gii';
% outname = 'Craddock_350_R.func.gii';
% hem = 'R';



bufsize=16384;
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;

maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii'];
mask = gifti(maskname);
mask = mask.cdata;
medialinds = find(mask);

parcels = gifti(parcelname); parcels = parcels.cdata;
parcels(medialinds) = 0;

parcelIDs = unique(parcels); parcelIDs(parcelIDs==0) = [];

for i = 1:length(parcelIDs)
    %find which verticies meet the threshold criteria
data_inthresh = find(parcels==parcelIDs(i));

%initialize the metric keeping track of unique cluster identifiers
clusteredmetric = zeros(size(parcels));

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

clustersize = zeros(length(uniqueclustervals),1);
for clustervalnum = 1:length(uniqueclustervals)
    clustersize(clustervalnum) = nnz(clusteredmetric==uniqueclustervals(clustervalnum));
end
for clustervalnum = 1:length(uniqueclustervals)
    if nnz(clusteredmetric==uniqueclustervals(clustervalnum)) ~= max(clustersize);
        parcels(clusteredmetric==uniqueclustervals(clustervalnum)) = 0;
    end
end



end








zeroinds = find(parcels==0);
zeroinds = setdiff(zeroinds,medialinds);

while ~isempty(zeroinds)
tempparcels = parcels;
for vert = zeroinds(:)'
   vertneighs = neighbors(vert,2:7); vertneighs(isnan(vertneighs)) = [];
   neighvals = parcels(vertneighs); neighvals(neighvals==0) = [];
   if ~isempty(neighvals)
       tempparcels(vert) = mode(neighvals);
   end
end
parcels = tempparcels;
zeroinds = find(parcels==0);
zeroinds = setdiff(zeroinds,medialinds);
end

save(gifti(single(parcels)),outname)

