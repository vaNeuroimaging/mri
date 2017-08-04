function Removeverts_and_recluster_data(filename,removevertsfilename,combineIDs)
%Removeverts_and_recluster_data(filename,removevertsfilename,combineIDs)

label = gifti(filename); label = label.cdata;

if ischar(removevertsfilename)
    [a b c] = textread(removevertsfilename,'%s%s%s','delimiter',' ');
    charverts = c(strmatch('VERTEX',a));
    for i = 1:length(charverts)
        whiteloc = find(isspace(charverts{i}));
        if isempty(whiteloc)
            verts(i) = str2num(charverts{i})+1;
        else
            verts(i) = str2num(charverts{i}(1:(whiteloc(1)-1))) + 1;
        end
    end
else
    verts = removevertsfilename;
end

label(verts) = 0;


bufsize=16384;
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;



clusteredlabel = zeros(size(label));

data_inparcels = find(label);

for vertex = data_inparcels'
    
    %find the neighbors of this vertex
    vertexneighbors = neighbors(vertex,:);
    
    %find which of those neighbors also pass the thresholds
    vertexneighbors_inthresh = intersect(data_inparcels,vertexneighbors);
    
    %find if those neighbors have already been assigned different cluster values
    uniqueneighborvals = unique(clusteredlabel(vertexneighbors_inthresh));
    uniqueneighborvals(uniqueneighborvals==0) = [];
    
    %if no neighbors have cluster identifiers, assign them the number of this vertex as a unique cluster identifier 
    if isempty(uniqueneighborvals)
        clusteredlabel(vertexneighbors_inthresh) = vertex;
    %if there is only one previous cluster identifier present, make all the neighbors that value 
    elseif length(uniqueneighborvals)==1
        clusteredlabel(vertexneighbors_inthresh) = uniqueneighborvals;
    %if there are multiple cluster identifier values in the neighborhood, merge them into one 
    else
        for valuenum = 2:length(uniqueneighborvals)
            clusteredlabel(clusteredlabel==uniqueneighborvals(valuenum)) = uniqueneighborvals(1);
        end
    end
    
end

parcelIDs = unique(clusteredlabel); parcelIDs(parcelIDs==0) = [];
for i = 1:length(parcelIDs)
    if nnz(clusteredlabel==parcelIDs(i)) < 10;
        clusteredlabel(clusteredlabel==parcelIDs(i)) = 0;
    end
end





for combinenum = 1:size(combineIDs,1)
    
    target = min(combineIDs(combinenum,:));
    source = max(combineIDs(combinenum,:));
    
    borderverts = find(clusteredlabel==0);
    for vertnum = borderverts(:)';
        vertneighs = neighbors(vertnum,2:7); vertneighs(isnan(vertneighs)) = [];
        neighvals = unique(clusteredlabel(vertneighs))'; neighvals(neighvals==0) = [];
        if (length(neighvals) == 2) && (all(neighvals == [target source]));
            clusteredlabel(vertnum) = target;
        end
    end
    
    
    clusteredlabel(clusteredlabel==source) = target;
end
        
        



save(gifti(single(clusteredlabel)),[filename(1:end-9) '_tweaked.func.gii'])
