function parcel_creator_cifti(edgeciftiname,filestem,threshperc)
%parcel_creator_cifti(edgeciftiname,filestem,threshperc)
%
%generates discrete parcels from a boundary map cifti file created by the
%surface_parcellation function.
%
% 'edgeciftiname' is the name of the boundary map
%
% 'filestem' is the desired output name of the parcels
%
% 'threshperc' specifies the threshold for merging neighboring parcels.
% This is a value representing a percentile value of the edgemap. If the
% median values of borders between parcels are below this threshold, they
% will be merged together. Testing suggests that percentiles between ~.35
% and .50 work best, but check how well the resulting parcels fit your
% edgemap in order to determine an optimal value for your parcellation.
%
% Note that several more parameters are specified at the top of the script
% that govern the width of borders between parcels and the minimum sizes
% of parcels. These generally should not be changed unless there is a good
% reason.
%
%EMG 06/24/15

%% Set up parameters

%minima will not be found in values higher than the (minimathreshperc)
%percentile of the edgemap
minimathreshperc = .75;

%vertices with edgemap values higher than the (edgevalthreshperc)
%percentile will be removed from parcels
edgevalthreshperc = .75;

%smallest allowed parcel; smaller parcels will get merged with neighbors
minparcelsize = 30;

%smallest allowed parcel with no neighbors; smaller parcels will be deleted
minisolatedparcelsize = 10;


%% Load data

%Get node neighbors
global neighbors
neighbors = cifti_neighbors(edgeciftiname);


%Load edgemap
cifti_template = ft_read_cifti_mod(edgeciftiname);
edgemetric = cifti_template.data;

%Define various thresholds
sortedcorticalmetric = sort(edgemetric,'ascend');
minimathresh = sortedcorticalmetric(round(numel(sortedcorticalmetric)*minimathreshperc));
edgevalthresh = sortedcorticalmetric(round(numel(sortedcorticalmetric)*edgevalthreshperc));
mergethresh = sortedcorticalmetric(round(numel(sortedcorticalmetric)*threshperc));

%% Get all local minima


edgemetric= edgemetric - min(edgemetric);

disp('Finding minima')

%make a copy of the edgemap to use
metric = edgemetric;


%Get minima
minimametric = zeros(size(metric));
for i = 1:length(metric)
    
    %get this point's neighbors
    nodeneigh = neighbors(i,1:end);
    nodeneigh(isnan(nodeneigh)) = [];
    
    %get the minimum edge values of this point and its neighbors
    [minval mini] = min(metric(nodeneigh));
    minindices = find(metric(nodeneigh)==minval);
    
    %if this point is the minimum
    if minval == metric(i)
        %add this point to the minima metric being built
        minimametric(i) = 1;
        %if there were multiple clustered minima (a "basin")
        if length(minindices) > 1
            %make the rest of them not minima
            minindices(logical(minindices==1)) = [];
            metric(nodeneigh(minindices)) = minval+.00001;
            
        end
    end
    
end

%remove minima with high edge values
minimametric(edgemetric>minimathresh) = 0;

%remove all but one minimum in each basin of zeros
if any(edgemetric==0)
    zerobasins = cifti_cluster(edgemetric,0,0,2);
    for basinnum = 1:size(zerobasins,2)
        basininds = find(zerobasins(:,basinnum));
        minimametric(basininds) = 0;
        minimametric(basininds(1)) = 1;
    end
end

clear metric

%% Grow parcels from minima using watershed technique

%Remove NaNs from edgemap (probably not any)
edgemetric(isnan(edgemetric)) = 0;

%Initialize the parcels and put unique values at all local minima
label = zeros(size(minimametric));
labelnum = nnz(minimametric);
labelpos = find(minimametric==1);
[ign sortorder] = sort(edgemetric(labelpos));
for j = 1:labelnum;
    label(labelpos(j)) = sortorder(j);
end

%Initialize a variable keeping track of final border positions
watershedzone = zeros(size(label));

%Find the unique edgemap values, which are the iterations used for
%watershed parcel growing
hiter = unique(edgemetric);

%Iterate through the edgemap values
for i = 1:length(hiter);
    
    string{i} = ['Growing parcels through ' num2str(i) ' out of ' num2str(length(hiter)) ' values'];
    if i==1; fprintf('%s',string{i}); else fprintf([repmat('\b',1,length(string{i-1})) '%s'],string{i}); end
    
    % Take values in edgemap less than current iteration
    maskmetrics = edgemetric<hiter(i); 
    maskmetrics = maskmetrics & ~label>0 & ~watershedzone;
    
    maskpos = find(sum(maskmetrics,2)>0);
    maskpos = maskpos(randperm(length(maskpos)));
    
    for m = 1:length(maskpos) %For all nodes at this threshold
        
        %get node neighbors and labels of those neighbors
        nodeneigh = neighbors(maskpos(m),2:end);
        nodeneigh(isnan(nodeneigh)) = [];
        nodeneighlab = label(nodeneigh);
        
        %Find minimum value other than 0 among neighbors
        minfindnodeneighlab = nodeneighlab;
        minfindnodeneighlab(nodeneighlab==0) = 100000; 
        minnodeneighlab = min(minfindnodeneighlab,[],1);
        
        %Find maximum value other than 0 among neighbors
        maxfindnodeneighlab = nodeneighlab;
        maxfindnodeneighlab(nodeneighlab==0) = -100000;
        maxnodeneighlab = max(maxfindnodeneighlab,[],1);
       
        %If min and max differ (i.e. two or more neighbor parcels), it's a
        %border
        maskinthismetric = maskmetrics(maskpos(m),:);
        watershed_nodes = (minnodeneighlab~=maxnodeneighlab) & (minnodeneighlab~=100000) & maskinthismetric;
        label(maskpos(m),watershed_nodes) = 0;
        watershedzone(maskpos(m),watershed_nodes) = 1;
        
        %If min and max are the same but different from 0, make the node
        %that value
        next_to_water = (minnodeneighlab==maxnodeneighlab) & (minnodeneighlab~=100000) & maskinthismetric;
        label(maskpos(m),next_to_water) = minnodeneighlab(next_to_water);
       
    end
end

disp(' ')



%% Merge watersheds based on edgemap values between them

fprintf('%s',['Merging parcels: .....'])

%Save a copy of the unaltered parcels
origlabel = label;

%Get the vertices that are watershed zones, i.e. borders between parcels
global borderindices
borderindices = find(label==0)';

%Get a list of the unique parcels
watersheds = unique(label);
watersheds(watersheds==0) = [];

%Find the vertices that are borders between each pair of parcels
[watershedborders, adjacentwatersheds, watershedshaveborders] = pairwise_border_verts(watersheds,label);

%Set up variable keeping track of which parcels still exist (haven't been merged)
watersheds_exist = true(length(watersheds),1);

%Set up variable tracking the number of merges
merges = 0;

%Set up variable tracking whether we're finished merging
done = 0;

%Merge pairs of parcels iteratively, starting with the pair that has the
%smallest median edge map value between them, and going until the smallest
%median edgemap value is above the merge threshold
while done ==0

    %Find which parcel pairs have borders
    [wateri_inds, waterj_inds] = find(triu(watershedshaveborders,1));
    
    %For each pair of parcels that still exist, get the median edgemap
    %value of the border between them
    medians = [];
    ijmat = [];
    for i = 1:length(wateri_inds);
        wateri = wateri_inds(i);
        waterj = waterj_inds(i);
        if watersheds_exist(wateri) && watersheds_exist(waterj) && (length(watershedborders{wateri,waterj}) > 1)
            medians(end+1,1) = median(edgemetric(watershedborders{wateri,waterj}));
            ijmat(end+1,:) = [wateri waterj];
        end
    end
    
    
    %Find the smallest median border value
    [smallest_median, smallestind] = min(medians);
    
    %If the smallest median is below the merge threshold, merge them
    if smallest_median < mergethresh
        merges = merges+1;
        string{merges} = [num2str(merges) ' parcels merged'];
        if merges==1; fprintf('%s',string{merges}); else fprintf([repmat('\b',1,length(string{merges-1})) '%s'],string{merges}); end
        
        %Figure out which two parcels are being merged
        wateri = ijmat(smallestind,1); waterj = ijmat(smallestind,2);
        
        %Make the verts in one parcel the label value of the verts in the other
        label(label==mean(label(origlabel==watersheds(waterj)))) = mean(label(origlabel==watersheds(wateri)),1);
        
        %Keep track that one of the parcels doesn't exist anymore
        watersheds_exist(waterj) = false;
        
        %Appropriately update the existence of borders
        watershedshaveborders(wateri,:) = (watershedshaveborders(wateri,:) | watershedshaveborders(waterj,:));
        watershedshaveborders(:,wateri) = (watershedshaveborders(:,wateri) | watershedshaveborders(:,waterj));
        watershedshaveborders(waterj,:) = false;
        watershedshaveborders(:,waterj) = false;
        
        %Appropriately update the borders
        for k = 1:length(watersheds)
            watershedborders{wateri,k} = unique([watershedborders{wateri,k} watershedborders{waterj,k} watershedborders{k,wateri} watershedborders{k,waterj}]);
            watershedborders{k,wateri} = watershedborders{wateri,k};
        end
    else
        done = 1;
    end
end

%Remove borders that are now between merged parcels
label = remove_merged_borders(label);

clear watershedborders adjacentwatersheds edgeval_withinwater
 
disp(' ')

%% Merge small neighboring parcels together

disp('Merging small parcels')

%Get the vertices that are watershed zones, i.e. borders between parcels
borderindices = find(label==0)';

%Get a list of the unique parcels
watersheds = unique(label);
watersheds(watersheds==0) = [];

%Find the vertices that are borders between each pair of parcels
[watershedborders, adjacentwatersheds, watershedshaveborders] = pairwise_border_verts(watersheds,label);

%For each parcel
for waternum = 1:length(watersheds)
    
    %If the parcel exists but is smaller than the minimum parcel size
    if (nnz(label==watersheds(waternum)) > 0) && (nnz(label==watersheds(waternum)) < minparcelsize);
        
        %Initialize variable keeping track of which neighboring parcel has
        %the smallest median border
        minwaternum = 0;
        
        %Initialize variable keeping track of the smallest median border
        %among nighboring parcels
        minedge = 1;
        
        %For each neighboring parcel
        for i = 1:length(adjacentwatersheds{waternum})
            thisadjacentwatershedindex = find(watersheds == adjacentwatersheds{waternum}(i));
            
            %If the median edgemap value is below the previous minimum
            if median(edgemetric(watershedborders{waternum,thisadjacentwatershedindex}),1) < minedge
                %Then this is the new minimum
                minwaternum = thisadjacentwatershedindex;
                minedge = mean(edgemetric(watershedborders{waternum,thisadjacentwatershedindex}),1);
            end
        end
        
        %If a minimum was found (i.e. there are neighbors)
        if minwaternum
            %Then merge the parcelw with the neighboring parcel that had
            %the minimum median border
            label(label==watersheds(waternum)) = watersheds(minwaternum);
        end
    end
end

%Remove borders that are now between merged parcels
label = remove_merged_borders(label);


clear adjacentwatersheds edgeval_withinwater watershedborders thiswatershedindex adjacentwatersheds watershedborders

%% Remove high edge values and merge new small parcels

disp('Eliminating high edge value vertices and eliminating new small parcels')

%Remove vertices with edgemap values above the threshold from the parcels
label(edgemetric>edgevalthresh) = 0;

%Make each discrete, contiguous cluster of label values a new parcel
label = discrete_clusters(label);

%Get the vertices that are watershed zones, i.e. borders between parcels
borderindices = find(label==0)';

%Get a list of the unique parcels
watersheds = unique(label);
watersheds(watersheds==0) = [];

%Find the vertices that are borders between each pair of parcels
[watershedborders, adjacentwatersheds, watershedshaveborders] = pairwise_border_verts(watersheds,label);

%Save a copy of the unaltered parcels
origlabel = label;

%For each parcel
for waternum = 1:length(watersheds)
    
    %If the parcel exists but is smaller than the minimum parcel size
    if (nnz(label==watersheds(waternum)) > 0) && (nnz(label==watersheds(waternum)) < minparcelsize);
        
        %Initialize variable keeping track of which neighboring parcel has
        %the smallest median border
        minwaternum = 0;
        
        %Initialize variable keeping track of the smallest median border
        %among nighboring parcels
        minedge = 1;
        
        %For each neighboring parcel
        for i = 1:length(adjacentwatersheds{waternum})
            thisadjacentwatershedindex = find(watersheds == adjacentwatersheds{waternum}(i));
            
            %If the median edgemap value is below the previous minimum
            if ~isempty(watershedborders{waternum,thisadjacentwatershedindex}) && (median(edgemetric(watershedborders{waternum,thisadjacentwatershedindex}),1) < minedge)% && (nnz(edgemetric(watershedborders{waternum,thisadjacentwatershedindex}) < thresh) > 2)
                %Then this is the new minimum
                minwaternum = thisadjacentwatershedindex;
                minedge = mean(edgemetric(watershedborders{waternum,thisadjacentwatershedindex}),1);
            end
        end
        
        %If a minimum was found (i.e. there are neighbors)
        if minwaternum
            %Then merge the parcelw with the neighboring parcel that had
            %the minimum median border
            label(label==watersheds(waternum)) = mean(label(origlabel==watersheds(minwaternum)));
            
            %And appropriately update the borders
            label(watershedborders{waternum,minwaternum}) = mean(label(origlabel==watersheds(minwaternum)));
        end
    end
end

%Remove borders that are now between merged parcels
label = remove_merged_borders(label);

%Remove vertices with edgemap values above the threshold from the parcels one last time
label(edgemetric>edgevalthresh) = 0;

%Split parcels that are joined only by a single vertex
verts_inparcels = find(label);
for vert = verts_inparcels(:)'
    vertneighs = neighbors(vert,2:end); vertneighs(isnan(vertneighs)) = [];
    nonzerovertneighs = vertneighs(label(vertneighs)>0);
    if length(nonzerovertneighs)==2
        neighofneigh = neighbors(nonzerovertneighs(1),2:end); neighofneigh(isnan(neighofneigh))=[];
        if isempty(intersect(neighofneigh,nonzerovertneighs(2)))
            label(vert) = 0;
        end
    end
end

%Make each discrete, contiguous cluster of label values a new parcel
label = discrete_clusters(label);

%% Delete isolated parcels smaller than a minimum size


%Get a list of the unique parcels
watersheds = unique(label); 
watersheds(watersheds==0) = [];

%Delete parcels smaller than the minimum size
for watershed = watersheds'
    if nnz(label==watershed) < minisolatedparcelsize
        label(label==watershed) = 0;
    end
end


%% Save final parcels

cifti_template.data = label;
ft_write_cifti_mod([filestem '_edgethresh_' num2str(round(threshperc*100)/100)],cifti_template);

numparcels = nnz(unique(label));
disp(['Final number of parcels: ' num2str(numparcels)])
end


function [watershedborders, adjacentwatersheds, watershedshaveborders]= pairwise_border_verts(watersheds,label)
%find which parcels are adjacent and which vertices are on the borders
%between those parcels
global borderindices neighbors

adjacentwatersheds = cell(length(watersheds),1);
watershedborders = cell(length(watersheds));
watershedshaveborders = false(length(watersheds));


for bordervertex = borderindices
    
    borderneighs = neighbors(bordervertex,2:end);
    borderneighs(isnan(borderneighs)) = [];
    borderneighvals = label(borderneighs);
    borderneighvals(borderneighvals==0) = [];
    
    watershedneighbors = unique(borderneighvals);
    
    for waterneighbornum = 1:length(watershedneighbors)
        
        thiswatershedindex(waterneighbornum) = find(watersheds==watershedneighbors(waterneighbornum));
        
        otherwatershedneighbors = watershedneighbors;
        otherwatershedneighbors(otherwatershedneighbors==watershedneighbors(waterneighbornum)) = [];
        
        adjacentwatersheds{thiswatershedindex(waterneighbornum)} = unique([adjacentwatersheds{thiswatershedindex(waterneighbornum)} otherwatershedneighbors']);
    end
    
    if length(watershedneighbors) == 2 && length(borderneighvals)>2
        
        watershedborders{thiswatershedindex(1),thiswatershedindex(2)} = unique([watershedborders{thiswatershedindex(1),thiswatershedindex(2)} bordervertex]);
        watershedborders{thiswatershedindex(2),thiswatershedindex(1)} = unique([watershedborders{thiswatershedindex(2),thiswatershedindex(1)} bordervertex]);
        
        watershedshaveborders(thiswatershedindex(1),thiswatershedindex(2)) = true;
        watershedshaveborders(thiswatershedindex(2),thiswatershedindex(1)) = true;
        
    end
    
end
watershedshaveborders(diag(true(length(watershedshaveborders)))) = false;
end

function label = remove_merged_borders(label)
%delete the border vertices between merged parcels
global borderindices neighbors
for bordervertex = borderindices
        
        borderneighs = neighbors(bordervertex,2:end);
        borderneighs(isnan(borderneighs)) = [];
        borderneighvals = label(borderneighs);
        borderneighvals(borderneighvals==0) = [];
        if length(unique(borderneighvals)) == 1;
            label(bordervertex) = unique(borderneighvals);
        end
end
end

function label = discrete_clusters(label)
%define new parcels that may have been created after removing high edgemap
%values
global neighbors

%initialize the metric keeping track of unique cluster identifiers
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

label = clusteredlabel;
end

function outputmetric = cifti_cluster(metric,minthresh,maxthresh,minsize)
%find discrete clusters of vertices meeting a specific criterion

global neighbors


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
clustercount = 1;
for clusternum = 1:length(uniqueclustervals)
    
    if length(find(clusteredmetric==uniqueclustervals(clusternum))) > minsize
        outputmetric(:,clustercount) = (clusteredmetric == uniqueclustervals(clusternum));
        clustercount = clustercount + 1;
    end
    
end

end



