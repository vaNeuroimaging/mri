function parcel_creator_onlyedge(edgemetricname,outputdir,filestem,hem,threshperc)
%parcel_creator(edgemetricname,outputdir,filestem,hem,threshperc)

%% Set up parameters and load data

%minima will not be found in values higher than the (minimathreshperc)
%percentile of the edgemap
minimathreshperc = .75;

%smallest allowed parcel; smaller parcels will get merged with neighbors
minparcelsize = 20;

%smallest allowed parcel with no neighbors; smaller parcels will be deleted
minisolatedparcelsize = 20;

%vertices with edgemap values higher than the (edgevalthreshperc) will be
%removed from parcels
edgevalthreshperc = .75;

%Get node neighbors
global neighbors
bufsize=16384;
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;

%Load medial mask
medialmaskdata = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii']);
medialmaskdata = medialmaskdata.cdata;

%Load an eroded medial mask
erodedmedialmaskdata = gifti(['/data/cn4/laumannt/subcortical_mask/' hem '.atlasroi_erode3.32k_fs_LR.shape.gii']);
erodedmedialmaskdata = ~erodedmedialmaskdata.cdata;

%Get cortical indices
corticalindices = find(medialmaskdata==0); 

%Get medial indices
medialindices = find(medialmaskdata);

%Load edgemap
edgemetric = gifti(edgemetricname); edgemetric = edgemetric.cdata;

%Define various thresholds
sortedcorticalmetric = sort(edgemetric(corticalindices),'ascend');
minimathresh = sortedcorticalmetric(round(numel(sortedcorticalmetric)*minimathreshperc));
edgevalthresh = sortedcorticalmetric(round(numel(sortedcorticalmetric)*edgevalthreshperc));
mergethresh = sortedcorticalmetric(round(numel(sortedcorticalmetric)*threshperc));

%% Get all local minima


disp('Finding minima')

%make a copy of the edgemap to use
metric = edgemetric;

%nodes in eroded medial mask can't be minima
metric(erodedmedialmaskdata) = 100000;

%Get minima
minimametric = zeros(size(metric));
for i = 1:length(metric)
    
    %get this point's neighbors
    nodeneigh = neighbors(i,1:7);
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
%remove minima in an eroded medial mask
minimametric(logical(erodedmedialmaskdata)) = 0;

%remove all but one minimum in each basin of zeros
if nnz(edgemetric==0)>2
    zerobasins = metric_cluster(edgemetric,0,0,2);
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
hiter = unique(edgemetric(corticalindices));

%Iterate through the edgemap values
for i = 1:length(hiter);
    
    string{i} = ['Running watershed iteration ' num2str(i) ' out of ' num2str(length(hiter))];
    if i==1; fprintf('%s',string{i}); else fprintf([repmat('\b',1,length(string{i-1})) '%s'],string{i}); end
    
    % Take values in edgemap less than current iteration
    maskmetrics = edgemetric<hiter(i); 
    maskmetrics = maskmetrics & ~label>0 & ~watershedzone;
    
    maskpos = find(sum(maskmetrics,2)>0);
    maskpos = maskpos(randperm(length(maskpos)));
    
    for m = 1:length(maskpos) %For all nodes at this threshold
        
        %get node neighbors and labels of those neighbors
        nodeneigh = neighbors(maskpos(m),2:7);
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

%Remove values in medial wall
label(logical(medialmaskdata)) = 0;



%% Merge watersheds based on edgemap values between them

fprintf('%s',['Merging watersheds: .....'])

%Save a copy of the unaltered parcels
origlabel = label;

%Get the vertices that are watershed zones, i.e. borders between parcels
global borderindices
borderindices = intersect(find(label==0),corticalindices)';

%Get a list of the unique parcels
watersheds = unique(label);
watersheds(watersheds==0) = [];

%Find the vertices that are borders between each pair of parcels
[watershedborders adjacentwatersheds] = pairwise_border_verts(watersheds,label);

%Set up variable keeping track of which parcels still exist (haven't been merged)
watersheds_exist = ones(length(watersheds),1);

%Set up variable tracking the number of merges
merges = 0;

%Set up variable tracking whether we're finished merging
done = 0;

%Merge pairs of parcels iteratively, starting with the pair that has the
%smallest median edge map value between them, and going until the smallest
%median edgemap value is above the merge threshold
while done ==0
    medians = [];
    ijmat = [];
    
    %For each pair of parcels that still exist, get the median edgemap
    %value of the border between them
    for wateri = 1:length(watersheds)
        if logical(watersheds_exist(wateri))
            for waterj = 1:length(watersheds)
                if logical(watersheds_exist(waterj))
                    if wateri < waterj
                        if length(watershedborders{wateri,waterj}) > 1
                            medians(end+1,:) = median(edgemetric(watershedborders{wateri,waterj}));
                            ijmat(end+1,:) = [wateri waterj];
                        end
                    end
                end
            end
        end
    end
    
    %Find the smallest median
    [smallest_median smallestind] = min(medians);
    
    %If the smallest median is below the merge threshold, merge them
    if smallest_median < mergethresh
        merges = merges+1;
        string{merges} = [num2str(merges) ' watersheds merged'];
        if merges==1; fprintf('%s',string{merges}); else fprintf([repmat('\b',1,length(string{merges-1})) '%s'],string{merges}); end
        
        %Figure out which two parcels are being merged
        wateri = ijmat(smallestind,1); waterj = ijmat(smallestind,2);
        
        %Make the verts in one parcel the label value of the verts in the other
        label(label==mean(label(origlabel==watersheds(waterj)))) = mean(label(origlabel==watersheds(wateri)),1);
        
        %Keep track that one of the parcels doesn't exist anymore
        watersheds_exist(waterj) = 0;
        
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
borderindices = intersect(find(label==0),corticalindices)';

%Get a list of the unique parcels
watersheds = unique(label);
watersheds(watersheds==0) = [];

%Find the vertices that are borders between each pair of parcels
[watershedborders adjacentwatersheds] = pairwise_border_verts(watersheds,label);

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





%% Save final parcels

label = label==0;
label(medialindices) = 0;


save(gifti(single(label)),[outputdir '/' filestem 'thinedges.func.gii']);

numparcels = nnz(unique(label));
disp(['Final number of parcels: ' num2str(numparcels)])
end


function [watershedborders adjacentwatersheds]= pairwise_border_verts(watersheds,label)

global borderindices neighbors

for waternum = 1:length(watersheds)
    
    adjacentwatersheds{waternum} = [];
    %
    for waternum2 = 1:length(watersheds)
        watershedborders{waternum,waternum2} = [];
    end
    
end


for bordervertex = borderindices
    
    borderneighs = neighbors(bordervertex,2:7);
    borderneighs(isnan(borderneighs)) = [];
    borderneighvals = label(borderneighs);
    borderneighvals(borderneighvals==0) = [];
    
    watershedneighbors = unique(borderneighvals);
    
    for waterneighbornum = 1:length(watershedneighbors)
        
        thiswatershedindex(waterneighbornum) = find(watersheds==watershedneighbors(waterneighbornum));
        
        otherwatershedneighbors = watershedneighbors;
        otherwatershedneighbors(otherwatershedneighbors==watershedneighbors(waterneighbornum)) = [];
        
        %watershedborders{thiswatershedindex(waterneighbornum)} = unique([watershedborders{thiswatershedindex(waterneighbornum)} bordervertex]);
        adjacentwatersheds{thiswatershedindex(waterneighbornum)} = unique([adjacentwatersheds{thiswatershedindex(waterneighbornum)} otherwatershedneighbors']);
    end
    
    if length(watershedneighbors) == 2 && length(borderneighvals)>2
        
        watershedborders{thiswatershedindex(1),thiswatershedindex(2)} = unique([watershedborders{thiswatershedindex(1),thiswatershedindex(2)} bordervertex]);
        watershedborders{thiswatershedindex(2),thiswatershedindex(1)} = unique([watershedborders{thiswatershedindex(2),thiswatershedindex(1)} bordervertex]);
        
    end
    
end

end

function label = remove_merged_borders(label)
global borderindices neighbors
for bordervertex = borderindices
        
        borderneighs = neighbors(bordervertex,2:7);
        borderneighs(isnan(borderneighs)) = [];
        borderneighvals = label(borderneighs);
        borderneighvals(borderneighvals==0) = [];
        if length(unique(borderneighvals)) == 1;
            label(bordervertex) = unique(borderneighvals);
        end
end
end

function label = discrete_clusters(label)
global borderindices neighbors

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



