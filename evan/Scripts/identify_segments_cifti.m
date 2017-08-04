function segments = identify_segments_cifti(ciftifile,neighbors,distances)
%segments = identify_segments_cifti(ciftifile,[neighbors],[distances])


datastruct = ft_read_cifti_mod(ciftifile);
data = datastruct.data;
if ~exist('neighbors') || isempty(neighbors)
    neighbors = cifti_neighbors(ciftifile);
end

        
evalc(['!wb_command -cifti-smoothing ' ciftifile ' 5 0.01 COLUMN Tempsmooth.dtseries.nii -left-surface /data/cn/data1/scripts/CIFTI_RELATED/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.L.midthickness.32k_fs_LR.surf.gii -right-surface /data/cn/data1/scripts/CIFTI_RELATED/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.R.midthickness.32k_fs_LR.surf.gii']);
datasmooth = ft_read_cifti_mod('Tempsmooth.dtseries.nii'); datasmooth = datasmooth.data;
%delete('Tempsmooth.dtseries.nii')

datasmooth = datasmooth + data;

%Get minima
minima = zeros(size(datasmooth));


temp_datasmooth = datasmooth;
for i = 1:length(datasmooth)
    
    %get this point's neighbors
    nodeneigh = neighbors(i,1:end);
    nodeneigh(isnan(nodeneigh)) = [];
    
    %get the minimum edge values of this point and its neighbors
    [minval mini] = min(temp_datasmooth(nodeneigh));
    minindices = find(temp_datasmooth(nodeneigh)==minval);
    
    %if this point is the minimum
    if minval == temp_datasmooth(i)
        %add this point to the minima metric being built
        minima(i) = 1;
        %if there were multiple clustered minima (a "basin")
        if length(minindices) > 1
            %make the rest of them not minima
            minindices(logical(minindices==1)) = [];
            temp_datasmooth(nodeneigh(minindices)) = minval+.00001;
            
        end
    end
    
end

%remove all but one minimum in each basin of zeros
if any(datasmooth==0)
    zerobasins = cifti_cluster(datasmooth,0,0,2,neighbors);
    for basinnum = 1:size(zerobasins,2)
        basininds = find(zerobasins(:,basinnum));
        minima(basininds) = 0;
        minima(basininds(floor(length(basininds)/2))) = 1;
    end
end

minima(datasmooth==2) = 0;

%delete('Tempsmooth.dtseries.nii');

 


%Initialize the parcels and put unique values at all local minima
label = zeros(size(minima));
labelnum = nnz(minima);
labelpos = find(minima==1);
[ign sortorder] = sort(datasmooth(labelpos));
for j = 1:labelnum;
    label(labelpos(j)) = sortorder(j);
end

%Initialize a variable keeping track of final border positions
watershedzone = zeros(size(label));

%Find the unique edgemap values, which are the iterations used for
%watershed parcel growing
hiter = unique(datasmooth);

%Iterate through the edgemap values
for i = 1:length(hiter);
    
    %string{i} = ['Growing parcels through ' num2str(i) ' out of ' num2str(length(hiter)) ' values'];
    %if i==1; fprintf('%s',string{i}); else fprintf([repmat('\b',1,length(string{i-1})) '%s'],string{i}); end
    
    % Take values in edgemap less than current iteration
    maskmetrics = datasmooth<hiter(i); 
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

%disp(' ')


%fprintf('%s',['Merging parcels: .....'])

%Save a copy of the unaltered parcels
origlabel = label;

%Get the vertices that are watershed zones, i.e. borders between parcels
%global borderindices
borderindices = find(label==0)';

%Get a list of the unique parcels
watersheds = unique(label);
watersheds(watersheds==0) = [];

%Find the vertices that are borders between each pair of parcels
%[watershedborders, adjacentwatersheds, watershedshaveborders] = pairwise_border_verts(watersheds,label);


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
    %minvals = [];
    medians = [];
    ijmat = [];
    for i = 1:length(wateri_inds);
        wateri = wateri_inds(i);
        waterj = waterj_inds(i);
        if watersheds_exist(wateri) && watersheds_exist(waterj) && (length(watershedborders{wateri,waterj}) > 1)
            %minvals(end+1,1) = min(data(watershedborders{wateri,waterj}));
            medians(end+1,1) = median(data(watershedborders{wateri,waterj}));
            ijmat(end+1,:) = [wateri waterj];
        end
    end
    
    
    %Find the smallest median border value
    [smallest_min, smallestind] = min(medians);%minvals);
    
    %If the smallest median is below the merge threshold, merge them
    if smallest_min == 0
        merges = merges+1;
        %string{merges} = [num2str(merges) ' parcels merged'];
        %if merges==1; fprintf('%s',string{merges}); else fprintf([repmat('\b',1,length(string{merges-1})) '%s'],string{merges}); end
        
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
borderindices_thatarezeros = intersect(borderindices,find(data==0));
for bordervertex = borderindices_thatarezeros(:)'
        
        borderneighs = neighbors(bordervertex,2:end);
        borderneighs(isnan(borderneighs)) = [];
        borderneighvals = label(borderneighs);
        borderneighvals(borderneighvals==0) = [];
        if length(unique(borderneighvals)) == 1
            label(bordervertex) = unique(borderneighvals);
        end
end

clear watershedborders adjacentwatersheds edgeval_withinwater
 
%disp(' ')

borders_tolabel = find(data .* (label==0));
neighborlabel_tracker = zeros(length(borders_tolabel),2);
labeledsegments = zeros(size(label));

for indnum = 1:length(borders_tolabel)
    ind = borders_tolabel(indnum);
    vertneighs = neighbors(ind,2:end); vertneighs(isnan(vertneighs)) = [];
    neighborlabels = unique(label(vertneighs)); neighborlabels(neighborlabels==0) = [];
    
    if (~isempty(neighborlabels)) && (nnz(neighborlabels) <= 2)
        neighborlabel_tracker(indnum,:) = neighborlabels(:)';
    else
        neighborlabel_tracker(indnum,:) = NaN;
    end
end

uniqueneighborlabels = unique(neighborlabel_tracker,'rows'); 
notnan_rows = ~any(isnan(uniqueneighborlabels),2);
uniqueneighborlabels = uniqueneighborlabels(notnan_rows,:);

skeleton_segments = zeros(size(data,1),1);

for seglabelnum = 1:size(uniqueneighborlabels,1)
    matchingvertinds = all(neighborlabel_tracker==repmat(uniqueneighborlabels(seglabelnum,:),size(neighborlabel_tracker,1),1),2);
    labeledsegments(borders_tolabel(matchingvertinds)) = seglabelnum;
    
    clustered = cifti_cluster(labeledsegments,seglabelnum,seglabelnum,3,neighbors);
    
    for clusnum = 1:size(clustered,2)
        thislabel = max(skeleton_segments(:)) + 1;
        skeleton_segments = skeleton_segments + (clustered(:,clusnum) .* thislabel);
    end
            
end

if ~exist('distances') 
    distances = smartload('/data/cn/data1/scripts/CIFTI_RELATED/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Cifti_geo_distances_uint8.mat');
end

skeletonized_inds = find(skeleton_segments);

datainds = find(data);

[~, mindist_inds] = min(distances(datainds,skeletonized_inds),[],2);

segments = zeros(size(data));
segments(datainds) = skeleton_segments(skeletonized_inds(mindist_inds));
    
        
    











