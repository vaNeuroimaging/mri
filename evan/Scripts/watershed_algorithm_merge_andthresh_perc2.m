function watershed_algorithm_merge_andthresh_perc2(edgemetricname,outputdir,filestem,hem,threshperc)
%watershed_algorithm_merge_andthresh_perc(edgemetricname,outputdir,filestem,hem,threshperc)

% edgemetricname = '120_L_wateredgethresh.func.gii';
% outputdir = './';
% filestem = '120_L_wateredgethresh3_';
% hem = 'L';
% threshperc = .4;

%stepnum = 30000;
fracmaxh = 1;
neighdist = 1;
%minimathresh = .1;
minimathreshperc = .75;
%percentvertices_insideparcel = 1;
edgethresh_tobeinsideparcel = .1;
edgeoutlierthresh = 100;%2.5;
edgeval_hysteresis = 100;%.15;%.165;
minparcelsize = 20;
minisolatedparcelsize = 20;



%edgevalthresh = .13;
edgevalthreshperc = 3/4;
%edgevalmin = .05;
%threshperc = 2/3;
%startingthresh = .09;



bufsize=16384;
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;

%%

medialmaskdata = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii']);
medialmaskdata = medialmaskdata.cdata;
corticalindices = find(medialmaskdata==0);
medialindices = find(medialmaskdata);

erodedmedialmaskdata = gifti(['/data/cn4/laumannt/subcortical_mask/' hem '.atlasroi_erode3.32k_fs_LR.shape.gii']);
erodedmedialmaskdata = ~erodedmedialmaskdata.cdata;

goodverts = gifti(['/data/cn4/evan/ROIs/all_meanimage_' hem '.func.gii']); goodverts = goodverts.cdata > 700;

%%


disp('Finding minima')

metric = gifti(edgemetricname);
metric = metric.cdata;

sortedcorticalmetric = sort(metric(corticalindices),'ascend');
minimathresh = sortedcorticalmetric(round(numel(sortedcorticalmetric)*minimathreshperc));
edgevalthresh = sortedcorticalmetric(round(numel(sortedcorticalmetric)*edgevalthreshperc));
thresh = sortedcorticalmetric(round(numel(sortedcorticalmetric)*threshperc));

metric(1:12) = 100000;

metric(erodedmedialmaskdata) = 100000;

%[temp metric] = upper_completion(metric);
clear temp;

%save(gifti(single(metric)),[outputdir '/avgcorrofcorr_smooth2.55_allgrad_L_smooth2.55_edge_avg_uc.func.gii']);


minimametric = zeros(size(metric));
for i = 1:length(metric)
    
    nodeneigh = i;
    newneigh = nodeneigh;
    curneigh = newneigh;
    
    for n = 1:neighdist
        
        for t = 1:length(curneigh)
            newneigh = [newneigh neighbors(curneigh(t),2:7)];
            newneigh(isnan(newneigh)) = [];
        end
        curneigh = setdiff(newneigh,nodeneigh);
        nodeneigh = union(nodeneigh,newneigh);
        
    end
    
    ind = find(nodeneigh==i); %Which nodeneigh is original point
    [minval mini] = min(metric(nodeneigh));
    minindices = find(metric(nodeneigh)==minval);
    
    
    if minval == metric(i)
        minimametric(i) = 1;
        if length(minindices) > 1
            
            minindices(logical(minindices==ind)) = [];
            metric(nodeneigh(minindices)) = minval+.00001;

        end
    end
    
end

clear metric


%%

%multithresh_label = zeros(length(minimametric),length(threshs));

%for threshnum = 1:length(threshs)
%    thresh = threshs(threshnum);

%disp(num2str(thresh))


edgemetric = gifti(edgemetricname); edgemetric = edgemetric.cdata;
origedgemetric = edgemetric;
%edgemetric(logical((edgemetric > edgeval_hysteresis) .* (edgemetric < .4))) = .4;
edgemetric(logical(edgemetric > edgeval_hysteresis)) = edgemetric(logical(edgemetric > edgeval_hysteresis)) + .1;
%minimametric = gifti(minimametricname); minimametric = minimametric.cdata;
minimametric(1:12) = 0;
minimametric(edgemetric>minimathresh) = 0;
minimametric(medialindices) = 0;

edgemetric(isnan(edgemetric)) = 0;

sortedge = unique(sort(edgemetric,'descend'));
%[sortedge sortedgepos] = sort(edgemetric,'descend');
%Label initial markers with unique value
labelnum = nnz(minimametric);
labelpos = find(minimametric==1);
label = zeros(size(minimametric));
%randomlabels = randperm(labelnum);
[ign sortorder] = sort(edgemetric(labelpos));
for j = 1:labelnum;
    %label(labelpos(j)) = randomlabels(j);
    label(labelpos(j)) = sortorder(j);
end

% minh = sortedge(1);
% maxh = sortedge(end);
% 
% stoph = sortedge(round(length(sortedge)*fracmaxh));
% step = (maxh-minh)/stepnum;
% hiter = minh:step:stoph;

hiter = unique(edgemetric(corticalindices));
%hiter = hiter(hiter<edgevalthresh);

watershedzone = zeros(size(label));

%for i = 1:length(sortedge)
for i = 1:length(hiter);
    
    string{i} = ['Running watershed iteration ' num2str(i) ' out of ' num2str(length(hiter))];
    if i==1; fprintf('%s',string{i}); else fprintf([repmat('\b',1,length(string{i-1})) '%s'],string{i}); end
        
    maskmetrics = edgemetric<hiter(i); % Take values in metric less than current iteration    
    maskmetrics = maskmetrics & ~label>0 & ~watershedzone;
    
    maskpos = find(sum(maskmetrics,2)>0);
    maskpos = maskpos(randperm(length(maskpos)));
    
    for m = 1:length(maskpos) %For all nodes at this threshold
        
        nodeneigh = neighbors(maskpos(m),2:7);
        maskinthismetric = maskmetrics(maskpos(m),:);
        %nodeneigh = neighbors(sortedgepos(i),2:7);
        nodeneigh(isnan(nodeneigh)) = [];
        nodeneighlab = label(nodeneigh);
        
        minfindnodeneighlab = nodeneighlab;
        minfindnodeneighlab(nodeneighlab==0) = 100000; 
        minnodeneighlab = min(minfindnodeneighlab,[],1);
        
        %Find maximum value other than 0 among neighbors
        maxfindnodeneighlab = nodeneighlab;
        maxfindnodeneighlab(nodeneighlab==0) = -100000;
        maxnodeneighlab = max(maxfindnodeneighlab,[],1);
       
        %If min and max differ (i.e. two or more neighbor water), watershed
        %zone
        watershed_nodes = (minnodeneighlab~=maxnodeneighlab) & (minnodeneighlab~=100000) & maskinthismetric;
        label(maskpos(m),watershed_nodes) = 0;
        watershedzone(maskpos(m),watershed_nodes) = 1;
        
        %If min and max the same but different from 0, add to neighbor
        %water
        next_to_water = (minnodeneighlab==maxnodeneighlab) & (minnodeneighlab~=100000) & maskinthismetric;
        label(maskpos(m),next_to_water) = minnodeneighlab(next_to_water);
        

    end
end

disp(' ')

label(erodedmedialmaskdata) = 0;

%label(edgemetric>edgevalthresh) = 0;

origlabel = label;



%%



label = origlabel;


    fprintf('%s',['Merging watersheds: .....'])
%     string{iteration} = ['Merging watersheds: iteration ' num2str(iteration) '.....'];
%     if iteration==1; fprintf('%s',string{iteration}); else fprintf([repmat('\b',1,length(string{iteration-1})) '%s'],string{iteration}); end
    

watersheds = unique(label);
watersheds(watersheds==0) = [];
for waternum = 1:length(watersheds)
    
    adjacentwatersheds{waternum} = [];
%    edgeval_withinwater(waternum) = mean(edgemetric(label==watersheds(waternum)),1);
     valueswithinthiswatershed = sort(origedgemetric(label==watersheds(waternum)));
     valueswithinthiswatershed(valueswithinthiswatershed > edgethresh_tobeinsideparcel) = [];
     %edgeval_withinwater(waternum) = mean(valueswithinthiswatershed(1:(round(length(valueswithinthiswatershed)*percentvertices_insideparcel))),1);
     edgeval_withinwater(waternum) = mean(valueswithinthiswatershed);
    
    for waternum2 = 1:length(watersheds)
        watershedborders{waternum,waternum2} = [];
    end
    
end

borderindices = intersect(find(label==0),corticalindices)';
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
merges = 0;
done = 0;
watersheds_exist = ones(length(watersheds),1);

while done ==0
medians = [];
ijmat = [];
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

[smallest_median smallestind] = min(medians);
% if any(watershedborders{ijmat(smallestind,1),ijmat(smallestind,2)}==28342)
%     1;
% end

if smallest_median < thresh
    merges = merges+1;
    string{merges} = [num2str(merges) ' watersheds merged'];
    if merges==1; fprintf('%s',string{merges}); else fprintf([repmat('\b',1,length(string{merges-1})) '%s'],string{merges}); end
    
    wateri = ijmat(smallestind,1); waterj = ijmat(smallestind,2); 
    label(label==mean(label(origlabel==watersheds(waterj)))) = mean(label(origlabel==watersheds(wateri)),1);
    watersheds_exist(waterj) = 0;
    for k = 1:length(watersheds)
        watershedborders{wateri,k} = unique([watershedborders{wateri,k} watershedborders{waterj,k} watershedborders{k,wateri} watershedborders{k,waterj}]);
        watershedborders{k,wateri} = watershedborders{wateri,k};
    end
else
    done = 1;
end
end
    
    
            
    %fprintf('%s',[num2str(merges) ' watersheds merged'])
    %string{iteration} = [string{iteration} num2str(merges) ' watersheds merged'];
    disp(' ')
    
    for bordervertex = borderindices
        
        borderneighs = neighbors(bordervertex,2:7);
        borderneighs(isnan(borderneighs)) = [];
        borderneighvals = label(borderneighs);
        borderneighvals(borderneighvals==0) = [];
        if length(unique(borderneighvals)) == 1;
            label(bordervertex) = unique(borderneighvals);
        end
    end
    
    clear watershedborders adjacentwatersheds edgeval_withinwater
    
%     if merges > 0
%         save(gifti(label),[outputdir '/' filestem 'watershedmerge_iter' num2str(iteration) '.func.gii']);
%         
%     end
  

%disp(' ')
disp('Merging small parcels')

for repeatthis = 1:2

watersheds = unique(label);
watersheds(watersheds==0) = [];
for waternum = 1:length(watersheds)
    
    adjacentwatersheds{waternum} = [];
    %edgeval_withinwater(waternum) = mean(edgemetric(label==watersheds(waternum)),1);
    valueswithinthiswatershed = sort(origedgemetric(label==watersheds(waternum)));
    valueswithinthiswatershed(valueswithinthiswatershed > edgethresh_tobeinsideparcel) = [];
    %edgeval_withinwater(waternum) = mean(valueswithinthiswatershed(1:(round(length(valueswithinthiswatershed)*percentvertices_insideparcel))),1);
    edgeval_withinwater(waternum) = mean(valueswithinthiswatershed);
    
    for waternum2 = 1:length(watersheds)
        watershedborders{waternum,waternum2} = [];
    end
    
end

borderindices = intersect(find(label==0),corticalindices)';
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
    
    if length(watershedneighbors) == 2
        
        watershedborders{thiswatershedindex(1),thiswatershedindex(2)} = unique([watershedborders{thiswatershedindex(1),thiswatershedindex(2)} bordervertex]);
        watershedborders{thiswatershedindex(2),thiswatershedindex(1)} = unique([watershedborders{thiswatershedindex(2),thiswatershedindex(1)} bordervertex]);
        
    end
    
end


for waternum = 1:length(watersheds)
    
    if (length(find(label==watersheds(waternum))) > 0) && (length(find(label==watersheds(waternum))) < minparcelsize);
        minedge = 1;
        minwaternum = 0;
        for i = 1:length(adjacentwatersheds{waternum})
            thisadjacentwatershedindex = find(watersheds == adjacentwatersheds{waternum}(i));
            if mean(edgemetric(watershedborders{waternum,thisadjacentwatershedindex}),1) < minedge
                minwaternum = thisadjacentwatershedindex;
                minedge = mean(edgemetric(watershedborders{waternum,thisadjacentwatershedindex}),1);
            end
        end
        if minwaternum
            label(label==watersheds(waternum)) = watersheds(minwaternum);
%         else
%             label(label==watersheds(waternum)) = 0;
        end
    end
end

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

clear adjacentwatersheds edgeval_withinwater watershedborders thiswatershedindex adjacentwatersheds watershedborders

disp('Eliminating high edge value vertices and eliminating new small parcels')

label(edgemetric>edgevalthresh) = 0;



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

watersheds = unique(label);
watersheds(watersheds==0) = [];
for waternum = 1:length(watersheds)
    
    adjacentwatersheds{waternum} = [];
    %edgeval_withinwater(waternum) = mean(edgemetric(label==watersheds(waternum)),1);
    valueswithinthiswatershed = sort(origedgemetric(label==watersheds(waternum)));
    valueswithinthiswatershed(valueswithinthiswatershed > edgethresh_tobeinsideparcel) = [];
    %edgeval_withinwater(waternum) = mean(valueswithinthiswatershed(1:(round(length(valueswithinthiswatershed)*percentvertices_insideparcel))),1);
    edgeval_withinwater(waternum) = mean(valueswithinthiswatershed);
    
    for waternum2 = 1:length(watersheds)
        watershedborders{waternum,waternum2} = [];
    end
    
end

borderindices = intersect(find(label==0),corticalindices)';
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
    
    if length(watershedneighbors) == 2
        
        watershedborders{thiswatershedindex(1),thiswatershedindex(2)} = unique([watershedborders{thiswatershedindex(1),thiswatershedindex(2)} bordervertex]);
        watershedborders{thiswatershedindex(2),thiswatershedindex(1)} = unique([watershedborders{thiswatershedindex(2),thiswatershedindex(1)} bordervertex]);
        
    end
    
end

templabel = label;

for waternum = 1:length(watersheds)
    
    if any(find(label==watersheds(waternum))==21752)
        1;
    end
    
    if (length(find(label==watersheds(waternum))) > 0) && (length(find(label==watersheds(waternum))) < minparcelsize);
        minedge = 1;
        minwaternum = 0;
        for i = 1:length(adjacentwatersheds{waternum})
            thisadjacentwatershedindex = find(watersheds == adjacentwatersheds{waternum}(i));
            if ~isempty(watershedborders{waternum,thisadjacentwatershedindex}) && (median(edgemetric(watershedborders{waternum,thisadjacentwatershedindex}),1) < minedge)% && (nnz(edgemetric(watershedborders{waternum,thisadjacentwatershedindex}) < thresh) > 2)
                minwaternum = thisadjacentwatershedindex;
                minedge = mean(edgemetric(watershedborders{waternum,thisadjacentwatershedindex}),1);
            end
        end
        if minwaternum
            label(label==watersheds(waternum)) = mean(label(templabel==watersheds(minwaternum)));
            label(watershedborders{waternum,minwaternum}) = mean(label(templabel==watersheds(minwaternum)));
%         else
%             label(label==watersheds(waternum)) = 0;
        end
    end
end
clear templabel

for bordervertex = borderindices
    
    borderneighs = neighbors(bordervertex,2:7);
    borderneighs(isnan(borderneighs)) = [];
    borderneighvals = label(borderneighs);
    borderneighvals(borderneighvals==0) = [];
    if length(unique(borderneighvals)) == 1;
        label(bordervertex) = unique(borderneighvals);
    end
end




label(edgemetric>edgevalthresh) = 0;

verts_inparcels = find(label);
for vert = verts_inparcels(:)'
    vertneighs = neighbors(vert,2:7); vertneighs(isnan(vertneighs)) = [];
    nonzerovertneighs = vertneighs(label(vertneighs)>0);
    if length(nonzerovertneighs)==2
        neighofneigh = neighbors(nonzerovertneighs(1),2:7); neighofneigh(isnan(neighofneigh))=[];
        if isempty(intersect(neighofneigh,nonzerovertneighs(2)))
            label(vert) = 0;
        end
    end
end

%%
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


watersheds = unique(label); watersheds(watersheds==0) = [];

for watershed = watersheds'
    if nnz(label==watershed) < minisolatedparcelsize
        label(label==watershed) = 0;
    end
end


label(medialindices) = 0;
%%


save(gifti(label),[outputdir '/' filestem 'watershedmerge_' num2str(round(threshperc*100)/100) '.func.gii']);

numparcels = nnz(unique(label));
disp(['Final number of parcels: ' num2str(numparcels)])



