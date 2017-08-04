function watershed_algorithm_submatchtogroup_merge2(edgemetricname,outputdir,hem,subname)
%watershed_algorithm_submatchtogroup_merge2(edgemetricname,outputdir,hem,subname)

% edgemetricname = '120_L_wateredgethresh.func.gii';
% outputdir = './';
% filestem = '120_L_wateredgethresh3_';
% hem = 'L';
% threshperc = .4;

%stepnum = 30000;
fracmaxh = 1;
neighdist = 1;
%minimathresh = .1;
minimathreshperc = .5;
%percentvertices_insideparcel = 1;
edgethresh_tobeinsideparcel = .1;
edgeoutlierthresh = 100;%2.5;
edgeval_hysteresis = 100;%.15;%.165;
minparcelsize = 40;
minisolatedparcelsize = 15;
%edgevalthresh = .13;
edgevalthreshperc = 3/4;
%edgevalmin = .05;
%threshperc = 2/3;
%startingthresh = .09;

max_centroiddistance = 20;

min_correl = .05;

expandifneeded = 1;



bufsize=16384;
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;

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
%hiter = hiter(hiter<=edgevalthresh);

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

subtmask = ['/data/cn4/laumannt/longRestingState/Subjects/BeckySubs/Interleave/FCPROCESS_NEW/gradients_cifti_erode_concat_wateredge/' subname '/allsubs_total_tmask.txt'];

parcelfilename = ['/data/cn4/evan/RestingState/FC_Mapping_120/120_' hem '_wateredgethresh_watershedmerge_0.45.func.gii'];
parcels = gifti(parcelfilename); parcels = parcels.cdata;
origparcels = parcels;

randomizedparcels = gifti(['/data/cn4/evan/RestingState/FC_Mapping_120/120_' hem '_wateredgethresh_watershedmerge_0.45_randomized.func.gii']);
randomizedparcels = randomizedparcels.cdata;

parcelcorrelpatternfile = ['/data/cn4/evan/RestingState/FC_Mapping_120/Group_parcelcorrelpatterns_' hem '.mat'];

parcelIDs = unique(parcels);
parcelIDs(parcelIDs==0)=[];

load(parcelcorrelpatternfile)
clear neighbors
bufsize=16384;
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;

maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii']; mask = gifti(maskname); mask = ~mask.cdata;
    parcels = parcels(logical(mask));
    maskL = gifti('/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii');ncortexLverts = nnz(maskL.cdata==0);
    
    load(['/data/cn4/evan/fsaverage_LR32k/Surface_distances_' hem '.mat']);
surfdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
midthick = gifti([surfdir '/Conte69.' hem '.midthickness.32k_fs_LR.surf.gii']);
sphere = gifti([surfdir '/Conte69.' hem '.sphere.32k_fs_LR.surf.gii']);
group_parcel_centroid_vert = zeros(length(parcelIDs),1);

for parcelnum = 1:length(parcelIDs)
    ind = find(origparcels==parcelIDs(parcelnum));
    randomizedIDs(parcelnum) = mean(randomizedparcels(ind));
    parcelindices{parcelnum} = find(parcels==parcelIDs(parcelnum)) + (strcmp(hem,'R') * ncortexLverts);
    
    
    meanX = mean(sphere.vertices(ind,1));
    meanY = mean(sphere.vertices(ind,2));
    meanZ = mean(sphere.vertices(ind,3));
    
    coord = [meanX meanY meanZ];
    sphere_coords = [sphere.vertices(ind,1) sphere.vertices(ind,2) sphere.vertices(ind,3)];
    
    rep_coord = repmat(coord, [size(sphere_coords,1) 1]);
    
    dist_coord = sum((sphere_coords-rep_coord).^2,2).^(1/2);
    [y indval] = min(dist_coord);
    group_parcel_centroid_vert(parcelnum) = ind(indval);
end


subassignmap = zeros(nnz(mask),1);


subparcels = label;

origsubparcels = subparcels;
    subparcels = subparcels(logical(mask));
    subparcelIDs = unique(subparcels); subparcelIDs(subparcelIDs==0)=[];
    
    groupassignmentIDs = zeros(length(parcelIDs),1);
    subassignmentIDs = zeros(length(subparcelIDs),1);
    
    evalc(['!wb_command -cifti-convert -to-gifti-ext /data/hcp-zfs/home/laumannt/LFRS_parcellation/' subname '/' subname '_timeseries_normalwall.dtseries.nii Temp.func.gii']);
    subtimecourse = gifti('Temp.func.gii'); subtimecourse = subtimecourse.cdata;
    delete('Temp.func.gii*');
        tmask = load(subtmask);
        subtimecourse = subtimecourse(:,logical(tmask));
    subtimecourse(isnan(subtimecourse)) = 0;
    
     groupvsubcorrel = zeros(length(parcelIDs),length(subparcelIDs));
    
    sub_parcel_centroid_vert = zeros(length(subparcelIDs),1);
    
    for subparcelnum = 1:length(subparcelIDs)
        string{subparcelnum} = ['  Parcel number ' num2str(subparcelnum) ' out of ' num2str(length(subparcelIDs))];
        if subparcelnum==1; fprintf('%s',string{subparcelnum}); else fprintf([repmat('\b',1,length(string{subparcelnum-1})) '%s'],string{subparcelnum}); end
        
        subparcelindices{subparcelnum} = find(subparcels==subparcelIDs(subparcelnum)) + (strcmp(hem,'R') * ncortexLverts);
        
        parceltimecourse = mean(subtimecourse(subparcelindices{subparcelnum},:),1);
        subparcelcorrelpattern = paircorr_mod(parceltimecourse',subtimecourse');
        subparcelcorrelpattern(isnan(subparcelcorrelpattern)) = 0;
        
        ind = find(origsubparcels==subparcelIDs(subparcelnum));
        meanX = mean(sphere.vertices(ind,1));
        meanY = mean(sphere.vertices(ind,2));
        meanZ = mean(sphere.vertices(ind,3));
        
        coord = [meanX meanY meanZ];
        sphere_coords = [sphere.vertices(ind,1) sphere.vertices(ind,2) sphere.vertices(ind,3)];
        
        rep_coord = repmat(coord, [size(sphere_coords,1) 1]);
        
        dist_coord = sum((sphere_coords-rep_coord).^2,2).^(1/2);
        [y indval] = min(dist_coord);
        sub_parcel_centroid_vert(subparcelnum) = ind(indval);
        
        thisparcel_distances = geo_distances(group_parcel_centroid_vert,sub_parcel_centroid_vert(subparcelnum));
        
        
        for parcelnum = 1:length(parcelIDs)
            indicestoavoid = union(find(geo_distances(:,group_parcel_centroid_vert(parcelnum))<20), find(geo_distances(:,sub_parcel_centroid_vert(subparcelnum))<20));
            indicestocorrelate = logical(ones(length(subparcelcorrelpattern),1)); indicestocorrelate(indicestoavoid) = 0;
            groupvsubcorrel(parcelnum,subparcelnum) = paircorr_mod(parcelcorrelpatterns(parcelnum,indicestocorrelate)',subparcelcorrelpattern(indicestocorrelate)');
        end
        
        
        thisparcel_groupvsubcorrel = groupvsubcorrel(:,subparcelnum);
        thisparcel_groupvsubcorrel(thisparcel_distances>max_centroiddistance) = -Inf;
        thisparcel_groupvsubcorrel(thisparcel_groupvsubcorrel<min_correl) = -Inf;
        
        [correl groupparcel] = max(thisparcel_groupvsubcorrel);
        if expandifneeded
            expanded_centroiddistance = max_centroiddistance;
            while correl==-Inf
                expanded_centroiddistance = expanded_centroiddistance + 5;
                thisparcel_groupvsubcorrel = groupvsubcorrel(:,subparcelnum);
                thisparcel_groupvsubcorrel(thisparcel_distances>expanded_centroiddistance) = -Inf;
                %if expanded_centroiddistance < 50
                    thisparcel_groupvsubcorrel(thisparcel_groupvsubcorrel<min_correl) = -Inf;
                %end
                [correl groupparcel] = max(thisparcel_groupvsubcorrel);
            end
        elseif correl==-Inf
            error(['No parcels within ' num2str(max_centroiddistance) ' vertices!'])
        end
        subassignment(subparcelnum) = groupparcel;

        if groupassignmentIDs(groupparcel)
            subassignmentIDs(subparcelnum) = randomizedIDs(groupparcel);
            groupassignment{groupparcel} = [groupassignment{groupparcel} subparcelnum];
        else
            subassignmentIDs(subparcelnum) = randomizedIDs(groupparcel);
            groupassignmentIDs(groupparcel) = randomizedIDs(groupparcel);
            groupassignment{groupparcel} = subparcelnum;
        end
        
    end
    disp(' ')
    
    for parcelnum = 1:length(subparcelIDs)
        subassignmap(subparcels==subparcelIDs(parcelnum)) = subassignmentIDs(parcelnum);
    end
    temp = zeros(32492,1); temp(logical(mask)) = subassignmap; subassignmap = temp;
    %%
    
    edgethresh_subassignmap = subassignmap;
    edgethresh_subassignmap(edgemetric>edgevalthresh) = 0;
    
for waternum = 1:length(parcelIDs)
    
    for waternum2 = 1:length(parcelIDs)
        watershedborders{waternum,waternum2} = [];
    end
    
    assignval(waternum) = mean(subassignmap(origlabel==parcelIDs(waternum)));
    
end

borderindices = intersect(find(subassignmap==0),corticalindices)';
for bordervertex = borderindices
        
        borderneighs = neighbors(bordervertex,2:7);
        borderneighs(isnan(borderneighs)) = [];
        borderneighvals = origparcels(borderneighs);
        borderneighvals(borderneighvals==0) = [];
        
        watershedneighbors = unique(borderneighvals);
        
        for waterneighbornum = 1:length(watershedneighbors)
            
            thiswatershedindex(waterneighbornum) = find(parcelIDs==watershedneighbors(waterneighbornum));
                        
        end
        
        if length(watershedneighbors) == 2 && length(borderneighvals)>2
        
            watershedborders{thiswatershedindex(1),thiswatershedindex(2)} = unique([watershedborders{thiswatershedindex(1),thiswatershedindex(2)} bordervertex]);
            watershedborders{thiswatershedindex(2),thiswatershedindex(1)} = unique([watershedborders{thiswatershedindex(2),thiswatershedindex(1)} bordervertex]);
        
        end        
end

for wateri = 1:length(parcelIDs)
    for waterj = 1:length(parcelIDs)
        if wateri < waterj
        
        if (length(watershedborders{wateri,waterj}) > 0) && (assignval(wateri)==assignval(waterj))
            
            edgethresh_subassignmap(watershedborders{wateri,waterj}) = assignval(wateri);
            
        end
        end
    end
end 
    
    
    zeroinds = find((edgethresh_subassignmap==0) .* (mask==1));
    
    for ind = zeroinds(:)'
        indneighs = neighbors(ind,2:7); indneighs(isnan(indneighs)) = [];
        indneighs_2away = neighbors(indneighs,2:7); indneighs_2away(isnan(indneighs_2away)) = [];
        
        nonzeroneighvals = unique(edgethresh_subassignmap(indneighs_2away)); nonzeroneighvals(nonzeroneighvals==0) = [];
        if length(nonzeroneighvals) == 1
            edgethresh_subassignmap(ind) = nonzeroneighvals;
        end
    end


    label = edgethresh_subassignmap;
    



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
    
    sizes(waternum) = nnz(label==watersheds(waternum));
    IDs(waternum) = mean(edgethresh_subassignmap(label==watersheds(waternum)));
    
end

[ign sorti] = sort(sizes,'descend');
sortedIDs = IDs(sorti);
sortedclusters = watersheds(sorti);

for waternum = 1:length(sortedIDs)
    if any(sortedIDs(1:waternum-1) == sortedIDs(waternum));
        edgethresh_subassignmap(label==sortedclusters(waternum)) = 0;
    end
end
        
   


label(medialindices) = 0;

filestem = subname;

save(gifti(edgethresh_subassignmap),[outputdir '/' filestem '_' hem '_watershedmatchtogroup2.func.gii']);



