max_centroiddistance = 50;

%min_correlval = 0;%.25;

hem = 'R';

iscifti = 1;

outputfolder = '/data/cn4/evan/RestingState/FC_Mapping_120/Cohorts/'; %location watersheds and match data will be written to


parcelfilename1 = '/data/cn4/evan/RestingState/FC_Mapping_120/Cohorts/C1_Rwatershedmerge_0.45.func.gii';

cormapfilename1 = '/data/cn4/evan/RestingState/FC_Mapping_120/Cohorts/C1_avg_corr_R_normalwall.mat';



parcelfilename2 = '/data/cn4/evan/RestingState/FC_Mapping_120/Cohorts/C2_Rwatershedmerge_0.45.func.gii';

cormapfilename2 = '/data/cn4/evan/RestingState/FC_Mapping_120/Cohorts/C2_avg_corr_R_normalwall.mat';





parcels1 = gifti(parcelfilename1); parcels1 = parcels1.cdata;
origparcels1 = parcels1;
parcelIDs1 = unique(parcels1);
parcelIDs1(parcelIDs1==0)=[];

parcels2 = gifti(parcelfilename2); parcels2 = parcels2.cdata;
origparcels2 = parcels2;
parcelIDs2 = unique(parcels2); parcelIDs2(parcelIDs2==0)=[];

bufsize=16384;
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;

if iscifti==1
    maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii']; mask = gifti(maskname); mask = ~mask.cdata;
    parcels1 = parcels1(logical(mask));
    parcels2 = parcels2(logical(mask));
    maskL = gifti('/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii');ncortexLverts = nnz(maskL.cdata==0);
elseif iscifti==2
    maskname = ['/data/cn4/laumannt/subcortical_mask/' hem '.atlasroi_erode3.32k_fs_LR.shape.gii']; mask = gifti(maskname); mask = mask.cdata;
    parcels1 = parcels1(logical(mask));
    parcels2 = parcels2(logical(mask));
    maskL = gifti('/data/cn4/laumannt/subcortical_mask/L.atlasroi_erode3.32k_fs_LR.shape.gii');ncortexLverts = nnz(maskL.cdata);
    maskR = gifti('/data/cn4/laumannt/subcortical_mask/R.atlasroi_erode3.32k_fs_LR.shape.gii');ncortexRverts = nnz(maskR.cdata);
else
    parcels1 = parcels_upsampled;
    ncortexLverts = length(parcels1);
end

load(['/data/cn4/evan/fsaverage_LR32k/Surface_distances_' hem '.mat']);

surfdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
sphere = gifti([surfdir '/Conte69.' hem '.sphere.32k_fs_LR.surf.gii']);

disp('Calculating parcel connectivities for map 1')

%cormap1 = gifti(cormapfilename1); cormap1 = cormap1.cdata;
load(cormapfilename1);
cormap1 = avg_corr_R; clear avg_corr_R

parcel1_centroid_vert = zeros(length(parcelIDs1),1);
parcelcorrelpattern1 = zeros(length(parcelIDs1),size(cormap1,2));
for parcelnum = 1:length(parcelIDs1)
    %disp(['  Parcel number ' num2str(parcelnum)])
    
    parcelindices1{parcelnum} = find(parcels1==parcelIDs1(parcelnum)); %+ (strcmp(hem,'R') * ncortexLverts);
     
    parcelcorrelpattern1(parcelnum,:) = mean(cormap1(parcelindices1{parcelnum},:),1);
    
    ind = find(origparcels1==parcelIDs1(parcelnum));
    meanX = mean(sphere.vertices(ind,1));
    meanY = mean(sphere.vertices(ind,2));
    meanZ = mean(sphere.vertices(ind,3));
    
    coord = [meanX meanY meanZ];
    sphere_coords = [sphere.vertices(ind,1) sphere.vertices(ind,2) sphere.vertices(ind,3)];
    
    rep_coord = repmat(coord, [size(sphere_coords,1) 1]);
    
    dist_coord = sum((sphere_coords-rep_coord).^2,2).^(1/2);
    [y indval] = min(dist_coord);
    parcel1_centroid_vert(parcelnum) = ind(indval);
end

parcelcorrelpattern1(isnan(parcelcorrelpattern1)) = 0;

clear cormap1

disp('Calculating parcel connectivities for map 2')

%cormap2 = gifti(cormapfilename2); cormap2 = cormap2.cdata;
load(cormapfilename2);
cormap2 = avg_corr_R; clear avg_corr_R

parcel2_centroid_vert = zeros(length(parcelIDs2),1);
parcelcorrelpattern2 = zeros(length(parcelIDs2),size(cormap2,2));
for parcelnum = 1:length(parcelIDs2)
    %disp(['  Parcel number ' num2str(parcelnum)])
    
    parcelindices2{parcelnum} = find(parcels2==parcelIDs2(parcelnum)); %+ (strcmp(hem,'R') * ncortexLverts);
     
    parcelcorrelpattern2(parcelnum,:) = mean(cormap2(parcelindices2{parcelnum},:),1);
    
    ind = find(origparcels2==parcelIDs2(parcelnum));
    meanX = mean(sphere.vertices(ind,1));
    meanY = mean(sphere.vertices(ind,2));
    meanZ = mean(sphere.vertices(ind,3));
    
    coord = [meanX meanY meanZ];
    sphere_coords = [sphere.vertices(ind,1) sphere.vertices(ind,2) sphere.vertices(ind,3)];
    
    rep_coord = repmat(coord, [size(sphere_coords,1) 1]);
    
    dist_coord = sum((sphere_coords-rep_coord).^2,2).^(1/2);
    [y indval] = min(dist_coord);
    parcel2_centroid_vert(parcelnum) = ind(indval);
end

clear cormap2

onevtwo_pctoverlap = zeros(length(parcelIDs1),length(parcelIDs2));
for i = 1:length(parcelIDs1)
    for j = 1:length(parcelIDs2)
        numoverlapverts = length(intersect(parcelindices1{i},parcelindices2{j}));
        %onevtwo_pctoverlap(i,j) = mean([numoverlapverts/length(parcelindices1{i}),numoverlapverts/length(parcelindices2{j})]);
        onevtwo_pctoverlap(i,j) = numoverlapverts;
    end
end

 
 disp('Conducting assignments and saving outputs')
 
 onevtwocorrel = paircorr_mod(parcelcorrelpattern1',parcelcorrelpattern2');

onevtwo_distances = geo_distances(parcel1_centroid_vert,parcel2_centroid_vert);

%onevtwocost = ((max(max(onevtwo_pctoverlap))-onevtwo_pctoverlap) / (max(max(onevtwo_pctoverlap)))) + (1-(onevtwocorrel+1)/2);
onevtwocost = (onevtwo_distances ./ (max(max(onevtwo_distances))))*.5 + ((max(max(onevtwo_pctoverlap))-onevtwo_pctoverlap) / (max(max(onevtwo_pctoverlap))))*.5 + (1-(onevtwocorrel+1)/2);
%onevtwocost(onevtwo_pctoverlap==0) = Inf;
%onevtwocost = (1-(onevtwocorrel+1)/2);
onevtwocost(onevtwo_distances > max_centroiddistance) = Inf;
%onevtwocost(onevtwocorrel < min_correlval) = Inf;

%calculate the assignments
[oneassignment, twoassignment, cost1, cost2] = munkres_mult(onevtwocost,4);

randomizedIDs = randperm(length(parcelIDs1)+length(parcelIDs2));

assignment1IDs = zeros(length(parcelIDs1),1);
assignment2IDs = zeros(length(parcelIDs2),1);

% for i=1:length(assignment1IDs)
%     if oneassignment(i) > 0
%     if assignment2IDs(oneassignment(i))==0
%        assignment1IDs(i) = randomizedIDs(i);
%        assignment2IDs(oneassignment(i)) = randomizedIDs(i);
%     else
%         assignment1IDs(i) = assignment2IDs(oneassignment(i));
%     end
%     end
% end
% 
% for i=1:length(assignment2IDs)
%     if twoassignment(i) > 0
%     if assignment1IDs(twoassignment(i))==0
%        assignment2IDs(i) = randomizedIDs(i+length(parcelIDs1));
%        assignment1IDs(twoassignment(i)) = randomizedIDs(i+length(parcelIDs1));
%     else
%        assignment2IDs(i) = assignment1IDs(twoassignment(i));
%     end
%     end
% end
%        
% for i = 1:length(randomizedIDs)
%      if nnz(assignment1IDs==randomizedIDs(i)) > 0 && nnz(assignment2IDs==randomizedIDs(i)) ==0
%          assignment2IDs(oneassignment(assignment1IDs==randomizedIDs(i))) = randomizedIDs(i);
%      end
%      if nnz(assignment2IDs==randomizedIDs(i)) > 0 && nnz(assignment1IDs==randomizedIDs(i)) ==0
%          assignment1IDs(twoassignment(assignment2IDs==randomizedIDs(i))) = randomizedIDs(i);
%      end
% end   


for i = 1:length(parcelIDs1)
    if assignment1IDs(i)==0
        done = 0;
        oneinds = [i];
        twoinds = [];
        while done == 0;
            thesetwoinds = unique([twoinds oneassignment(oneinds) find(ismember(twoassignment,oneinds))]);
            theseoneinds = unique([oneinds twoassignment(thesetwoinds) find(ismember(oneassignment,thesetwoinds))]);
            if isequal(oneinds,theseoneinds) && isequal(twoinds,thesetwoinds)
                done = 1;
                assignment1IDs(theseoneinds) = randomizedIDs(i);
                assignment2IDs(thesetwoinds) = randomizedIDs(i);
            else
                oneinds = theseoneinds;
                twoinds = thesetwoinds;
            end
        end
    end
end

orig_assignment1IDs = assignment1IDs;
orig_assignment2IDs = assignment2IDs;

for i = 1:length(randomizedIDs)
    if nnz(assignment1IDs==randomizedIDs(i)) > 1 && nnz(assignment2IDs==randomizedIDs(i)) > 1
        oneinds = find(assignment1IDs==randomizedIDs(i)); twoinds = find(assignment2IDs==randomizedIDs(i));
        [tempassignment ign] = munkres(onevtwocost(oneinds,twoinds));
        [temp1assignment temp2assignment ign ign] = munkres_mult(onevtwocost(oneinds,twoinds),4);
        
        for j = 2:length(tempassignment)
            if tempassignment(j)~=0
                unusedIDs = setdiff(randomizedIDs,assignment1IDs);
                random_unusedID = unusedIDs(randi(length(unusedIDs),1));
                assignment1IDs(oneinds(j)) = random_unusedID;
                assignment2IDs(twoinds(tempassignment(j))) = random_unusedID;
                oneassignment(oneinds(j)) = twoinds(tempassignment(j));
                twoassignment(twoinds(tempassignment(j))) = oneinds(j);
            end
        end
        
        for j = 2:length(tempassignment)
            if tempassignment(j)==0
                oneassignment(oneinds(j)) = twoinds(temp1assignment(j));
                assignment1IDs(oneinds(j)) = assignment2IDs(twoinds(temp1assignment(j)));
                
            end
        end
        
        for j = 1:length(twoinds)
            if ~any(tempassignment==j)
                twoassignment(twoinds(j)) = oneinds(temp2assignment(j));
                assignment2IDs(twoinds(j)) = assignment1IDs(oneinds(temp2assignment(j)));
            end
        end
    end
    
end

for i = 1:length(assignment1IDs);
    if ~any(assignment2IDs==assignment1IDs(i))
        possibleassigns = find(orig_assignment2IDs==orig_assignment1IDs(i));
        assignment1IDs(i) = assignment2IDs(possibleassigns(1));
    end
end

for i = 1:length(assignment2IDs);
    if ~any(assignment1IDs==assignment2IDs(i))
        possibleassigns = find(orig_assignment1IDs==orig_assignment2IDs(i));
        assignment2IDs(i) = assignment1IDs(possibleassigns(1));
        
    end
end





oneassignmap = zeros(nnz(mask),1);
onesimilaritymap = zeros(nnz(mask),1);
%onedistancemap = zeros(nnz(mask),1);
twoassignmap = zeros(nnz(mask),1);
twosimilaritymap = zeros(nnz(mask),1);
%twodistancemap = zeros(nnz(mask),1);

for parcelnum = 1:length(parcelIDs1)
    
    if (oneassignment(parcelnum) > 0) %&& (onevtwo_pctoverlap(parcelnum,oneassignment(parcelnum)) > 0)
        
        oneassignmap(parcels1==parcelIDs1(parcelnum)) = assignment1IDs(parcelnum);
        if isempty(find(twoassignment==parcelnum))
            onesimilaritymap(parcels1==parcelIDs1(parcelnum)) = onevtwocorrel(parcelnum,oneassignment(parcelnum));
        else
            onesimilaritymap(parcels1==parcelIDs1(parcelnum)) = mean(onevtwocorrel(parcelnum,(twoassignment==parcelnum)));
        end
        
    else
        
        oneassignmap(parcels1==parcelIDs1(parcelnum)) = -1;
        onesimilaritymap(parcels1==parcelIDs1(parcelnum)) = -1;
    end
end

for parcelnum = 1:length(parcelIDs2)
    
    if twoassignment(parcelnum) > 0  %&& (onevtwo_pctoverlap(twoassignment(parcelnum),parcelnum) > 0)
        
%        onedistancemap(parcels1==parcelIDs1(parcelnum)) = onevtwo_distances(parcelnum,assignment(parcelnum));
        
        twoassignmap(parcels2==parcelIDs2(parcelnum)) = assignment2IDs(parcelnum);
        
        if isempty(find(oneassignment==parcelnum))
            twosimilaritymap(parcels2==parcelIDs2(parcelnum)) = onevtwocorrel(twoassignment(parcelnum),parcelnum);
        else
            twosimilaritymap(parcels2==parcelIDs2(parcelnum)) = mean(onevtwocorrel((oneassignment==parcelnum),parcelnum));
        end
%        twodistancemap(parcels2==parcelIDs2(assignment(parcelnum))) = onevtwo_distances(parcelnum,assignment(parcelnum));
        
    else
        
        twoassignmap(parcels2==parcelIDs2(parcelnum)) = -1;
        twosimilaritymap(parcels2==parcelIDs2(parcelnum)) = -1;
        
        
%        onedistancemap(parcels1==parcelIDs1(parcelnum)) = -1;
        
    end
end





assgnoutput1 = zeros(length(mask),1); assgnoutput1(logical(mask),:) = oneassignmap; origassgnoutput1 = assgnoutput1;
assgnoutput2 = zeros(length(mask),1); assgnoutput2(logical(mask),:) = twoassignmap; origassgnoutput2 = assgnoutput2;

simoutput1 = zeros(length(mask),1); simoutput1(logical(mask),:) = onesimilaritymap;
simoutput2 = zeros(length(mask),1); simoutput2(logical(mask),:) = twosimilaritymap;


for vertnum = 1:size(assgnoutput1,1)
    vertneighs = neighbors(vertnum,2:7);
    vertneighs(isnan(vertneighs)) = [];
    if origassgnoutput1(vertnum)==-1 && ~any(origassgnoutput1(vertneighs)==0)
        assgnoutput1(vertnum)=0;
        simoutput1(vertnum) = 0;
    end
    if origassgnoutput2(vertnum)==-1 && ~any(origassgnoutput2(vertneighs)==0)
        assgnoutput2(vertnum)=0;
        simoutput2(vertnum) = 0;
    end
end


save(gifti(single(assgnoutput1)),[outputfolder 'C1_parcels_assignedtoC2_' hem  '.func.gii'])
save(gifti(single(assgnoutput2)),[outputfolder 'C2_parcels_assignedtoC1_' hem '.func.gii'])
save(gifti(single(simoutput1)),[outputfolder 'C1_parcels_assignedtoC2_similarity_' hem '.func.gii'])
save(gifti(single(simoutput2)),[outputfolder 'C2_parcels_assignedtoC1_similarity_' hem '.func.gii'])

%output = zeros(length(mask),1); output(logical(mask),:) = onesimilaritymap;


%output = zeros(length(mask),1); output(logical(mask),:) = onedistancemap;
%save(gifti(single(output)),[outputfolder 'Map1_parcels_assignedtosubs_distance_' hem 'distthresh_' num2str(max_centroiddistance) '.func.gii'])


%output = zeros(length(mask),1); output(logical(mask),:) = twosimilaritymap;


%output = zeros(length(mask),1); output(logical(mask),:) = twodistancemap;
%save(gifti(single(output)),[outputfolder 'Map2_parcels_assignedtogroup_distance_' hem 'distthresh_' num2str(max_centroiddistance) '.func.gii'])


