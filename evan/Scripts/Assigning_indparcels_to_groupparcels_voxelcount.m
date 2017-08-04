
hems = {'L','R'};
for hemnum = 1:2
hem = hems{hemnum};

subjectnums = [];

max_centroiddistance = 20;

datalist = '/data/cn4/evan/RestingState/Ind_variability/Subjects/LFRS_DATALIST.txt';
tmasklist = '/data/cn4/evan/RestingState/Ind_variability/Subjects/LFRS_TMASKLIST.txt';
iscifti = 1;

outputfolder = '/data/cn4/evan/RestingState/Ind_variability/Subjects/Assignments/'; %location match data will be written to

parcelfilename = ['/data/cn4/evan/RestingState/FC_Mapping_120/120_' hem '_wateredgethresh_watershedmerge_0.45.func.gii'];
parcels = gifti(parcelfilename); parcels = parcels.cdata;
origparcels = parcels;

parcelcorrelpatternfile = ['/data/cn4/evan/RestingState/FC_Mapping_120/Group_parcelcorrelpatterns_' hem '.mat'];

% baddata = gifti(['/data/cn4/evan/ROIs/all_meanimage_' hem '.func.gii']); baddata = baddata.cdata;
% baddata = baddata<750;

parcelIDs = unique(parcels);
parcelIDs(parcelIDs==0)=[];

% parcelIDs_incrap = unique(parcels .* baddata);
% for badparcelID = parcelIDs_incrap'
%     parcelIDs(parcelIDs == badparcelID) = [];
% end

load(parcelcorrelpatternfile)
clear neighbors
bufsize=16384;
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;


if iscifti==1
    maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii']; mask = gifti(maskname); mask = ~mask.cdata;
    parcels = parcels(logical(mask));
    maskL = gifti('/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii');ncortexLverts = nnz(maskL.cdata==0);
elseif iscifti==2
    maskname = ['/data/cn4/laumannt/subcortical_mask/' hem '.atlasroi_erode3.32k_fs_LR.shape.gii']; mask = gifti(maskname); mask = mask.cdata;
    parcels = parcels(logical(mask));
    maskL = gifti('/data/cn4/laumannt/subcortical_mask/L.atlasroi_erode3.32k_fs_LR.shape.gii');ncortexLverts = nnz(maskL.cdata);
    maskR = gifti('/data/cn4/laumannt/subcortical_mask/R.atlasroi_erode3.32k_fs_LR.shape.gii');ncortexRverts = nnz(maskR.cdata);
else
    parcels = parcels_upsampled;
    ncortexLverts = length(parcels);
end

load(['/data/cn4/evan/fsaverage_LR32k/Surface_distances_' hem '.mat']);
surfdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
midthick = gifti([surfdir '/Conte69.' hem '.midthickness.32k_fs_LR.surf.gii']);
sphere = gifti([surfdir '/Conte69.' hem '.sphere.32k_fs_LR.surf.gii']);
group_parcel_centroid_vert = zeros(length(parcelIDs),1);

for parcelnum = 1:length(parcelIDs)
    ind = find(origparcels==parcelIDs(parcelnum));
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
%group_center_coords = [midthick.vertices(indpos,1) midthick.vertices(indpos,2) midthick.vertices(indpos,3)];
%clear indpos

[subjects subdata] = textread(datalist,'%s%s');
if ~isempty(tmasklist)
    [subjects tmasks] = textread(tmasklist,'%s%s');
end

if ~isempty(subjectnums)
    subjects = subjects(subjectnums); subdata = subdata(subjectnums); tmasks = tmasks(subjectnums);
end


randomizedIDs = randperm(length(parcelIDs));

groupassignoutput = zeros(length(mask),length(subjects));
groupsimilarityoutput = zeros(length(mask),length(subjects));
groupdistanceoutput = zeros(length(mask),length(subjects));
subassignoutput = zeros(length(mask),length(subjects));
subsimilarityoutput = zeros(length(mask),length(subjects));
subdistanceoutput = zeros(length(mask),length(subjects));

%subject loop
for s = 1:length(subjects)
    
    groupassignmap = zeros(nnz(mask),1);
groupsimilaritymap = zeros(nnz(mask),1);
groupdistancemap = zeros(nnz(mask),1);
subassignmap = zeros(nnz(mask),1);
subsimilaritymap = zeros(nnz(mask),1);
subdistancemap = zeros(nnz(mask),1);
    
    disp(['Subject ' num2str(s)])
    
    subparcelfilename = ['/data/hcp-bluearc/home/laumannt/LFRS_parcellation/' subjects{s} '/' subjects{s} '_' hem '_watershedmerge_0.45.func.gii'];
    subparcels = gifti(subparcelfilename); subparcels = subparcels.cdata;
    origsubparcels = subparcels;
    subparcels = subparcels(logical(mask));
    subparcelIDs = unique(subparcels); subparcelIDs(subparcelIDs==0)=[];
    
    
    
    evalc(['!wb_command -cifti-convert -to-gifti-ext /data/hcp-bluearc/home/laumannt/LFRS_parcellation/' subjects{s} '/' subjects{s} '_timeseries_normalwall.dtseries.nii Temp.func.gii']);
    subtimecourse = gifti('Temp.func.gii'); subtimecourse = subtimecourse.cdata;
    delete('Temp.func.gii*');
    if ~isempty(tmasklist)
        tmask = load(tmasks{s});
        subtimecourse = subtimecourse(:,logical(tmask));
    end
    subtimecourse(isnan(subtimecourse)) = 0;
    
    groupvsub_voxelcountperc = zeros(length(parcelIDs),length(subparcelIDs));
    
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
        
        for parcelnum = 1:length(parcelIDs)
            indicestoavoid = union(find(geo_distances(:,group_parcel_centroid_vert(parcelnum))<20), find(geo_distances(:,sub_parcel_centroid_vert(subparcelnum))<20));
            indicestocorrelate = logical(ones(length(subparcelcorrelpattern),1)); indicestocorrelate(indicestoavoid) = 0;
            
            sortvals = sort(subparcelcorrelpattern(indicestocorrelate),'descend'); subthresh = sortvals(round(kden*length(sortvals)));
            sortvals = sort(parcelcorrelpatterns(parcelnum,indicestocorrelate),'descend'); groupthresh = sortvals(round(kden*length(sortvals)));
            
            intersect = (subparcelcorrelpattern(indicestocorrelate) > subthresh) .* (parcelcorrelpatterns(parcelnum,indicestocorrelate) > groupthresh);
            
            groupvsub_voxelcountperc(parcelnum,subparcelnum) = nnz(intersect) / nnz(indicestocorrelate);
            
        end
        
    end
    
    disp(' ')
    
    groupvsub_overlap = zeros(length(parcelIDs),length(subparcelIDs));
    for i = 1:length(parcelIDs)
        for j = 1:length(subparcelIDs)
            numoverlapverts = length(intersect(parcelindices{i},subparcelindices{j}));
            %groupvsub_overlap(i,j) = numoverlapverts;
            groupvsub_overlap(i,j) = mean([numoverlapverts / length(parcelindices{i}) , numoverlapverts / length(subparcelindices{j})]);
        end
    end
    
     %sub_center_coords = [midthick.vertices(indpos,1) midthick.vertices(indpos,2) midthick.vertices(indpos,3)];
     %clear indpos
    
    group_v_sub_distance = geo_distances(group_parcel_centroid_vert,sub_parcel_centroid_vert);
     
    groupvsubcost = (1-(groupvsub_voxelcountperc+1)/2) ;%+ .5*(1-groupvsub_overlap) + .5*(group_v_sub_distance ./ max_centroiddistance);
    %groupvsubcost = (1-(groupvsubcorrel+1)/2) + .5*((max(max(groupvsub_overlap))-groupvsub_overlap) / (max(max(groupvsub_overlap)))) + .5*(group_v_sub_distance ./ (max(max(group_v_sub_distance))));
    %groupvsubcost = (1-(groupvsubcorrel+1)/2);
    groupvsubcost(group_v_sub_distance > max_centroiddistance) = Inf;
    
    groupvsubcost(groupvsub_voxelcountperc < .1) = Inf;
    
    %groupvsubcost(groupvsubcorrel < min_correlval) = Inf;
    
    %calculate the assignments
    [groupassignment, subassignment, groupcost, subcost] = munkres_mult(groupvsubcost,6);
    
    
    randomizedIDs = randperm(length(parcelIDs)+length(subparcelIDs));

groupassignmentIDs = zeros(length(parcelIDs),1);
subassignmentIDs = zeros(length(subparcelIDs),1);

% for i=1:length(groupassignmentIDs)
%     if groupassignment(i) > 0
%     if subassignmentIDs(groupassignment(i))==0
%        groupassignmentIDs(i) = randomizedIDs(i);
%        subassignmentIDs(groupassignment(i)) = randomizedIDs(i);
%     else
%         groupassignmentIDs(i) = subassignmentIDs(groupassignment(i));
%     end
%     end
% end
% 
% for i=1:length(subassignmentIDs)
%     if subassignment(i) > 0
%     if groupassignmentIDs(subassignment(i))==0
%        subassignmentIDs(i) = randomizedIDs(i+length(parcelIDs));
%        groupassignmentIDs(subassignment(i)) = randomizedIDs(i+length(parcelIDs));
%     else
%        subassignmentIDs(i) = groupassignmentIDs(subassignment(i));
%     end
%     end
% end


for i = 1:length(parcelIDs)
    if groupassignment(i) > 0
    if groupassignmentIDs(i)==0
        done = 0;
        groupinds = [i];
        subinds = [];
        while done == 0;
            thesesubinds = unique([subinds groupassignment(groupinds) find(ismember(subassignment,groupinds))]);
            thesegroupinds = unique([groupinds subassignment(thesesubinds) find(ismember(groupassignment,thesesubinds))]);
            if isequal(groupinds,thesegroupinds) && isequal(subinds,thesesubinds)
                done = 1;
                groupassignmentIDs(thesegroupinds) = randomizedIDs(i);
                subassignmentIDs(thesesubinds) = randomizedIDs(i);
            else
                groupinds = thesegroupinds;
                subinds = thesesubinds;
            end
        end
    end
    else
        groupassignmentIDs(i) = -1;
    end
end
        




       
for i = 1:length(randomizedIDs)
%      if nnz(groupassignmentIDs==randomizedIDs(i)) > 0 && nnz(subassignmentIDs==randomizedIDs(i)) ==0
%          subassignmentIDs(groupassignment(groupassignmentIDs==randomizedIDs(i))) = randomizedIDs(i);
%      end
%      if nnz(subassignmentIDs==randomizedIDs(i)) > 0 && nnz(groupassignmentIDs==randomizedIDs(i)) ==0
%          groupassignmentIDs(subassignment(subassignmentIDs==randomizedIDs(i))) = randomizedIDs(i);
%      end
     
     if nnz(groupassignmentIDs==randomizedIDs(i)) > 1 && nnz(subassignmentIDs==randomizedIDs(i)) > 1
         groupinds = find(groupassignmentIDs==randomizedIDs(i)); subinds = find(subassignmentIDs==randomizedIDs(i));
         [tempassignment ign] = munkres(groupvsubcost(groupinds,subinds));
         [tempgassignment tempsassignment ign ign] = munkres_mult(groupvsubcost(groupinds,subinds),4);
         
         for j = 2:length(tempassignment)
             if tempassignment(j)~=0
                 unusedIDs = setdiff(randomizedIDs,groupassignmentIDs);
                 random_unusedID = unusedIDs(randi(length(unusedIDs),1));
                 groupassignmentIDs(groupinds(j)) = random_unusedID;
                 subassignmentIDs(subinds(tempassignment(j))) = random_unusedID;
                 groupassignment(groupinds(j)) = subinds(tempassignment(j));
                 subassignment(subinds(tempassignment(j))) = groupinds(j);
             end
         end
         
         for j = 1:length(tempassignment)
             if tempassignment(j)==0
                 groupassignment(groupinds(j)) = subinds(tempgassignment(j));
                 groupassignmentIDs(groupinds(j)) = subassignmentIDs(subinds(tempgassignment(j)));
                 
             end
         end
         
         for j = 1:length(subinds)
             if ~any(tempassignment==j)
                 subassignment(subinds(j)) = groupinds(tempsassignment(j));
                 subassignmentIDs(subinds(j)) = groupassignmentIDs(groupinds(tempsassignment(j)));
             end
         end
     end
         
end 





for parcelnum = 1:length(parcelIDs)
    
    if (groupassignment(parcelnum) > 0) %&& (onevtwo_pctoverlap(parcelnum,groupassignment(parcelnum)) > 0)
        
        groupassignmap(parcels==parcelIDs(parcelnum)) = groupassignmentIDs(parcelnum);
        if isempty(find(subassignment==parcelnum))
            groupsimilaritymap(parcels==parcelIDs(parcelnum)) = groupvsub_voxelcountperc(parcelnum,groupassignment(parcelnum));
            groupdistancemap(parcels==parcelIDs(parcelnum)) = group_v_sub_distance(parcelnum,groupassignment(parcelnum));
        else
            groupsimilaritymap(parcels==parcelIDs(parcelnum)) = mean(groupvsub_voxelcountperc(parcelnum,(subassignment==parcelnum)));
            groupdistancemap(parcels==parcelIDs(parcelnum)) = mean(group_v_sub_distance(parcelnum,(subassignment==parcelnum)));
        end
        
    else
        
        groupassignmap(parcels==parcelIDs(parcelnum)) = -1;
        groupsimilaritymap(parcels==parcelIDs(parcelnum)) = -1;
         groupdistancemap(parcels==parcelIDs(parcelnum)) = -1;
    end
end

for parcelnum = 1:length(subparcelIDs)
    
    if subassignment(parcelnum) > 0  %&& (onevtwo_pctoverlap(twoassignment(parcelnum),parcelnum) > 0)
        
        subassignmap(subparcels==subparcelIDs(parcelnum)) = subassignmentIDs(parcelnum);
        
        if isempty(find(groupassignment==parcelnum))
            subsimilaritymap(subparcels==subparcelIDs(parcelnum)) = groupvsub_voxelcountperc(subassignment(parcelnum),parcelnum);
            subdistancemap(subparcels==subparcelIDs(parcelnum)) = group_v_sub_distance(subassignment(parcelnum),parcelnum);
        else
            subsimilaritymap(subparcels==subparcelIDs(parcelnum)) = mean(groupvsub_voxelcountperc((groupassignment==parcelnum),parcelnum));
            subdistancemap(subparcels==subparcelIDs(parcelnum)) = mean(group_v_sub_distance((groupassignment==parcelnum),parcelnum));
        end
        
    else
        
        subassignmap(subparcels==subparcelIDs(parcelnum)) = -1;
        subsimilaritymap(subparcels==subparcelIDs(parcelnum)) = -1;
        subdistancemap(subparcels==subparcelIDs(parcelnum)) = -1;
        
    end
end

    
    
   subassignmap_temp = zeros(32492,1);
   subassignmap_temp(logical(mask)) = subassignmap;
   subassignmap = subassignmap_temp;
   groupassignmap_temp = zeros(32492,1);
   groupassignmap_temp(logical(mask)) = groupassignmap;
   groupassignmap = groupassignmap_temp;

% done = 0;
% iter = 0;
% while done==0
%     iter = iter+1;
%     %disp(['Iteration ' num2str(iter)])
%     done = 1;
% 
%     subassignments = unique(subassignmap); subassignments(subassignments==0) = [];
%     allvals = [1:max(subassignments)];
%     for i = 1:length(subassignments)
%         vertices = find(subassignmap==subassignments(i));
%         neighbors = find(sum((geo_distances(vertices,:)<10),1));
%         neighbors = setdiff(neighbors,vertices);
%         neighborvals = unique(subassignmap(neighbors)); neighborvals(neighborvals==0) = []; 
%         if any(abs(neighborvals-subassignments(i))<15)
%             inside_neighborvalrange = [];
%             for thisneighval = neighborvals(:)'
%                 inside_neighborvalrange = unique([inside_neighborvalrange [thisneighval-15 : thisneighval+30]]);
%             end
%             remainingvals = setdiff(allvals,[subassignments(:)' inside_neighborvalrange]);
%             newval = remainingvals(randi(length(remainingvals),1));
%             subassignmap(subassignmap==subassignments(i)) = newval;
%             groupassignmap(groupassignmap==subassignments(i)) = newval;
%             done = 0;
%         end
%     end
% end


    
    

     
    clear subparcelindices
    
    groupassignoutput(:,s) = groupassignmap;
    groupsimilarityoutput(logical(mask),s) = groupsimilaritymap;
    groupdistanceoutput(logical(mask),s) = groupdistancemap;
    subassignoutput(:,s) = subassignmap;
    subsimilarityoutput(logical(mask),s) = subsimilaritymap;
    subdistanceoutput(logical(mask),s) = subdistancemap;
    

    

    
end







save(gifti(single(groupassignoutput)),[outputfolder 'Group_parcels_assignedtosubs_count_' hem '.func.gii'])

save(gifti(single(groupsimilarityoutput)),[outputfolder 'Group_parcels_assignedtosubs_countsimilarity_' hem '.func.gii'])

%save(gifti(single(groupdistanceoutput)),[outputfolder 'Group_parcels_assignedtosubs_distance_' hem '.func.gii'])

save(gifti(single(subassignoutput)),[outputfolder 'Subject_parcels_assignedtogroup_count_' hem '.func.gii'])

save(gifti(single(subsimilarityoutput)),[outputfolder 'Subject_parcels_assignedtogroup_countsimilarity_' hem '.func.gii'])

%save(gifti(single(subdistanceoutput)),[outputfolder 'Subject_parcels_assignedtogroup_distance_' hem '.func.gii'])
end
    
    
