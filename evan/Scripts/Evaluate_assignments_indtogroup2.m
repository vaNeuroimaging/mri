datalist = '/data/cn4/evan/RestingState/Ind_variability/Subjects/LFRS_DATALIST.txt';
tmasklist = '/data/cn4/evan/RestingState/Ind_variability/Subjects/LFRS_TMASKLIST.txt';

outputfolder = '/data/cn4/evan/RestingState/Ind_variability/Subjects/Assignments/Matching_expand/';

distanceexclusion = 20;

hems = {'L','R'};

subjectnums = [];

centroid_distance_map_outputname = ['Assignment_median_centroid_distance2'];

%-------------------------------------------------------------------------

templatefile = '/data/hcp-zfs/home/laumannt/120_parcellation/cifti_timeseries_normalwall/vc25125_BOLD_LR_surf_subcort_32k_fsLR_smooth2.55.func.gii';
templatedata = gifti(templatefile);
bufsize=16384;
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;



[subjects subdata] = textread(datalist,'%s%s');
if ~isempty(tmasklist)
    [subjects tmasks] = textread(tmasklist,'%s%s');
end

if ~isempty(subjectnums)
    subjects = subjects(subjectnums); subdata = subdata(subjectnums); tmasks = tmasks(subjectnums);
end

groupparcelindices = [];
subparcelindices = [];

groupparcel_centroid = [];
subparcel_centroid = [];
centroid_distances = [];
centroid_map = [];
allcentroid_map = [];
centroid_of_centroids = [];
fociXYZ = [];
fociRGB = [];
fociname = [];
classname = [];

match = [];

counter = 1;

geo_distances_both = ones(32492*2) * 200;

for hemnum = 1:length(hems)
    hem = hems{hemnum};
    
    
    
    maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii']; mask = gifti(maskname); mask = ~mask.cdata;
    maskL = gifti('/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii');ncortexLverts = nnz(maskL.cdata==0);
    %maskR = gifti('/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.R.32k_fs_LR.func.gii');ncortexRverts = nnz(maskR.cdata==0);

    
    load(['/data/cn4/evan/fsaverage_LR32k/Surface_distances_' hem '.mat']);
    inds = [1:size(geo_distances)] + (size(geo_distances)*strcmp(hem,'R'));
    geo_distances_both(inds,inds) = geo_distances;
    
    sphere = gifti(['/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k/Conte69.' hem '.sphere.32k_fs_LR.surf.gii']);
    
    groupparcels = gifti(['/data/cn4/evan/RestingState/FC_Mapping_120/120_' hem '_wateredgethresh_watershedmerge_0.45_gooddata.func.gii']); 
    groupparcels = groupparcels.cdata;
    groupparcels_cifti = groupparcels(logical(mask));
    parcelIDs = unique(groupparcels); parcelIDs(parcelIDs==0) = [];
%     
%     minsize = 15;
% 
%     gooddataname = ['/data/cn4/evan/ROIs/all_meanimage_' hem '.func.gii'];
%     gooddata = gifti(gooddataname); gooddata = gooddata.cdata;
%     gooddata = gooddata>750;
    
    
    
    g2sassignments = gifti(['/data/cn4/evan/RestingState/Ind_variability/Subjects/Assignments/Matching_expand/Group_parcels_matchtosubs2_' hem '.func.gii']);
    g2sassignments = g2sassignments.cdata;
    g2sassignments_cifti = g2sassignments(logical(mask),:);
    s2gassignments = gifti(['/data/cn4/evan/RestingState/Ind_variability/Subjects/Assignments/Matching_expand/Subject_parcels_matchtogroup2_' hem '.func.gii']);
    s2gassignments = s2gassignments.cdata;
    s2gassignments_cifti = s2gassignments(logical(mask),:);
    
    thishem_centroiddistance_map = zeros(32492,1);
    thishem_centroiddistance_map_byparcel = zeros(32492,length(parcelIDs));
    thishem_allcentroid_map = zeros(32492,length(parcelIDs));
    thishem_match = ones(length(parcelIDs),length(subjects));
    
    
    thishem_fociname_group = [];
    thishem_fociname_mean = [];
    thishem_fociname_ind = [];
    thishem_indcentroids = [];
    for parcelnum = 1:length(parcelIDs)
        
        
        
        string{counter} = [hem 'hem, parcel number ' num2str(parcelnum)];
        if counter==1; fprintf('%s',string{counter}); else fprintf([repmat('\b',1,length(string{counter-1})) '%s'],string{counter}); end
        counter = counter+1;
        
        ind = find(groupparcels==parcelIDs(parcelnum));
        meanX = mean(sphere.vertices(ind,1));
        meanY = mean(sphere.vertices(ind,2));
        meanZ = mean(sphere.vertices(ind,3));
        coord = [meanX meanY meanZ];
        sphere_coords = [sphere.vertices(ind,1) sphere.vertices(ind,2) sphere.vertices(ind,3)];
        rep_coord = repmat(coord, [size(sphere_coords,1) 1]);
        dist_coord = sum((sphere_coords-rep_coord).^2,2).^(1/2);
        [y indval] = min(dist_coord);
        thishem_groupparcel_centroid(parcelnum,1) = ind(indval);
        
        thishem_groupparcelindices{parcelnum,1} = (find(groupparcels_cifti==parcelIDs(parcelnum))) + (strcmp('R',hem) * ncortexLverts);
        for s = 1:length(subjects)
            parcelmatchID = mean(g2sassignments_cifti(logical(groupparcels_cifti==parcelIDs(parcelnum)),s));
            
            ind = find(s2gassignments(:,s)==parcelmatchID);
            
            if isempty(ind)
                thishem_subparcel_centroid(parcelnum,s) = NaN;
                thishem_centroid_distances(parcelnum,s) = NaN;
                %thishem_match(parcelnum,s) = 0;
            else
                meanX = mean(sphere.vertices(ind,1));
                meanY = mean(sphere.vertices(ind,2));
                meanZ = mean(sphere.vertices(ind,3));
                coord = [meanX meanY meanZ];
                %sphere_coords = [sphere.vertices(ind,1) sphere.vertices(ind,2) sphere.vertices(ind,3)];
                rep_coord = repmat(coord, [size(sphere.vertices,1) 1]);
                dist_coord = sum((sphere.vertices-rep_coord).^2,2).^(1/2);
                [y indval] = min(dist_coord);
                thishem_subparcel_centroid(parcelnum,s) = indval;
                    
                thishem_centroid_distances(parcelnum,s) =  geo_distances(thishem_groupparcel_centroid(parcelnum,1),thishem_subparcel_centroid(parcelnum,s));

            end
            
            thishem_subparcelindices{parcelnum,s} = find(s2gassignments_cifti(:,s)==parcelmatchID) + (strcmp('R',hem) * ncortexLverts);
        end
        
        ind = thishem_subparcel_centroid(parcelnum,:);
        ind(isnan(ind)) = [];
        meanX = mean(sphere.vertices(ind,1));
        meanY = mean(sphere.vertices(ind,2));
        meanZ = mean(sphere.vertices(ind,3));
        coord = [meanX meanY meanZ];
        sphere_coords = [sphere.vertices(:,1) sphere.vertices(:,2) sphere.vertices(:,3)];
        rep_coord = repmat(coord, [size(sphere_coords,1) 1]);
        dist_coord = sum((sphere_coords-rep_coord).^2,2).^(1/2);
        [y indval] = min(dist_coord);
        thishem_centroid_of_centroids(parcelnum,1) = indval;
        if thishem_centroid_of_centroids(parcelnum,1)==thishem_groupparcel_centroid(parcelnum,1)
            thishem_centroid_of_centroids(parcelnum,1) = neighbors(thishem_centroid_of_centroids(parcelnum,1),2);
        end
        
        mean_distance = nanmean(thishem_centroid_distances(parcelnum,:));

        verts_within_mean = find(geo_distances(thishem_groupparcel_centroid(parcelnum,1),:)<=mean_distance);
        for vert = verts_within_mean(:)'
            vertneighs = neighbors(vert,2:7); vertneighs(isnan(vertneighs)) = [];
            if ~isempty(setdiff(vertneighs,verts_within_mean));
                thishem_centroiddistance_map(vert) = .5;
                thishem_centroiddistance_map_byparcel(vert,parcelnum) = .5;
                %thishem_allcentroid_map(vert,parcelnum) = .5;
            end
        end
        
        %thishem_allcentroid_map(thishem_groupparcel_centroid(parcelnum,1),parcelnum) = 3;
        %thishem_allcentroid_map(thishem_centroid_of_centroids(parcelnum),parcelnum) = 1;
        ind = thishem_subparcel_centroid(parcelnum,:);
        ind(isnan(ind)) = [];
        %thishem_allcentroid_map(ind,parcelnum) = 2;
        thishem_fociname_group{end+1,1} = [hem 'Parcel' num2str(parcelnum) '_parcelID' sprintf('%03i',parcelIDs(parcelnum)) '_group'];
        thishem_fociname_mean{end+1,1} = [hem 'Parcel' num2str(parcelnum) '_parcelID' sprintf('%03i',parcelIDs(parcelnum)) '_mean'];
        thishem_indcentroids = [thishem_indcentroids; ind'];
        for i = 1:length(ind)
            thishem_fociname_ind{end+1,1} = [hem 'Parcel' num2str(parcelnum) '_parcelID' sprintf('%03i',parcelIDs(parcelnum)) '_sub' num2str(i)];
        end
    end
    
    
    
    groupparcelindices = [groupparcelindices ; thishem_groupparcelindices];
    subparcelindices = [subparcelindices ; thishem_subparcelindices];
    
    groupparcel_centroid = [groupparcel_centroid ; thishem_groupparcel_centroid+(32492*strcmp(hem,'R'))];
    subparcel_centroid = [subparcel_centroid ; thishem_subparcel_centroid+(32492*strcmp(hem,'R'))];
    centroid_distances = [centroid_distances ; thishem_centroid_distances];
    
    %thishem_centroid_map(thishem_groupparcel_centroid(:,1)) = 2;
    %thishem_centroid_map(thishem_centroid_of_centroids) = 1;
    centroid_map = [centroid_map ; thishem_centroiddistance_map(logical(mask))];
    
%     thishem_allcentroidmap_output = zeros(size(templatedata.cdata,1),length(parcelIDs));
%     thishem_allcentroidmap_output((1:nnz(mask))+(ncortexLverts * strcmp(hem,'R')),:) = thishem_allcentroid_map(logical(mask),:);
%     allcentroid_map = [allcentroid_map thishem_allcentroidmap_output];
    
    fociRGB = [repmat([1 1 0],length(thishem_groupparcel_centroid),1) ; repmat([1 0 0],length(thishem_groupparcel_centroid),1) ; repmat([1 .5 0],length(thishem_indcentroids),1)];
    fociname = [thishem_fociname_group; thishem_fociname_mean; thishem_fociname_ind];
    classname = [repmat({'Group parcel centroid'},length(thishem_groupparcel_centroid),1); repmat({'Individual mean centroid'},length(thishem_groupparcel_centroid),1); repmat({'Individual parcel centroid'},length(thishem_indcentroids),1)];
    fociverts = [thishem_groupparcel_centroid; thishem_centroid_of_centroids; thishem_indcentroids];
    
    focifilemaker_workbench(fociverts,fociRGB,fociname,classname,[outputfolder '/' centroid_distance_map_outputname '_' hem '.foci'],hem);
    
    save(gifti(single(thishem_centroiddistance_map_byparcel)),[centroid_distance_map_outputname '_byparcels_' hem '.func.gii'])
    
    %match = [match ; thishem_match];
    
    clear thishem_groupparcelindices thishem_fociname_mean thishem_indcentroids thishem_subparcelindices thishem_groupparcel_centroid thishem_subparcel_centroid thishem_centroid_distances thishem_centroid_map thishem_nomatch thishem_centroid_of_centroids thishem_allcentroidmap_output thishem_allcentroid_map
    
end


centroid_map_output = zeros(size(templatedata.cdata,1),1);
centroid_map_output(1:size(centroid_map,1)) = centroid_map;
cd(outputfolder)
cifti_write_wHDR(centroid_map_output,templatefile,centroid_distance_map_outputname);

%cifti_write_wHDR(allcentroid_map,templatefile,[centroid_distance_map_outputname '_all']);

figure
hist(centroid_distances(:))

%%
groupparcel_correlmat = zeros(size(groupparcelindices,1),size(groupparcelindices,1),length(subjects));
subparcel_correlmat = zeros(size(groupparcelindices,1),size(groupparcelindices,1),length(subjects));
matchmat = zeros(size(groupparcelindices,1),size(groupparcelindices,1),length(subjects));

maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii']; mask = gifti(maskname); maskL = ~mask.cdata;
Power_L = gifti('/data/cn4/evan/ROIs/264_surfvert_ROIs_L_gooddata.func.gii'); Power_L = Power_L.cdata;
Power_L = Power_L(logical(maskL));
parcelIDs = unique(Power_L); parcelIDs(parcelIDs==0) = [];
powerparcelindices = [];
for i = 1:length(parcelIDs)
    powerparcelindices{end+1} = find(Power_L==parcelIDs(i));
end
maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.R.32k_fs_LR.func.gii']; mask = gifti(maskname); maskR = ~mask.cdata;
Power_R = gifti('/data/cn4/evan/ROIs/264_surfvert_ROIs_R_gooddata.func.gii'); Power_R = Power_R.cdata;
Power_R = Power_R(logical(maskR));
parcelIDs = unique(Power_R); parcelIDs(parcelIDs==0) = [];

for i = 1:length(parcelIDs)
    powerparcelindices{end+1} = find(Power_R==parcelIDs(i)) + ncortexLverts;
end

for s = 1:length(subjects)
    
    evalc(['!wb_command -cifti-convert -to-gifti-ext /data/hcp-zfs/home/laumannt/LFRS_parcellation/' subjects{s} '/' subjects{s} '_timeseries_normalwall.dtseries.nii Temp.func.gii']);
    subtimecourse = gifti('Temp.func.gii'); subtimecourse = subtimecourse.cdata;
    delete('Temp.func.gii*');
    if ~isempty(tmasklist)
        tmask = load(tmasks{s});
        subtimecourse = subtimecourse(:,logical(tmask));
    end
    subtimecourse(isnan(subtimecourse)) = 0;
    
    groupparceltimecourses = zeros(size(groupparcelindices,1),size(subtimecourse,2));
    subparceltimecourses = zeros(size(groupparcelindices,1),size(subtimecourse,2));
    
    for parcelnum = 1:size(groupparcelindices,1)
        
        groupparceltimecourses(parcelnum,:) = mean(subtimecourse(groupparcelindices{parcelnum},:),1);
        
        if ~isempty(subparcelindices{parcelnum,s})
        
            subparceltimecourses(parcelnum,:) = mean(subtimecourse(subparcelindices{parcelnum,s},:),1);
        else
            subparceltimecourses(parcelnum,:) = NaN;
        end
        
    end
    
    groupparcel_correlmat(:,:,s) = FisherTransform(paircorr_mod(groupparceltimecourses'));
    subparcel_correlmat(:,:,s) = FisherTransform(paircorr_mod(subparceltimecourses'));
    %matchmat(:,:,s) = match(:,s) * match(:,s)';
    
    clear groupparceltimecourses subparceltimecourses
    
    powerparceltimecourses = zeros(length(powerparcelindices),size(subtimecourse,2));
    for parcelnum = 1:length(powerparcelindices)
        
        powerparceltimecourses(parcelnum,:) = mean(subtimecourse(powerparcelindices{parcelnum},:),1);
    end
    powerparcel_correlmat(:,:,s) = FisherTransform(paircorr_mod(powerparceltimecourses'));
    clear powerparceltimecourses
end
%%
groupcorrelvals = [];
subcorrelvals = [];
powercorrelvals = [];

for subi = 1:length(subjects)
    for subj = 1:length(subjects)
        if subi > subj
            
            subimat = subparcel_correlmat(:,:,subi); subimat = subimat(:);
            subjmat = subparcel_correlmat(:,:,subj);  subjmat = subjmat(:);
            parcelsmatched = logical((~isnan(subimat)) .* (~isnan(subjmat)));
            notnan_centroidsi = subparcel_centroid(:,subi); notnan_centroidsi(isnan(notnan_centroidsi)) = 1;
            notnan_centroidsj = subparcel_centroid(:,subj); notnan_centroidsj(isnan(notnan_centroidsj)) = 1;
            xdistmati = geo_distances(notnan_centroidsi) > distanceexclusion;
            xdistmatj = geo_distances(notnan_centroidsj) > distanceexclusion;
            notnans = parcelsmatched .* xdistmati(:) .* xdistmatj(:);
            subcorrelvals(end+1) = corr2(subimat(notnans),subjmat(notnans));
            sub_subtosubcorrelmat(subi,subj) = corr2(subimat(notnans),subjmat(notnans));
            

            subimat = groupparcel_correlmat(:,:,subi); subimat = subimat(:);
            subjmat = groupparcel_correlmat(:,:,subj); subjmat = subjmat(:);
            xdistmati = geo_distances(groupparcel_centroid) > distanceexclusion;
            %xdistmatj = xdistmati;
            notnans = parcelsmatched .* xdistmati(:);
            groupcorrelvals(end+1) = corr2(subimat(notnans),subjmat(notnans));
            group_subtosubcorrelmat(subi,subj) = corr2(subimat(notnans),subjmat(notnans));
            
%             subimat = powerparcel_correlmat(:,:,subi); subimat = subimat(:);
%             subjmat = powerparcel_correlmat(:,:,subj); subjmat = subjmat(:);
%             
%             powercorrelvals(end+1) = corr2(subimat,subjmat);
%             power_subtosubcorrelmat(subi,subj) = corr2(subimat,subjmat);
        end
    end
end


mean_subparcel_correlmat = nanmean(subparcel_correlmat,3);
mean_groupparcel_correlmat = mean(groupparcel_correlmat,3);
%mean_power_correlmat = mean(powerparcel_correlmat,3);
%%
figure
parcel_correlmat_figmaker(mean_groupparcel_correlmat,'/data/cn4/evan/RestingState/FC_Mapping_120/120_LR_infomap/consensus_assignments_sephems_gooddata.txt',[-.6 .6],'Group')
figure
parcel_correlmat_figmaker(mean_subparcel_correlmat,'/data/cn4/evan/RestingState/FC_Mapping_120/120_LR_infomap/consensus_assignments_sephems_gooddata.txt',[-.6 .6],'Subject')

mean(groupcorrelvals)
figure
imagesc(group_subtosubcorrelmat,[.6,.85])
colormap hot
colorbar

mean(subcorrelvals)
figure
imagesc(sub_subtosubcorrelmat,[.6,.85])
colormap hot
colorbar

% mean(powercorrelvals)
% figure
% imagesc(power_subtosubcorrelmat,[.4,1])
% colormap hot
% colorbar
% 
% disp('Group vs Power')
% [H,P,CI,STATS] = ttest(FisherTransform(groupcorrelvals),FisherTransform(powercorrelvals))

disp('Sub vs Group')
[H,P,CI,STATS] = ttest(FisherTransform(subcorrelvals),FisherTransform(groupcorrelvals))
    
    
   
        
    
    