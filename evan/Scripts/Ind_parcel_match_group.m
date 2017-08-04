clear

make_avg_correlpatterns = 0;

LFRStmasklist = '/data/cn4/evan/RestingState/Ind_variability/Subjects/LFRS_TMASKLIST.txt';

HEMS = {'L','R'};

parcel_correlpatterns_file = '/data/cn4/evan/RestingState/FC_Mapping_120/All_parcel_correlpatterns_gooddata.dtseries.nii';

bufsize=16384;
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;

adjacent_dist = 5;

xdistance = 20;

for hem = 1:length(HEMS)
    thismedialwall = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' HEMS{hem} '.32k_fs_LR.func.gii']);
    thismedialwall = ~thismedialwall.cdata;
    if hem==1
        ncortexLverts = nnz(thismedialwall);
    else
        ncortexRverts = nnz(thismedialwall);
    end
    nsubcortverts = 7285;
end

parcelcount = 0;
for hem = 1:length(HEMS)
    thismedialwall = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' HEMS{hem} '.32k_fs_LR.func.gii']);
    thismedialwall = ~thismedialwall.cdata;
    
    thiswatershed = gifti(['/data/cn4/evan/RestingState/FC_Mapping_120/120_' HEMS{hem} '_wateredgethresh_watershedmerge_0.45_gooddata.func.gii']);
    watershed{hem} = thiswatershed.cdata;
    ciftispace_watershed{hem} = watershed{hem}(logical(thismedialwall));
    
    parcelIDs{hem} = unique(watershed{hem}); parcelIDs{hem}(parcelIDs{hem}==0) = [];
    
    load(['/data/cn4/evan/fsaverage_LR32k/Surface_distances_' HEMS{hem} '.mat']);
    geo_distances_byhem{hem} = geo_distances;
    clear geo_distances
    sphere = gifti(['/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k/Conte69.' HEMS{hem} '.sphere.32k_fs_LR.surf.gii']);
    
    for parcelnum=1:length(parcelIDs{hem})
        parcelcount = parcelcount+1;
        parcelindices{parcelcount} = find(ciftispace_watershed{hem}==parcelIDs{hem}(parcelnum)) + (strcmp(HEMS{hem},'R') * ncortexLverts);
        
        gifti_parcelindices{hem}{parcelnum} = find(watershed{hem}==parcelIDs{hem}(parcelnum));
        
        ind = gifti_parcelindices{hem}{parcelnum};
        meanX = mean(sphere.vertices(ind,1));
        meanY = mean(sphere.vertices(ind,2));
        meanZ = mean(sphere.vertices(ind,3));
        coord = [meanX meanY meanZ];
        sphere_coords = [sphere.vertices(ind,1) sphere.vertices(ind,2) sphere.vertices(ind,3)];
        rep_coord = repmat(coord, [size(sphere_coords,1) 1]);
        dist_coord = sum((sphere_coords-rep_coord).^2,2).^(1/2);
        [y indval] = min(dist_coord);
        centroidvert = ind(indval);
        
        if hem==1
            validcorrelverts{parcelcount} = [(geo_distances_byhem{hem}(logical(thismedialwall),centroidvert) > xdistance); ones(ncortexRverts,1); ones(nsubcortverts,1)];
        else
            validcorrelverts{parcelcount} = [ones(ncortexLverts,1); (geo_distances_byhem{hem}(logical(thismedialwall),centroidvert) > xdistance); ones(nsubcortverts,1)];
        end
        
        
        adjacent_verts = gifti_parcelindices{hem}{parcelnum}(:);
        for iter = 1:adjacent_dist
            vertneighbors = neighbors(adjacent_verts,2:7);
            vertneighbors = vertneighbors(:); vertneighbors(isnan(vertneighbors)) = [];
            adjacent_verts = unique([adjacent_verts; vertneighbors]);
        end
        
        adjacent_parcelIDs = unique(watershed{hem}(adjacent_verts)); adjacent_parcelIDs(adjacent_parcelIDs==0) = []; adjacent_parcelIDs(adjacent_parcelIDs==parcelIDs{hem}(parcelnum)) = [];
        
        for thisadjacentparcel = 1:length(adjacent_parcelIDs)
            adjacent_parcelnums{parcelcount}(thisadjacentparcel) = find(parcelIDs{hem}==adjacent_parcelIDs(thisadjacentparcel)) + (strcmp(HEMS{hem},'R') * length(parcelIDs{1}));
        end
        
        
        
    end
    
end
numparcels = parcelcount;
%% Create group average correlation patterns
if make_avg_correlpatterns
    
    tmaskfile = '/data/cn4/evan/RestingState/FC_Mapping_120/Old_concat/AllC_TMASKLIST.txt';
    [subjects tmasks] = textread(tmaskfile,'%s %s');
    
    mean_parcel_correlpatterns = zeros(66697,numparcels);
    
    for s = 1:length(subjects)
        disp(['Processing subject #' num2str(s) ': ' subjects{s}])
        tmask = load(tmasks{s});
        
        surf_timecourse = gifti(['/data/hcp-bluearc/home/laumannt/120_parcellation/cifti_timeseries_normalwall/' subjects{s} '_BOLD_LR_surf_subcort_32k_fsLR_smooth2.55.func.gii']);
        surf_timecourse = surf_timecourse.cdata;
        surf_timecourse = surf_timecourse(:,logical(tmask));
        surf_timecourse(isnan(surf_timecourse)) = 0;
        
        for parcelnum = 1:numparcels
            
            parcel_timecourse = mean(surf_timecourse(parcelindices{parcelnum},:),1);
            
            parcel_correlpattern = paircorr_mod(surf_timecourse',parcel_timecourse');
            parcel_correlpattern(isnan(parcel_correlpattern)) = 0;
            
            mean_parcel_correlpatterns(:,parcelnum) = mean_parcel_correlpatterns(:,parcelnum) + (parcel_correlpattern ./ length(subjects));
            
            clear parcel_correlpattern
        end
    end
    
    ciftitemplatefile = ['/data/hcp-bluearc/home/laumannt/120_parcellation/cifti_timeseries_normalwall/vc25125_BOLD_LR_surf_subcort_32k_fsLR_smooth2.55.func.gii'];
    
    cifti_write_wHDR(mean_parcel_correlpatterns,ciftitemplatefile,parcel_correlpatterns_file(1:end-13))
    
else
    mean_parcel_correlpatterns = cifti_read(parcel_correlpatterns_file);
end

%%

[LFRSsubjects LFRStmasks] = textread(LFRStmasklist,'%s%s');

for s = 1:length(LFRSsubjects)
    disp(['Running LFRS subject #' num2str(s) ': ' LFRSsubjects{s}])
    
    tmask = load(LFRStmasks{s});
    
    surf_timecourse = cifti_read(['/data/hcp-bluearc/home/laumannt/LFRS_parcellation/' LFRSsubjects{s} '/' LFRSsubjects{s} '_timeseries_normalwall.dtseries.nii']);
    
    surf_timecourse = surf_timecourse(:,logical(tmask));
    surf_timecourse(isnan(surf_timecourse)) = 0;
    
    groupparcel_subcorrelpatterns = zeros(66697,numparcels);
    
    for parcelnum = 1:numparcels
        parcel_timecourse = mean(surf_timecourse(parcelindices{parcelnum},:),1);
        
        parcel_correlpattern = paircorr_mod(surf_timecourse',parcel_timecourse');
        parcel_correlpattern(isnan(parcel_correlpattern)) = 0;
        
        groupparcel_subcorrelpatterns(:,parcelnum) = parcel_correlpattern;
        
        clear parcel_correlpattern
        
        grouptogroup_correls(s,parcelnum) = paircorr_mod(groupparcel_subcorrelpatterns(logical(validcorrelverts{parcelnum}),parcelnum),mean_parcel_correlpatterns(logical(validcorrelverts{parcelnum}),parcelnum));
    end
    
    for parcelnum = 1:numparcels
        for adjacent_parcelnum = 1:length(adjacent_parcelnums{parcelnum})
            verts_tocorrelate = logical(validcorrelverts{parcelnum}.*validcorrelverts{adjacent_parcelnums{parcelnum}(adjacent_parcelnum)});
            adjacent_correls(adjacent_parcelnum) = paircorr_mod(groupparcel_subcorrelpatterns(verts_tocorrelate,adjacent_parcelnums{parcelnum}(adjacent_parcelnum)),mean_parcel_correlpatterns(verts_tocorrelate,parcelnum));
        end
        grouptoadjacent_correls(s,parcelnum) = mean(adjacent_correls);
        clear adjacent_correls
    end
        
    
    
    
    for hem = 1:length(HEMS)
        
        load(['/data/cn4/evan/fsaverage_LR32k/Surface_distances_' HEMS{hem} '.mat']);
        sphere = gifti(['/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k/Conte69.' HEMS{hem} '.sphere.32k_fs_LR.surf.gii']);
        
        thismedialwall = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' HEMS{hem} '.32k_fs_LR.func.gii']);
        thismedialwall = ~thismedialwall.cdata;
        
        g2sassignments = gifti(['/data/cn4/evan/RestingState/Ind_variability/Subjects/Assignments/Matching_expand/Group_parcels_matchtosubs_' HEMS{hem} '.func.gii']);
        g2sassignments = g2sassignments.cdata;
        g2sassignments_cifti = g2sassignments(logical(thismedialwall),:);
        s2gassignments = gifti(['/data/cn4/evan/RestingState/Ind_variability/Subjects/Assignments/Matching_expand/Subject_parcels_matchtogroup_' HEMS{hem} '.func.gii']);
        s2gassignments = s2gassignments.cdata;
        s2gassignments_cifti = s2gassignments(logical(thismedialwall),:);
        
        for hem_parcelnum = 1:length(parcelIDs{hem})
            parcelnum = hem_parcelnum + (strcmp(HEMS{hem},'R') * length(parcelIDs{1}));
            matchID = mean(g2sassignments(gifti_parcelindices{hem}{hem_parcelnum},s));
            assignedindices = find(s2gassignments_cifti(:,s)==matchID) + (strcmp(HEMS{hem},'R') * ncortexLverts);
            
            if isempty(assignedindices)
                grouptoind_correls(s,parcelnum) = NaN;
            end
            
            parcel_timecourse = mean(surf_timecourse(assignedindices,:),1);
        
            parcel_correlpattern = paircorr_mod(surf_timecourse',parcel_timecourse');
            parcel_correlpattern(isnan(parcel_correlpattern)) = 0;
            
            gifti_assignedindices = find(s2gassignments(:,s)==matchID);
            ind = gifti_assignedindices;
            meanX = mean(sphere.vertices(ind,1));
            meanY = mean(sphere.vertices(ind,2));
            meanZ = mean(sphere.vertices(ind,3));
            coord = [meanX meanY meanZ];
            sphere_coords = [sphere.vertices(:,1) sphere.vertices(:,2) sphere.vertices(:,3)];
            rep_coord = repmat(coord, [size(sphere_coords,1) 1]);
            dist_coord = sum((sphere_coords-rep_coord).^2,2).^(1/2);
            [y indval] = min(dist_coord);
            centroidvert = indval;
            
            if hem==1
                assignedvalidcorrelverts = [(geo_distances_byhem{hem}(logical(thismedialwall),centroidvert) > 20); ones(ncortexRverts,1); ones(nsubcortverts,1)];
            else
                assignedvalidcorrelverts = [ones(ncortexLverts,1); (geo_distances_byhem{hem}(logical(thismedialwall),centroidvert) > 20); ones(nsubcortverts,1)];
            end
            
            verts_tocorrelate = validcorrelverts{parcelnum}.*assignedvalidcorrelverts;
            
            grouptoind_correls(s,parcelnum) = paircorr_mod(parcel_correlpattern(logical(verts_tocorrelate)),mean_parcel_correlpatterns(logical(verts_tocorrelate),parcelnum));
        end
    end
end

%%

figure
centers = [-.5: .1:.9];
[n,x] = hist(grouptoadjacent_correls(:),centers);
plot(centers,n/numel(grouptoadjacent_correls),'g');

hold on
[n x] = hist(grouptogroup_correls(:),centers);
plot(centers,n/numel(grouptogroup_correls),'b-');

[n x] = hist(grouptoind_correls(logical(~isnan(grouptoind_correls))),centers);
plot(centers,n./nnz(~isnan(grouptoind_correls)),'r-');
xlim([-.5 1]); ylim([0 .4])

[H,P,CI,STATS] = ttest2(FisherTransform(grouptoind_correls(logical(~isnan(grouptoind_correls)))),FisherTransform(grouptogroup_correls(:)))
    
    
