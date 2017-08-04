%function Match_group_parcels_to_ind_data(parcelfilename,allsubs_correlmatfile,datalist,tmasklist,hem,iscifti)
% Match_group_parcels_to_ind_data(parcelfilename,allsubs_correlmat,datalist,tmasklist,hem,iscifti)

parcelfilename = '120_R_crossthresh_watershedmerge.func.gii';
%allsubs_timecourse = 'allsubs_LR_timeseries.func.gii';
%allsubs_tmask = 'allsubs_total_tmask.txt';
allsubs_correlmatfile = 'allsubs_concat_corr_R.func.gii';
datalist = 'AllC_DATALIST.txt';
tmasklist = 'AllC_TMASKLIST.txt';
hem = 'R';
iscifti = 2;
max_dist = 50;


origparcels = gifti(parcelfilename); origparcels = origparcels.cdata;



%%

if iscifti==1
    maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii']; mask = gifti(maskname); mask = ~mask.cdata;
    parcels = origparcels(logical(mask));
    maskL = gifti('/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii');ncortexLverts = nnz(maskL.cdata==0);
elseif iscifti==2
    maskname = ['/data/cn4/laumannt/subcortical_mask/' hem '.atlasroi_erode3.32k_fs_LR.shape.gii']; mask = gifti(maskname); mask = mask.cdata;
    parcels = origparcels(logical(mask));
    maskL = gifti('/data/cn4/laumannt/subcortical_mask/L.atlasroi_erode3.32k_fs_LR.shape.gii');ncortexLverts = nnz(maskL.cdata);
else
    parcels = origparcels;
    ncortexLverts = length(parcels);
end

disp('Loading overall correlation matrix')

% alltimecourse = gifti(allsubs_timecourse); alltimecourse = alltimecourse.cdata;
% if ~isempty(allsubs_tmask)
%     alltmask = load(allsubs_tmask);
%     alltimecourse = alltimecourse(:,logical(alltmask));
% end

if ischar(allsubs_correlmatfile)
    allsubs_correlmat = gifti(allsubs_correlmatfile); allsubs_correlmat = allsubs_correlmat.cdata;
else
    allsubs_correlmat = allsubs_correlmatfile;
    clear allsubs_correlmatfile
end

parcelIDs = unique(parcels); parcelIDs(parcelIDs==0) = [];

groupparcelpatterns = zeros(length(parcelIDs),size(allsubs_correlmat,2));

for parcelnum = 1:length(parcelIDs)
    
    string{parcelnum} = ['Calculating group connectivity patterns for parcel number ' num2str(parcelnum) ' out of ' num2str(length(parcelIDs))];
    if parcelnum==1; fprintf('%s',string{parcelnum}); else fprintf([repmat('\b',1,length(string{parcelnum-1})) '%s'],string{parcelnum}); end
    
    parcelindices{parcelnum} = find(parcels==parcelIDs(parcelnum));
%     if strcmp(hem,'R')
%         rightparcelindices{parcelnum} = parcelindices{parcelnum}+ncortexLverts;
%     end
    %groupparcelpatterns(parcelnum,:) = paircorr_mod(mean(alltimecourse(parcelindices{parcelnum},:),1)',alltimecourse');
    groupparcelpatterns(parcelnum,:) = mean(allsubs_correlmat(parcelindices{parcelnum},:),1);
end

clear allsubs_correlmat

disp(' ')

%%

surfdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
sphere = gifti([surfdir '/Conte69.' hem '.sphere.32k_fs_LR.surf.gii']);

[phi theta r] = cart2sph(sphere.vertices(:,1), sphere.vertices(:,2),sphere.vertices(:,3));

for parcelnum = 1:length(parcelIDs)
    
    ind = find(origparcels==parcelIDs(parcelnum));
    
    meanX = mean(sphere.vertices(ind,1));
    meanY = mean(sphere.vertices(ind,2));
    meanZ = mean(sphere.vertices(ind,3));
    
    coord = [meanX meanY meanZ];
    sphere_coords = [sphere.vertices(ind,1) sphere.vertices(ind,2) sphere.vertices(ind,3)];
    
    rep_coord = repmat(coord, [size(sphere_coords,1) 1]);
    
    dist_coord = sum((sphere_coords-rep_coord).^2,2).^(1/2);
    [y indval] = min(dist_coord);
    indpos(parcelnum) = ind(indval);
    
end

metric = zeros(32492,1);
metric(indpos) = 1;

% medial_wall = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii']); medial_wall = medial_wall.cdata;
% metric(logical(medial_wall)) = 0;
 indpos = find(metric);

%Save out midthickness coordinates of centers
midthick = gifti([surfdir '/Conte69.' hem '.midthickness.32k_fs_LR.surf.gii']);
center_coords = [midthick.vertices(indpos,1) midthick.vertices(indpos,2) midthick.vertices(indpos,3)];
quickroifile(center_coords,'Temp_parcel_center.roi');
dist_mat = euclidean_distance('Temp_parcel_center.roi');
parcel_comparisons_to_keep = (dist_mat < max_dist);

%%


[subjects subdata] = textread(datalist,'%s%s');
if ~isempty(tmasklist)
    [subjects tmasks] = textread(tmasklist,'%s%s');
end

match_outputfile = zeros(length(mask),1);
similarity_outputfile = zeros(length(mask),length(subjects));

sub_match_count = zeros(length(subjects),1);

for subnum = 1:length(subjects)
    
    disp(['Subject ' num2str(subnum) ':'])
    disp('   Loading timecourse')
    
    evalc(['!wb_command -cifti-convert -to-gifti-ext ' subdata{subnum} ' Temp.func.gii']);
    subtimecourse = gifti('Temp.func.gii'); subtimecourse = subtimecourse.cdata;
    if ~isempty(tmasklist)
        tmask = load(tmasks{subnum});
        subtimecourse = subtimecourse(:,logical(tmask));
    end
    
    
    for parcelnum = 1:length(parcelIDs)
        
        string{parcelnum} = ['   Calculating subject connectivity patterns for parcel number ' num2str(parcelnum) ' out of ' num2str(length(parcelIDs))];
        if parcelnum==1; fprintf('%s',string{parcelnum}); else fprintf([repmat('\b',1,length(string{parcelnum-1})) '%s'],string{parcelnum}); end
    
        if strcmp(hem,'R')
            thisparcelpattern = paircorr_mod(mean(subtimecourse((parcelindices{parcelnum}+ncortexLverts),:),1)',subtimecourse');
        else
            thisparcelpattern = paircorr_mod(mean(subtimecourse(parcelindices{parcelnum},:),1)',subtimecourse');
        end
        notnanindices = logical(~isnan(thisparcelpattern));
        thisparcelsimilarity = paircorr_mod(thisparcelpattern(notnanindices)',groupparcelpatterns(:,notnanindices)');
        similarity_outputfile(origparcels==parcelIDs(parcelnum),subnum) = thisparcelsimilarity(parcelnum);
        [maxval maxi] = max(thisparcelsimilarity .* parcel_comparisons_to_keep(parcelnum,:));
        if any(maxi==parcelnum);
            match_outputfile(origparcels==parcelIDs(parcelnum)) = mean(match_outputfile(origparcels==parcelIDs(parcelnum))) + (1/length(subjects));
            sub_match_count(subnum) = sub_match_count(subnum)+1;
        end
    end
    disp(' ')
    
    disp(['   ' num2str(sub_match_count(subnum)) ' of ' num2str(length(parcelIDs)) ' parcels match'])
    
end

save(gifti(single(match_outputfile)),['Subject_matches2_to_' parcelfilename]);
save(gifti(single(mean(similarity_outputfile,2))),['Subject_avgcorrelation2_to_' parcelfilename]);


    
    
    

