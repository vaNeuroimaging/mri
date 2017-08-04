hem = 'L';

subjectnums = [];

datalist = '/data/cn4/evan/RestingState/Ind_variability/Subjects/LFRS_DATALIST.txt';
tmasklist = '/data/cn4/evan/RestingState/Ind_variability/Subjects/LFRS_TMASKLIST.txt';
iscifti = 2;

outputfolder = '/data/cn4/evan/RestingState/Ind_variability/Subjects/'; %location match data will be written to

parcelfilename = ['/data/cn4/evan/RestingState/FC_Mapping_120/120_' hem '_watershedmerge.func.gii'];
parcels = gifti(parcelfilename); parcels = parcels.cdata;
origparcels = parcels;

grouptosubs_assignments_filename = ['/data/cn4/evan/RestingState/Ind_variability/Subjects/Target250/Group_parcels_assignedtosubs_' hem '.func.gii'];
grouptosubs_assignments = gifti(grouptosubs_assignments_filename); grouptosubs_assignments = grouptosubs_assignments.cdata;

substogroup_assignments_filename = ['/data/cn4/evan/RestingState/Ind_variability/Subjects/Target250/Subject_parcels_assignedtogroup_' hem '.func.gii'];
substogroup_assignments = gifti(substogroup_assignments_filename); substogroup_assignments = substogroup_assignments.cdata;

parcelIDs = unique(parcels);
parcelIDs(parcelIDs==0)=[];



if iscifti==1
    maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii']; mask = gifti(maskname); mask = ~mask.cdata;
    parcels = parcels(logical(mask));
    substogroup_assignments = substogroup_assignments(logical(mask),:);
    grouptosubs_assignments = grouptosubs_assignments(logical(mask),:);
    maskL = gifti('/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii');ncortexLverts = nnz(maskL.cdata==0);
elseif iscifti==2
    maskname = ['/data/cn4/laumannt/subcortical_mask/' hem '.atlasroi_erode3.32k_fs_LR.shape.gii']; mask = gifti(maskname); mask = mask.cdata;
    parcels = parcels(logical(mask));
    substogroup_assignments = substogroup_assignments(logical(mask),:);
    grouptosubs_assignments = grouptosubs_assignments(logical(mask),:);
    maskL = gifti('/data/cn4/laumannt/subcortical_mask/L.atlasroi_erode3.32k_fs_LR.shape.gii');ncortexLverts = nnz(maskL.cdata);
    maskR = gifti('/data/cn4/laumannt/subcortical_mask/R.atlasroi_erode3.32k_fs_LR.shape.gii');ncortexRverts = nnz(maskR.cdata);
else
    parcels = parcels_upsampled;
    ncortexLverts = length(parcels);
end



[subjects subdata] = textread(datalist,'%s%s');
if ~isempty(tmasklist)
    [subjects tmasks] = textread(tmasklist,'%s%s');
end

if ~isempty(subjectnums)
    subjects = subjects(subjectnums); subdata = subdata(subjectnums); tmasks = tmasks(subjectnums);
end


sub_correlation_patterns = zeros(nnz(mask),length(parcelIDs),length(subjects));
group_correlation_patterns = zeros(nnz(mask),length(parcelIDs),length(subjects));


for s = 1:length(subjects)
    
    disp(['Subject ' subjects{s}])
    
    evalc(['!wb_command -cifti-convert -to-gifti-ext ' subdata{s} ' Temp.func.gii']);
    subtimecourse = gifti('Temp.func.gii'); subtimecourse = subtimecourse.cdata;
    if ~isempty(tmasklist)
        tmask = load(tmasks{s});
        subtimecourse = subtimecourse(:,logical(tmask));
    end
    subtimecourse(isnan(subtimecourse)) = 0;
    
    subparceltimecourses = zeros(length(parcelIDs),size(subtimecourse,2));
    groupparceltimecourses = zeros(length(parcelIDs),size(subtimecourse,2));
    
    for parcelnum = 1:length(parcelIDs)
        string{parcelnum} = ['  Parcel number ' num2str(parcelnum) ' out of ' num2str(length(parcelIDs))];
        if parcelnum==1; fprintf('%s',string{parcelnum}); else fprintf([repmat('\b',1,length(string{parcelnum-1})) '%s'],string{parcelnum}); end
        
        matchID = mean(grouptosubs_assignments(parcels==parcelIDs(parcelnum),s));
        
        subparcelindices = find(substogroup_assignments(:,s)==matchID) + strcmp('R',hem)*ncortexLverts;
        
        if isempty(subparcelindices)
            1;
        end
        
        subparceltimecourses(parcelnum,:) = mean(subtimecourse(subparcelindices,:),1);
        
        groupparcelindices = find(parcels==parcelIDs(parcelnum)) + strcmp('R',hem)*ncortexLverts;
        
        groupparceltimecourses(parcelnum,:) = mean(subtimecourse(groupparcelindices,:),1);
        
    end
    
    bigcorrelpatterns = paircorr_mod(subtimecourse',subparceltimecourses');
    
    sub_correlation_patterns(:,:,s) = bigcorrelpatterns([1:nnz(mask)] + strcmp('R',hem)*ncortexLverts,:);
    
    
    bigcorrelpatterns = paircorr_mod(subtimecourse',groupparceltimecourses');
    
    group_correlation_patterns(:,:,s) = bigcorrelpatterns([1:nnz(mask)] + strcmp('R',hem)*ncortexLverts,:);
    disp(' ')
end


output = zeros(length(mask),length(parcelIDs));
output(logical(mask),:) = mean(sub_correlation_patterns,3);

save(gifti(single(output)),'/data/cn4/evan/RestingState/Ind_variability/Subjects/Target250/Assigned_indparcel_avg_connectivity.func.gii')


output = zeros(length(mask),length(parcelIDs));
output(logical(mask),:) = mean(group_correlation_patterns,3);

save(gifti(single(output)),'/data/cn4/evan/RestingState/Ind_variability/Subjects/Target250/Assigned_groupparcel_avg_connectivity.func.gii')
        
        

