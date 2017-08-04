%function Match_group_parcels_to_ind_data(parcelfilename,allsubs_correlmatfile,datalist,tmasklist,hem,iscifti)
% Match_group_parcels_to_ind_data(parcelfilename,allsubs_correlmat,datalist,tmasklist,hem,iscifti)

parcelfilename = '120_L_crossthresh_watershedmerge.func.gii';
datalist = 'AllC_DATALIST.txt';
tmasklist = 'AllC_TMASKLIST.txt';
hem = 'L';
iscifti = 2;


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

parcelIDs = unique(parcels); parcelIDs(parcelIDs==0) = [];
for parcelnum = 1:length(parcelIDs)
    parcelindices{parcelnum} = find(parcels==parcelIDs(parcelnum)) + (strcmp(hem,'R') * ncortexLverts);
end





[subjects subdata] = textread(datalist,'%s%s');
if ~isempty(tmasklist)
    [subjects tmasks] = textread(tmasklist,'%s%s');
end



for subnum = 1:length(subjects)
    
    disp(['Subject ' num2str(subnum)])
    
    evalc(['!wb_command -cifti-convert -to-gifti-ext ' subdata{subnum} ' Temp.func.gii']);
    subtimecourse = gifti('Temp.func.gii'); subtimecourse = subtimecourse.cdata;
    if ~isempty(tmasklist)
        tmask = load(tmasks{subnum});
        subtimecourse = subtimecourse(:,logical(tmask));
    end
    subtimecourse(isnan(subtimecourse)) = 0;
    
    if subnum==1
        parcelcorrelpatterns = zeros(size(subtimecourse,1),length(subjects),length(parcelIDs));
    end
    
    for parcelnum = 1:length(parcelIDs)
        parceltimecourse = mean(subtimecourse(parcelindices{parcelnum},:),1);
        
        parcelcorrelpatterns(:,subnum,parcelnum) = paircorr_mod(subtimecourse',parceltimecourse');
    end
end

parcelcorrelpatterns(isnan(parcelcorrelpatterns)) = 0;

outputdata_firsteigval = zeros(size(mask));
outputdata30 = zeros(size(mask));
outputdata50 = zeros(size(mask));
outputdata70 = zeros(size(mask));
outputdata90 = zeros(size(mask));

for parcelnum = 1:length(parcelIDs)

    [ign ign2 eigvals_per] = PCA_reduction(squeeze(parcelcorrelpatterns(:,:,parcelnum)),'comps',2);
    totalexplained = cumsum(eigvals_per);
    
    outputdata_firsteigval(origparcels==parcelIDs(parcelnum)) = eigvals_per(1);
        
    numexplaining30 = min(find(totalexplained>30));
    outputdata30(origparcels==parcelIDs(parcelnum)) = numexplaining30;
    
    numexplaining50 = min(find(totalexplained>50));
    outputdata50(origparcels==parcelIDs(parcelnum)) = numexplaining50;
    
    numexplaining70 = min(find(totalexplained>70));
    outputdata70(origparcels==parcelIDs(parcelnum)) = numexplaining70;
    
    numexplaining90 = min(find(totalexplained>90));
    outputdata90(origparcels==parcelIDs(parcelnum)) = numexplaining90;
    
end

save(gifti(single(outputdata_firsteigval)),['First_eigval_sub_variance_' hem '.func.gii']);
save(gifti(single(outputdata30)),['NumComps_explaining_30pct_sub_variance_' hem '.func.gii']);
%save(gifti(single(outputdata50)),['NumComps_explaining_50pct_sub_variance_' hem '.func.gii']);
%save(gifti(single(outputdata70)),['NumComps_explaining_70pct_sub_variance_' hem '.func.gii']);
%save(gifti(single(outputdata90)),['NumComps_explaining_90pct_sub_variance_' hem '.func.gii']);


    
    
    

