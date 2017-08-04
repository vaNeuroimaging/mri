
tmasklist = '/data/hcp-bluearc/home/laumannt/120_parcellation/120_NEW_TMASKLIST.txt';
iscifti = 1;

[subjects tmasks] = textread(tmasklist,'%s%s');

hems = {'L','R'};

for s = 1:length(subjects)
    disp(subjects{s})
    subfilename = ['/data/hcp-bluearc/home/laumannt/120_parcellation/cifti_timeseries_normalwall/' subjects{s} '_BOLD_LR_surf_subcort_32k_fsLR_smooth2.55.func.gii'];
    subtimecourse = gifti(subfilename); subtimecourse = subtimecourse.cdata;
    tmask = load(tmasks{s});
    subtimecourse = subtimecourse(:,logical(tmask));
    subtimecourse(isnan(subtimecourse)) = 0;
    
    for hemnum = 1:length(hems)
        hem = hems{hemnum};
        
        parcelfilename = ['/data/cn4/evan/RestingState/FC_Mapping_120/120_' hem '_wateredgethresh_watershedmerge_0.45.func.gii'];
        parcels = gifti(parcelfilename); parcels = parcels.cdata;
        origparcels = parcels;
        
        maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii']; mask = gifti(maskname); mask = ~mask.cdata;
        parcels = parcels(logical(mask));
        maskL = gifti('/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii');ncortexLverts = nnz(maskL.cdata==0);
        
        parcelIDs = unique(parcels);
        parcelIDs(parcelIDs==0)=[];
        
        if s==1
            hem_parcelcorrelpatterns{hemnum} = zeros(length(parcelIDs),size(subtimecourse,1));
        end
        
        for parcelnum = 1:length(parcelIDs)
            ind = find(origparcels==parcelIDs(parcelnum));
            parcelindices{parcelnum} = find(parcels==parcelIDs(parcelnum)) + (strcmp(hem,'R') * ncortexLverts);
            parceltimecourse = mean(subtimecourse(parcelindices{parcelnum},:),1);
            parcelcorrelpattern = paircorr_mod(parceltimecourse',subtimecourse');
            parcelcorrelpattern(isnan(parcelcorrelpattern)) = 0;
            hem_parcelcorrelpatterns{hemnum}(parcelnum,:) = hem_parcelcorrelpatterns{hemnum}(parcelnum,:) + (parcelcorrelpattern/length(subjects));
        end
    end
end

for hemnum = 1:length(hems)
        hem = hems{hemnum};
        parcelcorrelpatterns = hem_parcelcorrelpatterns{hemnum};
        parcelcorrelpatternfile = ['/data/cn4/evan/RestingState/FC_Mapping_120/Group_parcelcorrelpatterns_' hem '.mat'];
        save(parcelcorrelpatternfile,'parcelcorrelpatterns','-v7.3')
end
            
            
