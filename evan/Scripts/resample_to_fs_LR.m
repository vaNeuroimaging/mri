%% resample areal distortion map from native to 32k fs_LR
dir = '/data/hcp-zfs/shared-nil/laumannt/120_parcellation/';
outputdir = '/data/cn4/laumannt/Poldrome/Poldrome_registration';
%[subjects tmasks] = textread([dir '/120_NEW_TMASKLIST.txt'],'%s%s');
subjects = {'sub013'};
%surfdir = '/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/';
surfdir = '/data/hcp-zfs/shared-nil/laumannt/Poldrome/shared_for_washu/freesurfer_washu/FREESURFER_fs_LR/';
HEMS = {'L';'R'};

for s = 1:length(subjects)
    
    subject = subjects{s};
    surfdirsub = [surfdir '/' subject];
    
    
    for hem = 1:2
        midsurf = [surfdirsub '/7112b_fs_LR/Native/' subject '.' HEMS{hem} '.midthickness.native.surf.gii'];
        midsurf_LR32k = [surfdirsub '/7112b_fs_LR/fsaverage_LR32k/' subject '.' HEMS{hem} '.midthickness.32k_fs_LR.surf.gii'];
        whitesurf = [surfdirsub '/7112b_fs_LR/Native/' subject '.' HEMS{hem} '.white.native.surf.gii'];
        pialsurf = [surfdirsub '/7112b_fs_LR/Native/' subject '.' HEMS{hem} '.pial.native.surf.gii'];
        nativedefsphere = [surfdirsub '/7112b_fs_LR/Native/' subject '.' HEMS{hem} '.sphere.reg.reg_LR.native.surf.gii'];
        outsphere = [surfdirsub '/7112b_fs_LR/fsaverage_LR32k/' subject '.' HEMS{hem} '.sphere.32k_fs_LR.surf.gii'];
        
        surfname = [subject '.' HEMS{hem} '.ArealDistortion'];
        system([ workbenchdir '/wb_command -metric-resample ' surfdirsub '/7112b_fs_LR/Native/' surfname '.native.shape.gii ' nativedefsphere ' ' outsphere ' ADAP_BARY_AREA ' outputdir '/' surfname '_32k_fs_LR.shape.gii -area-surfs ' midsurf ' ' midsurf_LR32k]);
        evalc(['!cp ' outputdir '/' surfname '_32k_fs_LR.shape.gii ' surfdirsub '/fsaverage_LR32k']);
   end
    surfname_L = [subject '.L.ArealDistortion_32k_fs_LR'];
    surfname_R = [subject '.R.ArealDistortion_32k_fs_LR'];
    
    system([workbenchdir '/wb_command -cifti-create-dense-timeseries ' outputdir '/' subject '.LR.ArealDistortion.32k_fs_LR.dtseries.nii -left-metric ' outputdir '/' surfname_L '.shape.gii -right-metric ' outputdir '/' surfname_R '.shape.gii'])
    
    
end

%% Bring all distortion maps into 1 file
outputdir = '/data/cn4/laumannt/Poldrome/Poldrome_registration';
[subjects tmasks] = textread([dir '/120_NEW_TMASKLIST.txt'],'%s%s');

for s = 1:length(subjects)
    subject = subjects{s};
    ciftistruct = ft_read_cifti_mod([outputdir '/' subject '.LR.ArealDistortion.32k_fs_LR.dtseries.nii']);
    
    areal(:,s) = ciftistruct.data;
    
    
end

ciftistruct.data = areal;
ft_write_cifti_mod([outputdir '/120subs_ArealDistortion.32k_fs_LR'],ciftistruct)

%% Compute zscore of areal distortion for poldrack data

ciftistruct = ft_read_cifti_mod([outputdir '/sub013.LR.ArealDistortion.32k_fs_LR.dtseries.nii']);
poldrome_areal = ciftistruct.data;
areal_std = std(areal,[],2);
areal_mean = mean(areal,2);

poldrome_areal_zscore = (poldrome_areal-areal_mean)./areal_std;
   
ciftistruct.data = poldrome_areal_zscore;
ft_write_cifti_mod([outputdir '/poldrome_ArealDistortion_zscore.32k_fs_LR'],ciftistruct);

%% Convert z to p-value

onetailed = 1-normcdf(abs(poldrome_areal_zscore));
twotailed = 2*onetailed;

p_fdr = FDR(twotailed,.05);
 ciftistruct.data = twotailed<p_fdr;%.*size(poldrome_areal,1);
 %ciftistruct.data = ciftistruct.data';
 ft_write_cifti_mod([outputdir '/poldrome_ArealDistortion_pval_FDRcorr.32k_fs_LR'],ciftistruct);

