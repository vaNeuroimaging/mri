subcort_mask = '/home/data/scripts/Resources/cifti_masks/subcortical_mask_LR_333_MNI.nii.gz';
medial_mask_L = '/home/data/scripts/Resources/cifti_masks/L.atlasroi.32k_fs_LR.shape.gii';
medial_mask_R = '/home/data/scripts/Resources/cifti_masks/R.atlasroi.32k_fs_LR.shape.gii';
systemresult = cell(0,2);
fsLRfolder = ['/home/data/subjects/' subject '/fs_LR/MNI/'];
funcvol = 'MAV067_snr_notmask';

workbenchdir = '/usr/local/workbench/bin_rh_linux64/';
for hem = 1:2
midsurf = [fsLRfolder '/Native/' subject '.' HEMS{hem} '.midthickness.native.surf.gii'];
midsurf_LR32k = [fsLRfolder '/fsaverage_LR32k/' subject '.' HEMS{hem} '.midthickness.32k_fs_LR.surf.gii'];
whitesurf = [fsLRfolder '/Native/' subject '.' HEMS{hem} '.white.native.surf.gii'];
pialsurf = [fsLRfolder '/Native/' subject '.' HEMS{hem} '.pial.native.surf.gii'];
nativedefsphere = [fsLRfolder '/Native/' subject '.' HEMS{hem} '.sphere.reg.reg_LR.native.surf.gii'];
outsphere = [fsLRfolder '/fsaverage_LR32k/' subject '.' HEMS{hem} '.sphere.32k_fs_LR.surf.gii'];
surfname = [funcvol '_' HEMS{hem}];
%disp(['Subject ' subject ': mapping ' HEMS{hem} ' hemisphere data to surface']);
[systemresult{end+1,1},systemresult{end+1,2}] = system([workbenchdir '/wb_command -volume-to-surface-mapping ' funcvol '.nii.gz ' midsurf ' ' surfname '.func.gii -ribbon-constrained ' whitesurf ' ' pialsurf]);
%disp(['Subject ' subject ': Deforming ' HEMS{hem} ' hemisphere timecourse to 32k fs_LR']);
[systemresult{end+1,1},systemresult{end+1,2}] = system([workbenchdir '/wb_command -metric-resample ' surfname '.func.gii ' nativedefsphere ' ' outsphere ' ADAP_BARY_AREA ' surfname '_dil10_32k_fs_LR.func.gii -area-surfs ' midsurf ' ' midsurf_LR32k]);
surfname_final{hem} = [surfname '_dil10_32k_fs_LR.func.gii'];
[systemresult{end+1,1},systemresult{end+1,2}] = system(['/usr/local/caret/bin_linux64/caret_command -file-convert -format-convert XML_BASE64 ' surfname_final{hem}]);
delete([surfname '.func.gii']);
end
[systemresult{end+1,1},systemresult{end+1,2}] = system([workbenchdir '/wb_command -cifti-create-dense-timeseries ' funcvol '_LR_surf_subcort_333_32k_fsLR.dtseries.nii -volume ' funcvol '.nii.gz ' subcort_mask ' -left-metric ' surfname_final{1} ' -roi-left ' medial_mask_L ' -right-metric ' surfname_final{2} ' -roi-right ' medial_mask_R]);
delete(surfname_final{1});
delete(surfname_final{2});