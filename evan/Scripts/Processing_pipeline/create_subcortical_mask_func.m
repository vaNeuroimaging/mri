function create_subcortical_mask_func(freesurfdir,maskdir)

workbenchdir = '/data/cn4/laumannt/workbench/bin_linux64/';
evalc(['!cp ' freesurfdir '/mri/wmparc.mgz ' maskdir]);
cd(maskdir)
evalc(['!cp /data/cn4/laumannt/subcortical_mask/FreeSurferSubcorticalLabelTableLut* ' maskdir])
evalc(['!cp /data/cn4/laumannt/subcortical_mask/*.atlasroi.32k_fs_LR.shape.gii ' maskdir])
evalc(['!mri_convert -rl ' freesurfdir '/mri/rawavg.mgz wmparc.mgz wmparc.nii']);
evalc(['!niftigz_4dfp -4 wmparc.nii wmparc']);
evalc(['!t4img_4dfp none wmparc wmparc_333 -O333'])
evalc(['!niftigz_4dfp -n wmparc_333 wmparc_333'])
evalc(['!' workbenchdir '/wb_command -volume-label-import wmparc_333.nii.gz FreeSurferSubcorticalLabelTableLut_nobrainstem_LR.txt subcortical_mask_LR_333.nii -discard-others -unlabeled-value 0'])
evalc(['!' workbenchdir '/wb_command -volume-label-import wmparc_333.nii.gz FreeSurferSubcorticalLabelTableLut_nobrainstem_sub_L_cbll_R.txt subcortical_mask_sub_L_cbll_R.nii -discard-others -unlabeled-value 0'])
evalc(['!' workbenchdir '/wb_command -volume-label-import wmparc_333.nii.gz FreeSurferSubcorticalLabelTableLut_nobrainstem_sub_R_cbll_L.txt subcortical_mask_sub_R_cbll_L.nii -discard-others -unlabeled-value 0'])
evalc('!rm wmparc.mgz')