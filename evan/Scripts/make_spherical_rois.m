function make_spherical_rois(roilist)

maskfile='/data/cn4/evan/ROIs/glm_atlas_mask_333.nii';

[ign roinum] = system(['cat ' roilist ' | wc -l']);
roinum = str2num(roinum);

for roi = 1:roinum
[ign roiname] = system(['awk -F '' '' ''NR==' num2str(s) '{printf $1}'' ' roilist]);
[ign xdim] = system(['awk -F '' '' ''NR==' num2str(s) '{printf $2}'' ' roilist]);
[ign ydim] = system(['awk -F '' '' ''NR==' num2str(s) '{printf $3}'' ' roilist]);
[ign zdim] = system(['awk -F '' '' ''NR==' num2str(s) '{printf $4}'' ' roilist]);


system(['fslmaths ' maskfile ' -roi ' xdim ' 1 ' ydim ' 1 ' zdim ' 1 0 1 -odt float -kernel sphere 10 -fmean -bin -mul ' maskfile ' ' pwd filesep roiname]);
system(['fslchfiletype NIFTI ' roiname 'nii.gz']);

end