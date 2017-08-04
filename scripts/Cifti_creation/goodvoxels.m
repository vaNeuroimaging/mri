function goodvoxels(funcvol, tmask, surfdir, outfile,subject)


neighsmooth = '5';
factor = .5;

ribbonname = [surfdir '/' subject '/7112b_fs_LR/Ribbon/ribbon_333.nii.gz'];
if ~exist(ribbonname)
    evalc(['!csh /data/cn/data1/scripts/CIFTI_RELATED/Cifti_creation/create_ribbon_singlesub.csh ' subject ' ' surfdir]);
end


data = load_nii([funcvol '.nii.gz']);
data_orig = data;
datamean = mean(data_orig.img(:,:,:,logical(tmask)),4);
data.img = datamean;
data.hdr.dime.dim(5) = 1;
meanname = [outfile subject '_mean'];
save_nii(data,[meanname '.nii']);
gzip([meanname '.nii'])
delete([meanname '.nii'])

datastd = std(data_orig.img(:,:,:,logical(tmask)),[],4);
data.img = datastd;
stdname = [outfile subject '_sd1'];
save_nii(data,[stdname '.nii']);
gzip([stdname '.nii'])
delete([stdname '.nii'])

clear data data_orig

covname = [outfile subject '_cov'];
evalc(['!fslmaths ' stdname ' -div ' meanname ' ' covname]);

evalc(['!fslmaths ' covname ' -mas ' ribbonname ' ' covname '_ribbon']);

Ribmean = num2str(str2num(evalc(['!fslstats ' covname '_ribbon -M'])));
evalc(['!fslmaths ' covname '_ribbon -div ' Ribmean ' ' covname '_ribbon_norm']);
evalc(['!fslmaths ' covname '_ribbon_norm -bin -s ' neighsmooth ' ' outfile subject '_SmoothNorm']);
evalc(['!fslmaths ' covname '_ribbon_norm -s ' neighsmooth ' -div ' outfile subject '_SmoothNorm -dilD ' covname '_ribbon_norm_s' neighsmooth]);
evalc(['!fslmaths ' covname ' -div ' Ribmean ' -div ' covname '_ribbon_norm_s' neighsmooth ' -uthr 1000 ' covname '_norm_modulate']);
evalc(['!fslmaths ' covname '_norm_modulate -mas ' ribbonname ' ' covname '_norm_modulate_ribbon']);

Final_Ribstd = str2num(evalc(['!fslstats ' covname '_norm_modulate_ribbon -S']));
Final_Ribmean = str2num(evalc(['!fslstats ' covname '_norm_modulate_ribbon -M']));

%Lower = Final_Ribmean - (Final_Ribstd * factor);
Upper = Final_Ribmean + (Final_Ribstd * factor);

evalc(['!fslmaths ' meanname ' -bin ' outfile subject '_mask']);
evalc(['!fslmaths ' covname '_norm_modulate -thr ' num2str(Upper) ' -bin -sub ' outfile subject '_mask -mul -1 ' outfile '/' subject '_goodvoxels']);

delete([meanname '.nii.gz'])
delete([stdname '.nii.gz'])
delete([covname '.nii.gz'])
delete([covname '_ribbon.nii.gz'])
delete([covname '_ribbon_norm.nii.gz'])
delete([outfile subject '_SmoothNorm.nii.gz'])
delete([covname '_ribbon_norm_s' neighsmooth '.nii.gz'])
delete([covname '_norm_modulate.nii.gz'])
delete([covname '_norm_modulate_ribbon.nii.gz'])
delete([outfile subject '_mask.nii.gz']);




