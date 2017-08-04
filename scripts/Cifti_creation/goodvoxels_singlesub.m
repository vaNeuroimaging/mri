function systemresult = goodvoxels_singlesub(funcvol, surfdir, outfile,subject,sequence_name,systemresult)


neighsmooth = '5';
factor = .5;

ribbonname = [surfdir '/Ribbon/ribbon_333.nii.gz'];
%if ~exist(ribbonname)
    [systemresult{end+1,1},systemresult{end+1,2}] = system(['csh /home/data/scripts/Cifti_creation/create_ribbon_singlesub.csh ' subject ' ' surfdir ' T1_biascorr_MNI.nii.gz 3']);
%end

meanname = [outfile sequence_name '_mean'];
stdname = [outfile sequence_name '_sd1'];

[systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths ' funcvol ' -Tmean ' meanname]);
[systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths ' funcvol ' -Tstd ' stdname]);

covname = [outfile sequence_name '_cov'];
[systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths ' stdname ' -div ' meanname ' ' covname]);

[systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths ' covname ' -mas ' ribbonname ' ' covname '_ribbon']);

Ribmean = num2str(str2num(evalc(['!fslstats ' covname '_ribbon -M'])));
[systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths ' covname '_ribbon -div ' Ribmean ' ' covname '_ribbon_norm']);
[systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths ' covname '_ribbon_norm -bin -s ' neighsmooth ' ' outfile sequence_name '_SmoothNorm']);
[systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths ' covname '_ribbon_norm -s ' neighsmooth ' -div ' outfile sequence_name '_SmoothNorm -dilD ' covname '_ribbon_norm_s' neighsmooth]);
[systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths ' covname ' -div ' Ribmean ' -div ' covname '_ribbon_norm_s' neighsmooth ' -uthr 1000 ' covname '_norm_modulate']);
[systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths ' covname '_norm_modulate -mas ' ribbonname ' ' covname '_norm_modulate_ribbon']);

Final_Ribstd = str2num(evalc(['!fslstats ' covname '_norm_modulate_ribbon -S']));
Final_Ribmean = str2num(evalc(['!fslstats ' covname '_norm_modulate_ribbon -M']));

%Lower = Final_Ribmean - (Final_Ribstd * factor);
Upper = Final_Ribmean + (Final_Ribstd * factor);

[systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths ' meanname ' -bin ' outfile sequence_name '_mask']);
[systemresult{end+1,1},systemresult{end+1,2}] = system(['fslmaths ' covname '_norm_modulate -thr ' num2str(Upper) ' -bin -sub ' outfile sequence_name '_mask -mul -1 ' outfile '/' sequence_name '_goodvoxels']);

warning off
delete([meanname '.nii.gz'])
delete([stdname '.nii.gz'])
delete([covname '.nii.gz'])
delete([covname '_ribbon.nii.gz'])
delete([covname '_ribbon_norm.nii.gz'])
delete([outfile sequence_name '_SmoothNorm.nii.gz'])
delete([covname '_ribbon_norm_s' neighsmooth '.nii.gz'])
delete([covname '_norm_modulate.nii.gz'])
delete([covname '_norm_modulate_ribbon.nii.gz'])
delete([outfile sequence_name '_mask.nii.gz']);




