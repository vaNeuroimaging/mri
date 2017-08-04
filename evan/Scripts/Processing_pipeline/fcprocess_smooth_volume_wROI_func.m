function fcprocess_smooth_volume_wROI_func(subject,funcdir,maskdir,smoothnum,outputdir)
% Function volumetrically smooths data within a mask, TOL 09/14

workbenchdir = '/data/cn4/evan/workbench/bin_linux64/';
    disp(['Volume smooothing, processing subject #' subject])
    funcvol = [outputdir '/' subject '_funcvol'];
    %[funcdir '/' subject '/' subject '_333_zmdt_resid_bpss_zmdt_g7'];
    %system(['niftigz_4dfp -n ' funcvol ' ' outputdir '/temp']);
    system([workbenchdir '/wb_command -volume-smoothing ' funcvol '.nii.gz ' num2str(smoothnum) ' ' funcvol '_wROI255.nii.gz -roi ' maskdir '/subcortical_mask_LR_333.nii'])
    %fslstr = [fslstr outputdir 'temp_sess' num2str(s) '_wROI255.nii.gz '];
    %sesdir = [funcdir '/' sessions{s}];
    %funcvol = [sesdir '/' sessions{s} '_333_zmdt_resid_ntrpl_bpss_zmdt'];
    %funcvol = [pth '/' name];
    %disp(['Volume smooothing, processing subject #' num2str(s)])
    %system([workbenchdir '/wb_command -volume-smoothing ' funcvol '.nii.gz ' num2str(smoothnum) ' ' funcvol '_wROI255.nii.gz -roi ' maskdir '/subcortical_mask_LR_333.nii'])

%end
%system(fslstr)
system('rm temp*.nii.gz')
%matlabpool close