SPMfolder = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Nonstationarity/dACC_Rest_vs_Nback/';

ROIroot = 'Nback_Rest';

Imagetomask = 'spmT_0002.hdr';

Masks = dir([SPMfolder ROIroot '*_roi.mat']);

FSLstring = ['fslmaths ' SPMfolder Imagetomask ' -mul 0 '];

for mask = 1:length(Masks)
    
    roi = maroi('load', [SPMfolder Masks(mask).name]);
    sp = mars_space([SPMfolder Imagetomask]);
    save_as_image(roi, [SPMfolder Masks(mask).name(1:end-8) '.nii'], sp);
    
    FSLstring = [FSLstring '-add ' SPMfolder Masks(mask).name(1:end-8) '.nii '];
    
end

eval(['!' FSLstring '-bin -mul ' SPMfolder Imagetomask ' ' SPMfolder Imagetomask(1:end-4) '_' ROIroot '_clusterthresh; fslchfiletype ANALYZE ' SPMfolder Imagetomask(1:end-4) '_' ROIroot '_clusterthresh']);
delete([SPMfolder ROIroot '*.nii'])



