warning off

subjects = {'112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','247','295','297','300','305','309','334','339','340','343','359','362','383','395','396','374'};
%{'112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','247','295','297','300','305','309','334','339','340','343','359','362','383','395','396','374'};
    %

%,'400','401','402','406','407','410','412','415','416','417','420'};

roipath = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Seed-based/DAT_DC_FirstRest/2mm/';

rois = {'dmPFC','laIns','ldlPFC','raIns','rdlPFC'};

EPItemplate = '/data/apps/spm/8/templates/EPI.nii';

for subject = 1:length(subjects)
    subjectname = subjects{subject};
    DTIpath = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subjectname '/DTI/'];
    
    eval(['!bet ' DTIpath  'DTIStudioFirstB0.hdr ' DTIpath 'DTIStudioFirstB0_bet']);
    gunzip([DTIpath 'DTIStudioFirstB0_bet.nii.gz']);
    delete([DTIpath 'DTIStudioFirstB0_bet.nii.gz']);

    eval(['!flirt -in ' DTIpath 'DTIStudioFirstB0_bet -ref /data/apps/fsl/4.1.7/data/standard/MNI152_T1_2mm_brain -omat ' DTIpath 'diffspace2standard']);
    
    eval(['!convert_xfm -inverse ' DTIpath 'diffspace2standard -omat ' DTIpath 'standard2diffspace']);
    
%     for roinum=1:length(rois)
%          roiname = rois{roinum};
%         eval(['!flirt -in ' roipath roiname '.nii -ref ' DTIpath 'DTIStudioFirstB0_bet -applyxfm -init ' DTIpath 'standard2diffspace -out ' DTIpath roiname '_diffspace -interp nearestneighbour']);
%         gunzip([DTIpath roiname '_diffspace.nii.gz'])
%         delete([DTIpath roiname '_diffspace.nii.gz']);
%     end
    
    %delete([DTIpath 'DC_diffspace.nii']);
    %eval(['!flirt -in /fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/diMartino/DC_2mm.nii -ref ' DTIpath 'DTIStudioFirstB0_bet -applyxfm -init ' DTIpath 'standard2diffspace -out ' DTIpath 'DC_diffspace -interp nearestneighbour']);
    eval(['!flirt -in /fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/Caudate_body_bilat.nii -ref ' DTIpath 'DTIStudioFirstB0_bet -applyxfm -init ' DTIpath 'standard2diffspace -out ' DTIpath 'Caudate_diffspace -interp nearestneighbour']);
    %eval(['!fslmaths ' DTIpath 'DC_diffspace -thr .9 -bin ' DTIpath 'DC_diffspace']);
    gunzip([DTIpath 'Caudate_diffspace.nii.gz']);
    delete([DTIpath 'Caudate_diffspace.nii.gz']);
    

    
end
