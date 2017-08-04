subjects = {'415'};
%{'101','102','113','118','120','122','125','127','132','138','147','150','151','154','156','159','160','161','162','166','172','181','182','187','202','207','211','214','215','221','225','229','232','233','242','250','254','255','272','274','112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','247','295','297','300','305','309','334','339','340','343','359','362','383','395','396','374','110','189','199','269'};
%'400','401','402','406','407','410','412','415','416','417','420'

for subject = 1:length(subjects)
    subjid = subjects{subject};
    disp(subjid)
    
    basedir = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subjid '/Struct/'];
    
%     structhdr = dir([basedir 'Siemens*.nii']);
%     %structimage = dir([basedir 'D*.img']);
%     
%     try rmdir([basedir 'FIRST_Hip_seg/'],'s');catch;end
%     mkdir([basedir 'FIRST_Hip_seg/']);
%     
%     copyfile([basedir structhdr(1).name],[basedir 'FIRST_Hip_seg/MRPAGE.nii']);
%     %copyfile([basedir structimage(1).name],[basedir 'FIRST_Hip_seg/MRPAGE.img']);
    
    cd([basedir 'FIRST_Hip_seg/']);
    
    eval(['!run_first_all -i ' basedir 'FIRST_Hip_seg/MPRAGE_reslice_step2.nii -s L_Hipp,R_Hipp -o Hip']);

    eval(['!fslsplit ' basedir 'FIRST_Hip_seg/Hip_all_fast_origsegs.nii.gz ' basedir 'FIRST_Hip_seg/Side -t']);

    copyfile([basedir 'FIRST_Hip_seg/Side0000.nii.gz'],[basedir 'FIRST_Hip_seg/' subjid '_L.nii.gz']);
    copyfile([basedir 'FIRST_Hip_seg/Side0001.nii.gz'],[basedir 'FIRST_Hip_seg/' subjid '_R.nii.gz']);
    
end
        
startup