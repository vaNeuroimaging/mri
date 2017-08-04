
subjects = {'101','102','110','113','118','120','122','125','127','132','138','147','150','151','154','156','159','160','161','162','172','181','187','189','199','202','207','211','214','215','221','225','229','232','233','242','250','254','255','269','272','274','112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','247','295','297','300','305','309','334','339','340','343','359','362','383','395','396','374','400','401','402','406','407','410','412','415','416','417','420'};
%
%'181','182'
%166

subjects_ordered = {'101','102','110','112','113','118','120','122','125','126','127','132','133','137','138','147','150','151','154','156','159','160','161','162','172','181','187','189','199','202','207','208','211','214','215','221','222','225','227','229','232','233','242','247','250','251','253','254','255','256','258','260','261','264','269','270','272','274','279','281','283','292','295','297','300','301','305','307','309','322','327','334','339','340','343','359','362','374','383','395','396','400','401','402','406','407','410','412','415','416','417','420'}; 
%{'112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','247','295','297','300','305','309','334','339','340','343','359','362','383','395','396','374','400','401','402','406','407','410','412','415','416','417','420'};
%
%'181','182',

addpath /fmri/data3/Evan/Gene-Rest-Nback/Scripts

baseTBSSpath = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/TBSS/BDNF_92/';


try rmdir(baseTBSSpath,'s'); catch; end
mkdir(baseTBSSpath);
cd(baseTBSSpath)

needsaflip = [1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
%

for subject = 1:length(subjects)
    
    subjid = subjects{subject};
    
    if needsaflip(subject)
        
        copyfile(['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subjid '/DTI/bothruns/dti_FA.nii.gz'], [baseTBSSpath subjid '_dti_FA_wrong.nii.gz']);
        
        eval(['!fslswapdim ' baseTBSSpath subjid '_dti_FA_wrong.nii.gz -x -y z ' baseTBSSpath subjid '_dti_FA.nii.gz']);
        
        delete([baseTBSSpath subjid '_dti_FA_wrong.nii.gz']);
        
    else
        
        copyfile(['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subjid '/DTI/bothruns/dti_FA.nii.gz'], [baseTBSSpath subjid '_dti_FA.nii.gz']);
        
    
     
    end
    
end

disp('Step 1: Preprocessing')

eval('!tbss_1_preproc *.nii.gz');

disp('Step 2: Registration and normalization (this is the longest step)')
 
eval('!tbss_2_reg -T');

disp('Step 3: Skeletonize mean image')
 
eval('!tbss_3_postreg -S');

disp('Step 4: Skeletonize individual images')
 
eval('!tbss_4_prestats .2');


copyfile([baseTBSSpath 'stats/all_FA_skeletonised.nii.gz'],[baseTBSSpath 'stats/indsubs_FA_skeletonized.nii.gz']);

eval(['!cd ' baseTBSSpath 'stats/; fslsplit ' baseTBSSpath 'stats/indsubs_FA_skeletonized.nii.gz Subject -t']);

files = dir([baseTBSSpath 'stats/Subject*.nii.gz']);

for subject = 1:length(subjects)
    
    subjid = subjects_ordered{subject};
    
    eval(['!fslchfiletype NIFTI ' baseTBSSpath 'stats/' files(subject).name ' ' baseTBSSpath 'stats/' subjid '_skeletonized']);
    
end

delete([baseTBSSpath 'stats/Subject*.nii.gz']);
delete([baseTBSSpath 'stats/indsubs_FA_skeletonized.nii.gz']);

 
rmpath /fmri/data3/Evan/Gene-Rest-Nback/Scripts
 
 cd /fmri/data3/Evan/Gene-Rest-Nback/Scripts