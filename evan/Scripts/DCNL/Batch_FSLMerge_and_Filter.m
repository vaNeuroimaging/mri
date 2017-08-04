subjects = {'110','189','199','269'};
    %'101','102','113','118','120','122','125','127','132','138','147','150','151','154','156','159','160','161','162','166','172','181','182','187','202','207','211','214','215','221','225','229','232','233','242','250','254','255','272','274','112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','247','295','297','300','305','309','334','339','340','343','359','362','383','395','396','374','400','401','402','406','407','410','412','415','416','417'};
%

%Names of runs to interrogate
runs = {'Nback'};

%Header of functional data to interrogate
header3d = 'sw';
header4d = 'smoothed.nii.gz';

for subject = 1:length(subjects)
    for run = 1:length(runs)
        disp([subjects{subject} ': ' runs{run}])
    
        data_dir = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/', subjects{subject}, '/SPM8/' runs{run} '/'];
    
        eval(['!fslmerge -t ' data_dir header4d ' ' data_dir header3d '*.img']);
    
        copyfile([data_dir header4d],[data_dir 'tempsmoothed.nii.gz']);
        
        eval(['!fslmaths ' [data_dir 'tempsmoothed.nii.gz'] ' -bptf 25 2.5 ' data_dir 'filt_' header4d]);
        eval(['!fslsplit ' data_dir 'filt_' header4d ' ' data_dir 'f' header3d ' -t']);
        gunzip({[data_dir 'f' header3d '*.nii.gz'],[data_dir header4d],[data_dir 'filt_' header4d]});
        delete([data_dir 'tempsmoothed.nii.gz']);
        
    end
end