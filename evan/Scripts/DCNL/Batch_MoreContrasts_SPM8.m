subjects = {'402','406','407','410','412','415','416','417','420'};
%{'101' '102'  '113' '118' '120' '122' '125' '127' '132' '138' '147' '150' '151' '154' '156' '159' '160' '161' '162' '166' '172' '181' '182' '187' '202' '207' '211' '214' '215' '221' '225' '229' '232' '233' '242' '250' '254' '255' '272' '274' '112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','281','247','295','297','300','305','309','334','339','340','343','359','362','383','395','396','374'};

for subject = 1:length(subjects)
    clear matlabbatch
    load /fmri/data3/Evan/Gene-Rest-Nback/Scripts/SPM8_Morecontrasts_Template.mat
    
    matlabbatch{1}.spm.stats.con.spmmat{1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subjects{subject} '/SPM8/Cond/SPM.mat'];
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.name = '2and3 vs Fixation';
    matlabbatch{1}.spm.stats.con.consess{1}.tcon.convec = [0 1 1 -2];
    
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.name = '2and3 vs 1';
    matlabbatch{1}.spm.stats.con.consess{2}.tcon.convec = [-2 1 1];
    
    batchfilename = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM8_MoreContrasts_Temp.mat'];
    save(batchfilename, 'matlabbatch');
    spm_jobman('run',batchfilename)
end