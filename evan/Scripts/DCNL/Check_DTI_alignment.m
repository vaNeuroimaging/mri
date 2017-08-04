subjects = {'101','102','113','118','120','122','125','127','132','138','147','150','151','154','156','159','160','161','162','172','181','182','187','202','207','211','214','215','221','225','229','232','233','242','250','254','255','272','274','112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','247','295','297','300','305','309','334','339','340','343','359','362','374','383','395','396'};

for sub = 1:length(subjects)
    subject = subjects{sub};
    
    DTIfolder = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subject '/DTI/'];
    copyfile([DTIfolder 'DTIStudioFA.hdr'],[DTIfolder 'TEMPFA.hdr']);
    copyfile([DTIfolder 'DTIStudioFA.img'],[DTIfolder 'TEMPFA.img']);
    copyfile([DTIfolder 'DTIStudioFirstB0.hdr'],[DTIfolder 'TEMPB0.hdr']);
    copyfile([DTIfolder 'DTIStudioFirstB0.img'],[DTIfolder 'TEMPB0.img']);
    
    load SPM8_Realign_only_template.mat
    
    matlabbatch{1}.spm.spatial.realign.estimate.data{1}{1} = [DTIfolder 'TEMPFA.img'];
    matlabbatch{1}.spm.spatial.realign.estimate.data{1}{2} = [DTIfolder 'TEMPB0.img'];
    
    batchfilename = [DTIfolder 'TEMPalignment.mat'];

    save(batchfilename, 'matlabbatch');

    spm_jobman('run',batchfilename)
    
    alignmentfile = dir([DTIfolder 'rp*.txt']);
    data = textread([DTIfolder alignmentfile.name]);
    
    disp(['Subject ' subject ': ' num2str(data(2,1)) ' ' num2str(data(2,2)) ' ' num2str(data(2,3)) ' ' num2str(data(2,4)) ' ' num2str(data(2,5)) ' ' num2str(data(2,6))])
    
    eval(['!rm ' DTIfolder 'TEMP*'])
    eval(['!rm ' DTIfolder 'rp*'])
end