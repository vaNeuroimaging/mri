

warning off

directory = pwd;

subjects = {'101','102','113','118','120','122','125','127','132','138','147','150','151','154','156','159','160','161','162','166','172','187','202','207','211','214','215','221','225','229','232','233','242','250','254','255','272','274'};

models = {'model5','model6','model7'};

load DCM_model_comparison.mat


for subject = 1:length(subjects)
    
    for modelnum = 1:length(models)
        model = models{modelnum};
        
        matlabbatch{1}.spm.stats.bms.bms_dcm.sess_dcm{subject}.mod_dcm{modelnum,1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subjects{subject} '/DCM_' model '.mat'];
        matlabbatch{1}.spm.stats.bms.bms_dcm.family_level.family(modelnum).family_name = model;
        matlabbatch{1}.spm.stats.bms.bms_dcm.family_level.family(modelnum).family_models = modelnum;
        
    end
    
    
    
end

save('/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/DCM/DCM_model_comparison.mat','matlabbatch');
spm_jobman('run','/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/DCM/DCM_model_comparison.mat')





