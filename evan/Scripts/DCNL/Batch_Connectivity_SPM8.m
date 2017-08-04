subjects = {'112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','247','295','297','300','305','309','334','339','340','343','359','362','383','395','396','374','400','401','402','406','407','410','412','415','416','417','420'};
%subjects = {'101','102','113','118','120','122','125','127','132','138','147','150','151','154','156','159','160','161','162','166','172','181','182','187','202','207','211','214','215','221','225','229','232','233','242','250','254','255','272','274','112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','247','295','297','300','305','309','334','339','340','343','359','362','383','395','396','374'};
%
roi_location = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/diMartino/';
roi_names = {'DC_bilat'};
%,'Putamen_bilat'

runs = {'FirstRest'};
%'Rest',
% conditions = {'1Back','2Back','3Back'};
% condtimes = {[49 93 137] [5 115 159] [27 71 181]};
% conddur = 11;
% 
% timepointstouse = [];
% for cond = 1:length(conditions)
%     for block = 1:3
%         timepointstouse = [timepointstouse condtimes{cond}(block):(condtimes{cond}(block)+conddur-1)];
%     end
% end
% timepointstouse = sort(timepointstouse);
% 
% 
% load('taskregressor.mat');
% taskregressor = taskregressor(timepointstouse,:);


header = 'fsw';


for subject = 1:length(subjects)
    
    
    for run = 1:length(runs)
        
        clear timecourse residual_timecourse residuals
        
        clear datamatrix imgfiles fullimgfiles motionparams gatingtimes residuals
        
        subjsid = subjects{subject};
        disp(subjsid)
        
        data_dir = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/', subjsid, '/SPM8/' runs{run} '/'];
        
        motionparamfile = dir([data_dir 'rp*.txt']);
        motionparams = textread([data_dir motionparamfile.name]);
        
        
        imgfiles = dir([data_dir header '*.nii']);
        
        if strcmp(runs{run},'Nback')
            motionparams = motionparams(timepointstouse,:);
            imgfiles = imgfiles(timepointstouse);
        end
        
        m = size(imgfiles, 1);
        for j=1:m
            fullimgfiles(j, :) = [data_dir imgfiles(j).name];
        end
        
        clear Y
        
        try
            CSF_rois = maroi('load_cell',['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/Ind_ROIs/' subjsid '_LV_roi.mat']);
        catch
            CSF_rois = maroi('load_cell',['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/Ind_ROIs/' subjsid '_CSF_roi.mat']);
        end
        [Y a b c] = getdata(CSF_rois{1}, fullimgfiles,'l');
        CSF_timecourse = mean(Y,2);
        
        clear Y
        
        
        WM_rois = maroi('load_cell',['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/Ind_ROIs/' subjsid '_WM_roi.mat']);
        [Y a b c] = getdata(WM_rois{1}, fullimgfiles,'l');
        WM_timecourse = mean(Y,2);
        
        clear Y
        
%         WB_rois = maroi('load_cell',['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/ICA79_wholebrainmask_roi.mat']);
%         [Y a b c] = getdata(WB_rois{1}, fullimgfiles,'l');
%         WB_timecourse = mean(Y,2);
        
        for roi_num = 1:length(roi_names)
            
            seedroi_name = [roi_location roi_names{roi_num} '_roi.mat'];
            
            seed_rois = maroi('load_cell',seedroi_name);
            [Y a b c] = getdata(seed_rois{1}, fullimgfiles,'l');
            seed_timecourse = mean(Y,2);
            
            clear Y
            
            clear matlabbatch
            %load SPM8_Seedbased_Connectivity_template
            load SPM8_Analysis_FirstRest
            %matlabbatch{1}.spm.stats.fmri_spec.dir{1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/', subjsid, '/SPM8/' roi_names{roi_num} '_' runs{run} '/'];
            matlabbatch{1}.spm.stats.fmri_spec.dir{1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/', subjsid, '/SPM8/FirstRest_forDCM/'];
            for scan = 1:length(seed_timecourse)
                matlabbatch{1}.spm.stats.fmri_spec.sess.scans{scan} = [fullimgfiles(scan,1:end-3) 'nii,1'];
            end
            matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1).val = seed_timecourse;
%             matlabbatch{1}.spm.stats.fmri_spec.sess.regress(2).val = WB_timecourse;
%             matlabbatch{1}.spm.stats.fmri_spec.sess.regress = matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1:2);
%             matlabbatch{1}.spm.stats.fmri_spec.sess.regress(2).val = WM_timecourse;
%             matlabbatch{1}.spm.stats.fmri_spec.sess.regress(3).val = CSF_timecourse;
%             
%             
%             if strcmp(runs{run},'Nback')
%                 matlabbatch{1}.spm.stats.fmri_spec.sess.regress(4).val = taskregressor(:,1);
%                 matlabbatch{1}.spm.stats.fmri_spec.sess.regress(4).name = 'Task1';
%                 matlabbatch{1}.spm.stats.fmri_spec.sess.regress(5) = matlabbatch{1}.spm.stats.fmri_spec.sess.regress(4);
%                 matlabbatch{1}.spm.stats.fmri_spec.sess.regress(5).val = taskregressor(:,2);
%                 matlabbatch{1}.spm.stats.fmri_spec.sess.regress(5).name = 'Task2';
%                 matlabbatch{1}.spm.stats.fmri_spec.sess.regress(6) = matlabbatch{1}.spm.stats.fmri_spec.sess.regress(4);
%                 matlabbatch{1}.spm.stats.fmri_spec.sess.regress(6).val = taskregressor(:,3);
%                 matlabbatch{1}.spm.stats.fmri_spec.sess.regress(6).name = 'Task3';
%             end
%             
%             regressorssofar = length(matlabbatch{1}.spm.stats.fmri_spec.sess.regress);
%             
%             for i=1:6
%                 matlabbatch{1}.spm.stats.fmri_spec.sess.regress(i+regressorssofar) = matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1);
%                 matlabbatch{1}.spm.stats.fmri_spec.sess.regress(i+regressorssofar).name = ['Motion' num2str(i)];
%                 matlabbatch{1}.spm.stats.fmri_spec.sess.regress(i+regressorssofar).val = motionparams(:,i);
%             end
            
%             ran = 0;
%             while ran ==0
                try
                    try;rmdir(matlabbatch{1}.spm.stats.fmri_spec.dir{1},'s');catch;end
                    mkdir(matlabbatch{1}.spm.stats.fmri_spec.dir{1});
                    
                    batchfilename = [matlabbatch{1}.spm.stats.fmri_spec.dir{1}(1:end-1) '.mat'];
                    
                    save(batchfilename, 'matlabbatch');
                    
                    spm_jobman('run',batchfilename)
                    
                    ran = 1;
                    
                catch
                    error(['Failure; retry  ' datestr(now)])
                end
            %end
        end
        
        
        
        
    end
end