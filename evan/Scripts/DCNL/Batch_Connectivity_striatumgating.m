subjects = {'101','102','113','118','120','122','125','127','132','138','147','150','151','154','156','159','160','161','162','166','172','181','182','187','202','207','211','214','215','221','225','229','232','233','242','250','254','255','272','274','112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327'};
    %'101','102','113','118','120','122','125','127','132','138','147','150','151','154','156','159','160','161','162','166','172','181','182','187','202','207','211','214','215','221','225','229','232','233','242','250','254','255','272','274'};
%

runs = {'Nback'};
conditions = {'1Back','2Back','3Back'};
condtimes = {[49 93 137] [5 115 159] [27 71 181]};
conddur = 11;

timepointstouse = [];
for cond = 1:length(conditions)
    for block = 1:3
        timepointstouse = [timepointstouse condtimes{cond}(block):(condtimes{cond}(block)+conddur-1)];
    end
end
timepointstouse = sort(timepointstouse);


load('taskregressor.mat');
taskregressor = taskregressor(timepointstouse,:);


% runs = {'Rest'};
% timepointstouse = [1:148];

header = 'sw';

gatingroi_name = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Activation/SPM8/TaskvRest_61/Striatum_bilat_p00001FWE_roi.mat';
seedroi_name = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Activation/SPM8/TaskvRest_61/Peak_SMA_roi.mat';

for subject = 1:length(subjects)

    clear timecourse residual_timecourse residuals

    clear datamatrix imgfiles fullimgfiles motionparams gatingtimes residuals;
    for run = 1:length(runs)

        subjsid = subjects{subject};


        data_dir = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/', subjsid, '/SPM8/' runs{run} '/'];

         motionparamfile = dir([data_dir 'rp*.txt']);
         motionparams = textread([data_dir motionparamfile.name]);
         motionparams = motionparams(timepointstouse,:);

        
        imgfiles = dir([data_dir header '*.hdr']);
        imgfiles = imgfiles(timepointstouse);
        
        m = size(imgfiles, 1);
        for j=1:m
            fullimgfiles(j, :) = [data_dir imgfiles(j).name];
        end
        
            
        
        clear Y
        

        try CSF_rois = maroi('load_cell',['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/Ind_ROIs/' subjsid '_LV_roi.mat']);
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
        
        
        seed_rois = maroi('load_cell',seedroi_name);
        [Y a b c] = getdata(seed_rois{1}, fullimgfiles,'l');
        seed_timecourse = mean(Y,2);
        
        clear Y
        
        gatingrois = maroi('load_cell', gatingroi_name);
            [Y a b c] = getdata(gatingrois{1}, fullimgfiles,'l');


        for voxelnum = 1:size(Y,2)
                [b, bint, r] = regress(Y(:,voxelnum),[CSF_timecourse WM_timecourse taskregressor ones(length(CSF_timecourse),1)]);
                %motionparams 
                residuals(:,voxelnum) = r;
               
        end
            
        gating_residual_timecourse = mean(residuals,2);
        
        sorted_gatingvals = sort(gating_residual_timecourse);
        
%         gatingtimes{1} = find(gating_residual_timecourse > sorted_gatingvals(66)); 
%         gatingtimes{2} = find(gating_residual_timecourse <= sorted_gatingvals(66) & gating_residual_timecourse > sorted_gatingvals(33)); 
%         gatingtimes{3} = find(gating_residual_timecourse <= sorted_gatingvals(33));
%         gatinglabels = {'Highthird','Medthird','Lowthird'};

        gatingtimes{1} = find(gating_residual_timecourse > sorted_gatingvals(49)); 
        gatingtimes{2} = find(gating_residual_timecourse <= sorted_gatingvals(49)); 
        gatinglabels = {'Highhalf','Lowhalf'};
        
        
        for division = 1:length(gatingtimes)
            clear matlabbatch
            load SPM8_Seedbased_Connectivity_template
            matlabbatch{1}.spm.stats.fmri_spec.dir{1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/', subjsid, '/SPM8/ActivatedStriatumgating_wTPN/' gatinglabels{division} '/'];
            
            for scan = 1:length(gatingtimes{division})
                matlabbatch{1}.spm.stats.fmri_spec.sess.scans{scan} = [fullimgfiles(gatingtimes{division}(scan),1:end-3) 'img,1'];
            end
            matlabbatch{1}.spm.stats.fmri_spec.sess.regress(1).val = seed_timecourse(gatingtimes{division});
            matlabbatch{1}.spm.stats.fmri_spec.sess.regress(2).val = WM_timecourse(gatingtimes{division});
            matlabbatch{1}.spm.stats.fmri_spec.sess.regress(3).val = CSF_timecourse(gatingtimes{division});
            matlabbatch{1}.spm.stats.fmri_spec.sess.regress(4).val = taskregressor(gatingtimes{division},1);
            matlabbatch{1}.spm.stats.fmri_spec.sess.regress(4).name = 'Task1';
            matlabbatch{1}.spm.stats.fmri_spec.sess.regress(5) = matlabbatch{1}.spm.stats.fmri_spec.sess.regress(4);
            matlabbatch{1}.spm.stats.fmri_spec.sess.regress(5).val = taskregressor(gatingtimes{division},2);
            matlabbatch{1}.spm.stats.fmri_spec.sess.regress(5).name = 'Task2';
            matlabbatch{1}.spm.stats.fmri_spec.sess.regress(6) = matlabbatch{1}.spm.stats.fmri_spec.sess.regress(4);
            matlabbatch{1}.spm.stats.fmri_spec.sess.regress(6).val = taskregressor(gatingtimes{division},3);
            matlabbatch{1}.spm.stats.fmri_spec.sess.regress(6).name = 'Task3';
            
%             for i=1:6
%                 matlabbatch{1}.spm.stats.fmri_spec.sess.regress(i+3).val = motionparams(gatingtimes{division},i);
%             end
            
            try;rmdir(matlabbatch{1}.spm.stats.fmri_spec.dir{1},'s');catch;end
            mkdir(matlabbatch{1}.spm.stats.fmri_spec.dir{1});
            
            batchfilename =  [matlabbatch{1}.spm.stats.fmri_spec.dir{1}(1:end-1) '.mat'];

            save(batchfilename, 'matlabbatch');

            spm_jobman('run',batchfilename)
                    
                    
        end
        
        
    end
end