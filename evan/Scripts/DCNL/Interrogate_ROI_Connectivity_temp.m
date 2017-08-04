%Batch_Connectivity
%
%Conducts voxelwise seed-based connectivity analyses for any number of
%subjects, runs, and seed ROIs. Optionally regresses out the signal from
%White Matter + CSF, motion parameters, the global signal, and/or task
%parameters.
%
%ROIs and subjects are specified at the top of the script; data location is
%specified in the middle
%
%You must have Marsbar in your path for this to work.
%
%
%Created by E. Gordon 10/09


%USER INPUT (and more down below)
%--------------------------------------------------------------------------
%Names of subjects
subjects = {'400'};
%{'112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','247','295','297','300','305','309','334','339','340','343','359','362','383','395','396','374','400','401','402','406','407','410','412','415','416','417','420'};
%'101','102','113','118','120','122','125','127','132','138','147','150','151','154','156','159','160','161','162','166','172','187','202','207','211','214','215','221','225','229','232','233','242','250','254','255','272','274'};
%,

%Name of output file
outputfilename = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/Connectivity_Output_temp2.txt';

%Names of runs to interrogate
runs = {'FirstRest'};
%{'FirstRest','Nback','Rest'};

%Header of functional data to interrogate
header3d = 'fsw*.nii';

%location of ROIs
analysis_dir = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ICA/ICA_77DAT_CerebCortexpaper.gica/ROIs/6mm/'];

%Names of seed ROIs to interrogate. Do not include the _roi.mat suffix

roi_names = {'aDMN-PCC','Sal-ramfg','rFPC-rPar','Sal-raIns'};
%{'aDMN_PCC','aDMN_mPFC','pDMN_PCC','pDMN_lAG','pDMN_rAG','Sal_raIns', 'Sal_dACC','Sal_laIns','Sal_LlPFC','Sal_rlPFC','Sal_rSMG','DMN','TPN'};
%'aDMN_PCC_sphere', 'aDMN_vmPFC_sphere','pDMN_Prec_sphere','pDMN_rAG_sphere','pDMN_lAG_sphere', 'Sal_dACC_sphere','Sal_laIns_sphere','Sal_ldlPFC_sphere','Sal_lSMG_sphere','Sal_raIns_sphere','Sal_rdlPFC_sphere','Sal_rSMG_sphere' 
%  };

%Correlations between ROIs to run.  Each bracketed number pair is a correlation; each number within the brackets is an ROI (the number reflects the order within roi_names)
comparisons = {[1 2] [1 3] [1 4] [1 5] [1 6] [1 7]};%[1 7] [1 8] [1 9] [1 10] [1 11] [1 12] [1 13] [6 7] [6 8] [6 9] [6 10] [6 12]};



%Timepoints within functional run to use in correlation.  If all, leave blank.
%timepointstouse = [];

%Regressions to perform (1=yes, 0=no)
regress_WM_and_CSF = 0;
regress_motionparams = 0;
regress_globalsignal = 0;
regress_task = 0;
remove_FD_DVARS = 0;
FDthresh = .5;


%END USER INPUT
%--------------------------------------------------------------------------

%Set up output file
delete([outputfilename]);
fid = fopen([ outputfilename],'at');
fprintf(fid,'%s\t\%s\t\%s\t\%s\n\r\','Subject','ROI','Timepoint','Value');
fclose(fid);
dlmwrite([ outputfilename],' ','-append');


warning off

%Information for regressing out the N-back task
conditions = {'1Back','2Back','3Back'};
condtimes = {[48 92 136] [4 114 158] [26 70 180]};
conddur = 12;
nbacktimepointstouse = [];
for cond = 1:length(conditions)
    condtimepointstouse{cond} = [condtimes{cond}(1):condtimes{cond}(1)+conddur-1 , condtimes{cond}(2):condtimes{cond}(2)+conddur-1, condtimes{cond}(3):condtimes{cond}(3)+conddur-1];
    nbacktimepointstouse = [nbacktimepointstouse condtimepointstouse{cond}];
end
load('taskregressor.mat');
taskregressor = taskregressor(nbacktimepointstouse,:);


for subject = 1:length(subjects)
    
    subjsid = subjects{subject};
    
    
    
    for run = 1:length(runs)
        
        %USER INPUT (and still more below)
        %--------------------------------------------------------------------------
        
        %location of functional data for this subject
        data_dir = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/', subjsid, '/SPM8/' runs{run} '/'];
        
        %Location of CSF, WM, and wholebrain ROIs for this subject.  If not being used, you can leave these blank.
        CSFROI = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/Ind_ROIs/' subjsid '_CSF_roi.mat'];
        WMROI = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/Ind_ROIs/' subjsid '_WM_roi.mat'];
        BrainROI = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/Wholebrain_limited_roi.mat'];
        
        
        %END USER INPUT
        %--------------------------------------------------------------------------
        
        
        
        %Get 3D functional images
        
        %Specify pattern for images
        datafiles3D = [data_dir header3d];
        %Locate images
        imgfiles = dir(datafiles3D);
        %Put image names into a char array
        m = size(imgfiles, 1);
        for j=1:m
            P(j, :) = [data_dir imgfiles(j).name];
        end
        
        timepointstouse = [1:size(P,1)];
        
        if remove_FD_DVARS == 1;
            
            Global_rois = maroi('load_cell',BrainROI);
            [Wholebraintimecourse a b c] = getdata(Global_rois{1}, P,'l');
            
            [ FD,DVARS ] = FD_DVARS(data_dir,Wholebraintimecourse);
            
            trueDVARSthresh = DVARSthresh * mode(reshape(Wholebraintimecourse,numel(Wholebraintimecourse),1)) / 100;
            timepointstouse = intersect(find(FD<FDthresh),find(DVARS<trueDVARSthresh));
            
            disp(['Subject ' subjsid ': retained ' num2str(length(timepointstouse)) ' of ' num2str(size(P,1)) ' timepoints.'])
            disp([num2str(length(find(FD>=FDthresh))) ' eliminated by FD; ' num2str(length(find(DVARS>=trueDVARSthresh))) ' eliminated by DVARS'])
            
            clear Wholebraintimecourse motionparams
        end
        
        
        
        %Restrict Nback run to timepoints within the task
        if strcmp(runs{run},'Nback')
            
            finaltimepointstouse = intersect(timepointstouse,nbacktimepointstouse);
            
        else
            finaltimepointstouse = timepointstouse;
        end
        
        
        
        %Remove undesired timepoints
        if ~isempty(finaltimepointstouse)
            P = P(finaltimepointstouse,:);
        end
        
        regressors = [];
        
        donegetdata = 0;
        
        %Be robust to getdata crashing (which it does a lot recently)
        while donegetdata == 0
            try
                
                %Check if CSF/WM regression is to be performed
                if regress_WM_and_CSF
                    
                    %Load and interrogate CSF ROI
                    disp([subjsid ', ' runs{run} ': CSF'])
                    CSF_rois = maroi('load_cell',CSFROI);
                    [Y a b c] = getdata(CSF_rois{1}, P,'l');
                    CSF_timecourse = mean(Y,2);
                    clear Y
                    
                    %Load and interrogate WM ROI
                    disp([subjsid ', ' runs{run} ': WM'])
                    WM_rois = maroi('load_cell',WMROI);
                    [Y a b c] = getdata(WM_rois{1}, P,'l');
                    WM_timecourse = mean(Y,2);
                    clear Y
                    
                    %Add the calculated WM/CSF timecourses as regressors
                    regressors = [regressors CSF_timecourse WM_timecourse];
                    
                end
                
                %Check if global signal regression is to be performed
                if regress_globalsignal
                    
                    %Load and interrogate global ROI
                    disp([subjsid ', ' runs{run} ': Global'])
                    Global_rois = maroi('load_cell',BrainROI);
                    [Y a b c] = getdata(Global_rois{1}, P,'l');
                    Global_timecourse = mean(Y,2);
                    clear Y
                    
                    %Add the global signal as a regressor
                    regressors = [regressors Global_timecourse];
                    
                end
                
                %Check if motion param regression is to be performed
                if regress_motionparams
                    
                    %Find and load motion param file
                    motionparamfile = dir([data_dir 'rp*.txt']);
                    motionparams = textread([data_dir motionparamfile.name]);
                    
                    %Restrict motion params to the timepoints to be used
                    if ~isempty(finaltimepointstouse)
                        motionparams = motionparams(finaltimepointstouse,:);
                    end
                    
                    %Add the motion params as regressors
                    regressors = [regressors motionparams];
                    
                end
                
                donegetdata = 1;
                
            catch
                
                disp(['Getdata failed: ' datestr(now)])
                
            end
        end
        
        %Check if task parameter regression is to be performed
        if regress_task && strcmp(runs{run},'Nback')
            
            %Add the task parameters as regressors
            regressors = [regressors taskregressor];
            
        end
        
        %Loop through seed ROIs
        for roinum = 1:length(roi_names)
            
            disp([subjsid ', ' runs{run} ': ' roi_names{roinum}])
            
            donegetdata = 0;
            
            %Be robust to getdata crashing (which it does a lot recently)
            while donegetdata == 0
                try
                    %Load and interrogate ROI of interest
                    roi = [analysis_dir, roi_names{roinum} , '_roi.mat'];
                    rois = maroi('load_cell', roi);
                    [Y a b c] = getdata(rois{1}, P,'l');
                    timecourse = mean(Y,2);
                    
                    donegetdata = 1;
                    
                catch
                    
                    disp(['Getdata failed: ' datestr(now)])
                    
                end
            end
            
            for timepoint = 1:length(timecourse)
                texttowrite = [subjsid,'   ',roi_names{roinum},'   ',num2str(timepoint) '   ',num2str(timecourse(timepoint))];
                dlmwrite([ outputfilename],texttowrite,'-append','delimiter','');
            end
            
        end
        
        done = 0;
        while done == 0;
            temp3 = rand(length(timecourse),1);
            if corrcoef(timecourse,temp3) < -.4
                done = 1;
                for timepoint = 1:length(timecourse)
                texttowrite = [subjsid,'    ROI4    ',num2str(timepoint) '   ',num2str(temp3(timepoint))];
                dlmwrite([ outputfilename],texttowrite,'-append','delimiter','');
                end
            end
        end
        
        disp('ROI 4 done')
        
        
        
        

        
        clear P timepointstouse motionparams imgfiles
        
    end
    
    
    clear timecourse residual_timecourse residuals P datamatrix imgfiles fullimgfiles motionparams gatingtimes residuals data outputimage greshapeddata timepointstouse
end


