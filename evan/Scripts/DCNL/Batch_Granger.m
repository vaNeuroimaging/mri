subjects = {'112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','247','295','297','300','305','309','334','339','340','343','359','362','383','395','396','374','400','401','402','406','407','410','412','415','416','417','420'};
    

runs = {'FirstRest'};%,'Rest','Nback'};

%Header of functional data to interrogate
header3d = 'fsw*.nii';
header4d = 'filt_smoothed.nii.gz';

%location of seed ROIs
roi_dir = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/Granger/'];


%Names of seed ROIs to interrogate. Do not include the _roi.mat suffix
roi_names = {'raIns','dACC','laIns','lPFC','aDMN_mPFC','pDMN_lAG','pDMN_PCC','pDMN_rAG'};

regress_WM_and_CSF = 1;
regress_motionparams = 1;
regress_globalsignal = 0;
regress_task = 1;
remove_FD = 0;
    FDthresh = .5;
    
outputfilename = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/Granger/GrangerOutput.txt';
    
delete([outputfilename]);
fid = fopen([outputfilename],'at');
fprintf(fid,'%s\t\%s\t\%s\n\r\','Subject','ROI_connection','Parameter');
fclose(fid);
dlmwrite([outputfilename],' ','-append');


for subject = 1:length(subjects)
    
    subjsid = subjects{subject};
        
    disp(subjsid)
    
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
        %------------------------------------------------------------------
        %--------
        
        
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
        
        if remove_FD == 1;
           
           FD = Calc_FD(data_dir);
           
           timepointstouse = find(FD<FDthresh);
           
           disp(['Subject ' subjsid ', ' runs{run} ': retained ' num2str(length(timepointstouse)) ' of ' num2str(size(P,1)) ' timepoints.'])
           clear motionparams 
        end
        
        finaltimepointstouse = timepointstouse;
        
         %Remove undesired timepoints
        if ~isempty(finaltimepointstouse)
            P = P(finaltimepointstouse,:);
        end
        
        
        %Get 4D functional data file
        
        %Specify data file name
        datafile4D =  [data_dir header4d];
        %Unzip it if it's a zipped nifti
        try gunzip(datafile4D); catch; end;
        %Load the unzipped data
        data = load_nii(datafile4D(1:end-3));
        %Save the data as a matrix
        datamatrix = data.img;
        %Remove undesired timepoints
        if ~isempty(finaltimepointstouse)
            datamatrix = datamatrix(:,:,:,finaltimepointstouse);
        end
        
        
        %Reshape the 4D data to be 2D (time x voxels)
        reshapeddata = double(reshape(datamatrix,size(datamatrix,1)*size(datamatrix,2)*size(datamatrix,3),size(datamatrix,4))');
        
        
        regressors = [];
        
        donegetdata = 0;
        
        %Be robust to getdata crashing (which it does a lot recently)
        while donegetdata == 0
            try
                
                %Check if CSF/WM regression is to be performed
                if regress_WM_and_CSF
                    
                    %Load and interrogate CSF ROI
                    CSF_rois = maroi('load_cell',CSFROI);
                    [Y a b c] = getdata(CSF_rois{1}, P,'l');
                    CSF_timecourse = mean(Y,2);
                    clear Y
                    
                    %Load and interrogate WM ROI
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
        
        seed_data = zeros(length(finaltimepointstouse), length(roi_names));
        
        for roinum = 1:length(roi_names)
            
             donegetdata = 0;
            
            %Be robust to getdata crashing (which it does a lot recently)
            while donegetdata == 0
                try
                    %Load and interrogate ROI of interest
                    roi = [roi_dir, roi_names{roinum} , '_roi.mat'];
                    rois = maroi('load_cell', roi);
                    [seed_voxels a b c] = getdata(rois{1}, P,'l');
                    
                    %Average across voxels to get a mean seed timecourse
                    seed_timecourse = mean(seed_voxels,2);
                    
                    donegetdata = 1;
                    
                catch
                    
                    disp(['Getdata failed: ' datestr(now)])
                    
                end
            end
            
            [B,BINT,R] = REGRESS(seed_timecourse,regressors);
            seed_data(:,roinum) = R;
            
        end
        
        result(subject) = cca_granger_regress(seed_data',3);
        
        for i = 1:length(roi_names)
            for j = 1:length(roi_names)
                if i~=j
                    
                    texttowrite = [subjects{subject},'   ',[roi_names{i} '_to_' roi_names{j}],'   ',num2str(result(subject).gc(j,i))];
                    dlmwrite(outputfilename,texttowrite,'-append','delimiter','');
                    
                end
            end
        end
        
        
    end
end
            
            