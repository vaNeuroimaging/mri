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
subjects = {'112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','247','295','297','300','305','309','334','339','340','343','359','362','383','395','396','374'};
%'101','102','113','118','120','122','125','127','132','138','147','150','151','154','156','159','160','161','162','166','172','181','182','187','202','207','211','214','215','221','225','229','232','233','242','250','254','255','272','274',

%Names of runs to interrogate
runs = {'FirstRest'};
%{'Rest','Nback'};

%Header of functional data to interrogate
header3d = 'fsw*.nii';
header4d = 'filt_smoothed.nii.gz';

%location of seed ROIs
roi_dir = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/diMartino/'];

%Names of seed ROIs to interrogate. Do not include the _roi.mat suffix
roi_names = {'DC_bilat','VSi_bilat','VRP_bilat','DCP_bilat'};
%{'aDMN-vmPFC','pDMN-Prec','rFPC-rdlPFC','lFPC-ldlPFC','Sal-ramfg','Striatum-l'};
%{'rFPC-rdlPFC','sal-SMA','DMN-Prec','str-r','DMN-vmPFC'};
%{'DC_bilat','VSi_bilat','VRP_bilat','DCP_bilat','PCC_DeLuca'};

%Timepoints within functional run to use in correlation.  If all, leave blank.
timepointstouse = [];

%Regressions to perform (1=yes, 0=no)
regress_WM_and_CSF = 1;
regress_motionparams = 1;
regress_globalsignal = 1;
regress_task = 1;

%Whole-brian binary mask image (to keep from doing correlations outside the brain).
brainmaskname = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/Wholebrain_limited.hdr';

%END USER INPUT
%--------------------------------------------------------------------------

warning off

brainmask = load_nii(brainmaskname);

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
    
    outputfolder = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/', subjsid, '/VoxelwiseConnectivity/'];
    
    
    for run = 1:length(runs)
        
        %USER INPUT (and still more below)
        %--------------------------------------------------------------------------
        
        %location of functional data for this subject
        data_dir = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/', subjsid, '/SPM8/' runs{run} '/'];
        
        %Location of CSF, WM, and wholebrain ROIs for this subject.  If not being used, you can leave these blank.
        CSFROI = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/Ind_ROIs/' subjsid '_CSF_roi.mat'];
        WMROI = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/Ind_ROIs/' subjsid '_WM_roi.mat'];
        BrainROI = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/Wholebrain_limited_roi.mat'];
        
        %folder (to be created) in which connectivity results for this subject/ROI will be written
        outputfolder = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/', subjsid, '/VoxelwiseConnectivity/'];
        
        %END USER INPUT
        %--------------------------------------------------------------------------
        
        
        
        %Restrict Nback run to timepoints within the task
        if strcmp(runs{run},'Nback')
            if isempty(timepointstouse)
                finaltimepointstouse = nbacktimepointstouse;
            else
                finaltimepointstouse = intersection(timepointstouse,nbacktimepointstouse);
            end
        else
            finaltimepointstouse = timepointstouse;
        end
        
        
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
        
        %Loop through seed ROIs
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
                    seed_timecourses(:,roinum) = mean(seed_voxels,2);
                    
                    donegetdata = 1;
                    
                catch
                    
                    disp(['Getdata failed: ' datestr(now)])
                    
                end
            end
        end
        
        [Q,~] = qr_GS([ones(size(seed_timecourses,1),1) seed_timecourses regressors]);
        Q = Q(:,2:end);
        
        
        
        for roinum = 1:length(roi_names)
            
            %USER INPUT
            %--------------------------------------------------------------------------
            
            %Name of image to be written
            outputimgname = [roi_names{roinum} '_' runs{run} '_orth.nii'];
            
            %END USER INPUT
            %--------------------------------------------------------------------------
            
            
            
            disp([subjsid ', ' runs{run} ': ' roi_names{roinum}])
            
            orthseed_timecourse = Q(:,roinum);
            orthregressors = [Q(:,1:roinum-1) Q(:,roinum+1:end)];
            
            %Perform a partial correlation between the seed timecourse and every other voxel in the brain, partialling out effects of the specified regressors
            
            %Check whether you have access to the graphics card (which can do this much faster)
            %             try a = glogical(1);
            %                 R = jacket_partialcorr(seed_timecourse, reshapeddata, regressors);
            %                 clear a
            %             %if not, do it the old way
            %             catch
            R = partialcorr(orthseed_timecourse, reshapeddata, orthregressors);
            %             end
            
            
            %Convert the Pearson's R values to Fisher's Z
            Fishervals = .5*(log(1+R)-log(1-R));
            
            %Reshape the time x voxels Z matrix back into a 3D image of the same dimensionality as the original data
            Fisherimg = reshape(Fishervals,size(datamatrix,1),size(datamatrix,2),size(datamatrix,3));
            
            %Clear a bunch of memory-hogging variables
            clear R Fishervals seed_voxels seed_timecourse residual residual_timecourses mean_residual_timecourse
            
            %Set the Z value of every voxel outside the brain to be zero
            Fisherimg(find(brainmask.img==0)) = 0;
            
            %Create the output directory
            try mkdir(outputfolder); catch; end
            
            %Delete any image of the same name already in there
            try delete([outputfolder outputimgname]); catch; end
            
            %Make the output nifti
            outputimage = make_nii(Fisherimg, data.hdr.dime.pixdim(2:4), data.hdr.hist.originator(1:3));
            
            %Save the output nifti
            save_nii(outputimage,[outputfolder outputimgname]);
            
            clear Fisherimg outputimage
            
        end
        
        clear seed_timecourses orthseed_timecourse orthregressors timecourse residual_timecourse residuals P datamatrix imgfiles fullimgfiles motionparams gatingtimes residuals data outputimage greshapeddata
        
    end
end



