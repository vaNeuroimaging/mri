%dual_regression.m
%
%Uses the output from a group ICA analysis to determine individual subject versions
%of the identified components within each subject's functional data.
%
%This is performed by spatially regressing each timepoint against the group
%components, and then temporally regressing the individual's whole functional
%timecourse against the resulting timecourses of beta weights.
%
%Importantly, the ordering of these resulting individual images (the
%con_000# images SPM spits out) are in the same order as the group
%components (so e.g., con_0004 matches IC#4 on the FSL report).
%
%You must have SPM8, Marsbar, and the FSL_from_matlab script in your path for this script to work
%
%Created by E. Gordon July 09.


%User-specified stuff (more down on line 99)
%--------------------------------------------------------------------------
subjects = {'101','102','113','118','120','122','125','127','132','138','147','150','151','154','156','159','160','161','162','166','172','181','182','187','202','207','211','214','215','221','225','229','232','233','242','250','254','255','272','274', '112' '126' '133' '137' '208' '222' '227' '251' '253' '256' '258' '260' '261' '264' '270' '279' '281' '283' '292' '301' '307' '322' '327' '247','295','297','300','305','309','334','339','340','343','359','362','383','395','396'};
%

%location of Group ICA results (i.e. the melodic_IC.nii.gz file)
componentpath = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ICA/TensorICA_76DAT.gica/groupmelodic.ica/'];

%Prefix of preprocessed functional images
prefix = 'Timepoint';

%location of SPM template file
SPM_analysistemplatefile = '/fmri/data3/Evan/Gene-Rest-Nback/Scripts/DualRegression_SPM8_template.mat';

%TR in seconds
TR = 2;

%Microtime resolution and onset
MicroRes = 16;
MicroOnset = 1;

conditions = {'1Back','2Back','3Back'};

%--------------------------------------------------------------------------

%If the Marbar wholebrain mask hasn't been created, make it
if isempty(dir([componentpath '../mask_roi.mat']));
    
    %Convert brain mask image crated by Melodic to Nifti format
    FSL_from_matlab(['fslchfiletype NIFTI ' componentpath '../mask.nii.gz']);
    
    %Build Marsbar mask from the brain mask image
    imgname = [componentpath '../mask.nii'];
    outputname = [componentpath '../mask_roi.mat'];
    o = []; d = [];
    [p f e] = fileparts(imgname);
    binf = 1;
    func = 'img >= 1';
    d = f; l = f;
    d = [d ' func: ' func]; l = [l '_f_' func];
    d = [d ' - binarized']; l = [l '_bin'];
    o = maroi_image(struct('vol', spm_vol(imgname), 'binarize', binf, 'func', func));
    o = maroi_matrix(o);
    o = descrip(o,d);  o = label(o,l);
    varargout = {saveroi(o, outputname)};
end

%load the wholebrain mask
%This is needed so that the variance in voxels outside the brain isn't entered into the regression model.
brainroi = maroi([componentpath '../mask_roi.mat']);


%If the component files haven't been split out and converted, do that
if isempty(dir([componentpath 'GroupComponent*.nii']));
    
    %Split melodic_IC file
    FSL_from_matlab(['fslsplit ' componentpath 'melodic_IC.nii.gz ' componentpath 'GroupComponent -t']);
    
    %Convert components to Nifti format
    componentnames = dir([componentpath 'GroupComponent*.nii.gz']);
    for i=1:length(componentnames);
        FSL_from_matlab(['fslchfiletype NIFTI ' componentpath componentnames(i).name]);
    end
end

%Get the component names and the number of components
componentnames = dir([componentpath 'GroupComponent*.nii']);
numcomponents = length(componentnames);

%Put names of components, padded with spaces, into a matrix so that marsbar can run them all at once in the same space as a functional image
P = [];
for component = 1:numcomponents
    P = [P;sprintf('%-150s',[componentpath componentnames(component).name])];
end

warning off

%subject loop
for subnum = 1:length(subjects)
    subj = subjects{subnum};
    
    
    for condition = 1:length(conditions)
        conditionname = conditions{condition};
        
        %User-specified stuff
        %--------------------------------------------------------------------------
        Rest_datafolder = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subj '/SPM8/Nback/' conditionname '_Temp/'];
        Outputfolder = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/SPM8/DualRegression_' conditionname '_ICA76/'];
        %----------------------------------------------------------------------
        
        try rmdir(Rest_datafolder,'s'); catch; end
        mkdir([Rest_datafolder]);
        FSL_from_matlab(['fslsplit ' Rest_datafolder '../' conditionname 'only_smoothed.nii.gz ' Rest_datafolder prefix ' -t']);
        files = dir([Rest_datafolder prefix '*.nii.gz']);
        if length(files)==0
            error('files not created')
        end
        for i = 1:length(files)
            FSL_from_matlab(['fslchfiletype ANALYZE ' Rest_datafolder files(i).name]);
        end
        
        
        
        %1st regression: spatial
        %-----------------------------------------------------------------------
        
        %find functional data for this subject
        Imagesinfolder = dir([Rest_datafolder prefix '*.img']);
        num_timepoints = length(Imagesinfolder);
        
        %For the first subject, get the data within the wholebrain ROI for all the components, interpolated into same space as the functional images.
        %Then save that 3D data as a vector within a Designmatrix for each component, for later regression.
        %Note that you only need to do this for the first subject unless your subjects are in different spaces for some reason.
        if subnum == 1
            charactermat = [sprintf('%-150s',[Rest_datafolder Imagesinfolder(1).name]);P];
            Y = getdata(brainroi,charactermat,'l');
            Designmatrix = [Y(2:end,:)',ones(size(Y,2),1)];
        end
        
        timecourse = zeros(num_timepoints,1);
        
        disp(['';''])
        
        %For each timepoint:
        for time = 1:num_timepoints
            
            %print progress to screen
            string{time} = ['Spatial regression: Subject ' subj ', condition ' conditionname ', timepoint ' num2str(time)];
            
            if time==1; fprintf('%s',string{time}); else; fprintf([repmat('\b',1,length(string{time-1})) '%s'],string{time}); end
            
            %Get subject's data within the ROI at that timepoint (the data is a vector)
            subjectdata = getdata(brainroi,[Rest_datafolder Imagesinfolder(time).name],'l');
            
            %Regress the data against all the components
            betavals = Designmatrix\(subjectdata');
            
            %Save the resulting beta values for that timepoint
            %These betas represent the strength/explanatory power of each group component at this timepoint
            timecourse(time,1:(length(betavals)-1)) = betavals(1:end-1)';
        end
        %-----------------------------------------------------------------------
        
        
        %2nd regression: temporal
        %-----------------------------------------------------------------------
        
        
        %load the template for the standard SPM temporal regression
        load(SPM_analysistemplatefile)
        
        %Feed all the subject-specific info into the template
        for scan = 1:num_timepoints
            matlabbatch{1}.spm.stats.fmri_spec.sess.scans{scan} = [Rest_datafolder Imagesinfolder(scan).name ',1'];
        end
        
        matlabbatch{1}.spm.stats.fmri_spec.timing.RT = TR;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t = MicroRes;
        matlabbatch{1}.spm.stats.fmri_spec.timing.fmri_t0 = MicroOnset;
        
        matlabbatch{1}.spm.stats.fmri_spec.dir{1} = Outputfolder;
        
        %make/clear directory for SPM data to be saved
        try rmdir(matlabbatch{1}.spm.stats.fmri_spec.dir{1},'s');catch;end
        mkdir(matlabbatch{1}.spm.stats.fmri_spec.dir{1})
        
        
        %Put in the beta value timecourses as temporal regressors
        for component = 1:numcomponents
            matlabbatch{1}.spm.stats.fmri_spec.sess.regress(component).val = timecourse(:,component);
            matlabbatch{1}.spm.stats.fmri_spec.sess.regress(component).name = ['Component ' component];
        end
        
        %save the batch file for the temporal regression
        nameofSPMjob = [matlabbatch{1}.spm.stats.fmri_spec.dir{1}(1:end-1) '.mat'];
        try delete(nameofSPMjob); catch; end
        save(nameofSPMjob, 'matlabbatch');
        
        %Call SPM5 and run the temporal regression
        spm_jobman('run',nameofSPMjob)
        
        clear matlabbatch
        %-----------------------------------------------------------------------
        
        rmdir(Rest_datafolder,'s')
        
    end
end


