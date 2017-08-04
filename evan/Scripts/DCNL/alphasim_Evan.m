% alphasim_Evan
%
% Calculates corrected p-values and cluster-thresholds for a given SPM model
% based on the alphasim utility within AFNI.
% 
% A GUI prompts the user to choose an SPM.mat file from an
% already-completed analysis. The script then chooses a random subject from
% that analysis, calculates the smoothness of that data, finds the mask
% that was used for the SPM analysis, and runs the AFNI 3dClustSim utility
% within that mask to simulate the data 10,000 times and establish the
% probability of false clusters appearing.
% 
% A table is spit out listing the correct cluster size to get corrected
% alphas at a number of uncorrected alphas.
%
% You must have SPM8 in your Matlab path, and you must have both FSL and
% AFNI loaded on your computer for this to work. 
%
%Created by E. Gordon January 2011.


%Desired alphas to test significance at.  Leave empty to use defaults.
alphas = [.05 .01 .005];


%Location of afni programs.  Can be left empty if afni is loaded into your path by default
afni_location = '/data/apps/afni/2009_12_31_1431/x86_64/';

%Prompt subject for an SPM.mat file
[spmmatfile, sts] = spm_select(1,'^SPM\.mat$','Select SPM.mat');
if sts==0
    return
end

%Load that file
load(spmmatfile);
groupSPM = SPM;
clear SPM

%delete any previous temporary files
delete([groupSPM.swd '/temp_alphasim*']);

%identify the mask used in the SPM analysis
maskfile = [groupSPM.swd '/mask.hdr'];

%change that mask file to a Nifti
eval(['!fslchfiletype NIFTI ' maskfile ' ' groupSPM.swd '/temp_alphasimmask']);

%Pick a random subject
randsubject = ceil(rand*groupSPM.nscan);

%if this is a second-level analysis
if ~isempty(strfind(groupSPM.xY.P{randsubject},'con_')) || ~isempty(strfind(groupSPM.xY.P{randsubject},'beta_'))
    
    disp('Concatenating scans from random subject...')
    
    secondlevel = 1;
    
    %Find and load the SPM.mat file of the 1st-level analysis from the random subject
    lastslashlocation = find(groupSPM.xY.P{randsubject}=='/');
    load([groupSPM.xY.P{randsubject}(1:lastslashlocation(end)) 'SPM.mat']);
    singlesubSPM = SPM;
    clear SPM
    
    %Set up to merge that subject's functional data into a temporary folder
    FSLstring = ['fslmerge -t ' groupSPM.swd '/temp_alphasimdata'];
    for scan = 1:singlesubSPM.nscan
        
        %Find each functional data file, robust to analyze or nifti format
        stringlength = strfind(singlesubSPM.xY.P(scan,:),'.img');
        if isempty(stringlength)
            stringlength = strfind(singlesubSPM.xY.P(scan,:),'.nii');
        end
        
        %Add the found image to the FSLmerge command
        FSLstring = [FSLstring ' ' singlesubSPM.xY.P(scan,1:(stringlength+3))];
        
    end
    
    %Execute the merge
    eval(['!' FSLstring]);
    
    %And convert the merged functional run to a Nifti
    eval(['!fslchfiletype NIFTI ' groupSPM.swd '/temp_alphasimdata']);
    
else %if this group-level analysis is not a second-level analysis
    disp('Saving scan from random subject...')
    
    %Check whether there is a ',1' at the end of the SPM entry, then convert the random subject's data into a Nifti in the temporary folder
    if strcmp(groupSPM.xY.P{randsubject}(end-1:end),',1')
        eval(['!fslchfiletype NIFTI ' groupSPM.xY.P{randsubject}(1:end-2) ' ' groupSPM.swd '/temp_alphasimdata']);
    else
        eval(['!fslchfiletype NIFTI ' groupSPM.xY.P{randsubject} ' ' groupSPM.swd '/temp_alphasimdata']);
    end
    secondlevel = 0;
end

%Convert the data and mask into AFNI format
disp('Converting scans into afni format...')
eval(['!' afni_location '3dcopy ' groupSPM.swd '/temp_alphasimdata.nii ' groupSPM.swd '/temp_alphasimdata_afni']);
eval(['!' afni_location '3dcopy ' groupSPM.swd '/temp_alphasimmask.nii ' groupSPM.swd '/temp_alphasimmask_afni']);

%Sometimes this gets saved as a +orig file, sometimes as a +tlrc file. Check which one is the case. I have no idea why this changes.
if ~isempty(dir([groupSPM.swd '/temp_alphasimdata_afni+orig*']))
    datainorigform = 1;
else
    datainorigform = 0;
end

%Use AFNI's 3dFWHMx script to calculate the smoothness of the data (detrending if it's a timeseries)
disp('Calculating smoothness...')
if secondlevel == 1
    if datainorigform
        t = evalc(['!' afni_location '3dFWHMx -detrend 20 ' groupSPM.swd '/temp_alphasimdata_afni+orig']);
    else
        t = evalc(['!' afni_location '3dFWHMx -detrend 20 ' groupSPM.swd '/temp_alphasimdata_afni+tlrc']);
    end
else
    if datainorigform
        t = evalc(['!' afni_location '3dFWHMx ' groupSPM.swd '/temp_alphasimdata_afni+orig']);
    else
        t = evalc(['!' afni_location '3dFWHMx ' groupSPM.swd '/temp_alphasimdata_afni+tlrc']);
    end
end

%Load the caluclated smoothness values into variables
for i = 25:-1:10
    allsmooths = str2num(t(end-i:end));
    if ~isempty(allsmooths)
        break
    end
end
xsmooth = num2str(allsmooths(1));
ysmooth = num2str(allsmooths(2));
zsmooth = num2str(allsmooths(3));

%display smoothness
disp(['X smoothness: ' xsmooth '; Y smoothness: ' ysmooth '; Z smoothness: ' zsmooth])

%Run alphasim with that smothness, within the right mask
disp('Running 3DClustSim...')

if isempty(alphas)
    resultsstring = evalc(['!' afni_location '3dClustSim -mask ' groupSPM.swd '/temp_alphasimmask_afni+orig -fwhmxyz ' xsmooth ' ' ysmooth ' ' zsmooth]);
else
    resultsstring = evalc(['!' afni_location '3dClustSim -athr ' num2str(alphas) ' -mask ' groupSPM.swd '/temp_alphasimmask_afni+orig -fwhmxyz ' xsmooth ' ' ysmooth ' ' zsmooth]);
end

%Display the alphasim results
disp(resultsstring)

%Write out the alphasim results to a file within the analysis folder
try delete([groupSPM.swd '/AlphasimOutput.txt']); catch; end
fid = fopen([groupSPM.swd '/AlphasimOutput.txt'],'at');
fprintf(fid,'%s',resultsstring);
fclose(fid);

%Delete the tmeporary folder
delete([groupSPM.swd '/temp_alphasim*']);

    
