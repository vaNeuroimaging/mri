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

subjects = {'208','222','258','227','253','256','270','279','292','301','305','309','334','343','359','362','374','383','395','415','416','420'};

%Header of functional data to interrogate
header3d = 'sw*.img';


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
randsubject = ceil(rand*length(subjects));
subjsid = subjects{randsubject};


disp('Concatenating scans from random subject...')

data_dir = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/', subjsid, '/SPM8/FirstRest/'];

datafiles3D = [data_dir header3d];
%Locate images
imgfiles = dir(datafiles3D);
if isempty(imgfiles)
    datafiles3D = [data_dir header3d(1:end-4) '.nii'];
    %Locate images
    imgfiles = dir(datafiles3D);
end


%Set up to merge that subject's functional data into a temporary folder
FSLstring = ['fslmerge -t ' groupSPM.swd '/temp_alphasimdata'];
for scan = 1:size(imgfiles, 1)
    
    %Add the found image to the FSLmerge command
    FSLstring = [FSLstring ' ' [data_dir imgfiles(scan).name]];
    
end

%Execute the merge
eval(['!' FSLstring]);

%And convert the merged functional run to a Nifti
eval(['!fslchfiletype NIFTI ' groupSPM.swd '/temp_alphasimdata']);



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
if datainorigform
    t = evalc(['!' afni_location '3dFWHMx -detrend 20 ' groupSPM.swd '/temp_alphasimdata_afni+orig']);
else
    t = evalc(['!' afni_location '3dFWHMx -detrend 20 ' groupSPM.swd '/temp_alphasimdata_afni+tlrc']);
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
disp('Running alphasim...')
resultsstring = evalc(['!' afni_location '3dClustSim -mask ' groupSPM.swd '/temp_alphasimmask_afni+orig -fwhmxyz ' xsmooth ' ' ysmooth ' ' zsmooth]);

%Display the alphasim results
disp(resultsstring)

%Write out the alphasim results to a file within the analysis folder
try delete([groupSPM.swd '/AlphasimOutput.txt']); catch; end
fid = fopen([groupSPM.swd '/AlphasimOutput.txt'],'at');
fprintf(fid,'%s',resultsstring);
fclose(fid);

%Delete the tmeporary folder
delete([groupSPM.swd '/temp_alphasim*']);

    
