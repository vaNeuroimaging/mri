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


finalgroupSPM = groupSPM;



%Prompt subject for an SPM.mat file
[spmmatfile, sts] = spm_select(1,'^SPM\.mat$','Select an SPM.mat with the same subjects, run on 1st-level SPMs');
if sts==0
    return
end

%Load that file
load(spmmatfile);
groupSPM = SPM;
clear SPM






%if this is a second-level analysis
    
    disp('Concatenating scans from random subject...')
    
    secondlevel = 1;
    
    %Find and load the SPM.mat file of the 1st-level analysis from the random subject
    lastslashlocation = find(groupSPM.xY.P{randsubject}=='/');
    load([groupSPM.xY.P{randsubject}(1:lastslashlocation(end)) 'SPM.mat']);
    singlesubSPM = SPM;
    clear SPM
    
    %Set up to merge that subject's functional data into a temporary folder
    FSLstring = ['fslmerge -t ' finalgroupSPM.swd '/temp_alphasimdata'];
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
    eval(['!fslchfiletype NIFTI ' finalgroupSPM.swd '/temp_alphasimdata']);
    


%Convert the data and mask into AFNI format
disp('Converting scans into afni format...')
eval(['!' afni_location '3dcopy ' finalgroupSPM.swd '/temp_alphasimdata.nii ' finalgroupSPM.swd '/temp_alphasimdata_afni']);
eval(['!' afni_location '3dcopy ' finalgroupSPM.swd '/temp_alphasimmask.nii ' finalgroupSPM.swd '/temp_alphasimmask_afni']);

%Sometimes this gets saved as a +orig file, sometimes as a +tlrc file. Check which one is the case. I have no idea why this changes.
if ~isempty(dir([finalgroupSPM.swd '/temp_alphasimdata_afni+orig*']))
    datainorigform = 1;
else
    datainorigform = 0;
end

%Use AFNI's 3dFWHMx script to calculate the smoothness of the data (detrending if it's a timeseries)
disp('Calculating smoothness...')
if secondlevel == 1
    if datainorigform
        t = evalc(['!' afni_location '3dFWHMx -detrend 20 ' finalgroupSPM.swd '/temp_alphasimdata_afni+orig']);
    else
        t = evalc(['!' afni_location '3dFWHMx -detrend 20 ' finalgroupSPM.swd '/temp_alphasimdata_afni+tlrc']);
    end
else
    if datainorigform
        t = evalc(['!' afni_location '3dFWHMx ' finalgroupSPM.swd '/temp_alphasimdata_afni+orig']);
    else
        t = evalc(['!' afni_location '3dFWHMx ' finalgroupSPM.swd '/temp_alphasimdata_afni+tlrc']);
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
disp('Running alphasim...')
resultsstring = evalc(['!' afni_location '3dClustSim -mask ' finalgroupSPM.swd '/temp_alphasimmask_afni+orig -fwhmxyz ' xsmooth ' ' ysmooth ' ' zsmooth]);

%Display the alphasim results
disp(resultsstring)

%Write out the alphasim results to a file within the analysis folder
try delete([finalgroupSPM.swd '/AlphasimOutput.txt']); catch; end
fid = fopen([finalgroupSPM.swd '/AlphasimOutput.txt'],'at');
fprintf(fid,'%s',resultsstring);
fclose(fid);

%Delete the tmeporary folder
delete([finalgroupSPM.swd '/temp_alphasim*']);