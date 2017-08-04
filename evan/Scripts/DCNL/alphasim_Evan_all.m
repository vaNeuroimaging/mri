

[spmmatfile, sts] = spm_select(1,'^SPM\.mat$','Select SPM.mat');
if sts==0
    return
end

load(spmmatfile);

groupSPM = SPM;
clear SPM

delete([groupSPM.swd '/temp_alphasim*']);

maskfile = [groupSPM.swd '/mask.hdr'];

eval(['!fslchfiletype NIFTI ' maskfile ' ' groupSPM.swd '/temp_alphasimmask']);
eval(['!module load afni; 3dcopy ' groupSPM.swd '/temp_alphasimmask.nii ' groupSPM.swd '/temp_alphasimmask_afni']);

try delete([groupSPM.swd '/AlphasimOutput.txt']); catch; end
fid = fopen([groupSPM.swd '/AlphasimOutput.txt'],'at');

for sub = 1:groupSPM.nscan
    
    %if this is a second-level analysis
    if ~isempty(strfind(groupSPM.xY.P{sub},'con_')) || ~isempty(strfind(groupSPM.xY.P{sub},'beta_'))
        
        disp(['Concatenating scans from subject ' num2str(sub) '...'])
        
        secondlevel = 1;
        
        lastslashlocation = find(groupSPM.xY.P{sub}=='/');
        load([groupSPM.xY.P{sub}(1:lastslashlocation(end)) 'SPM.mat']);
        singlesubSPM = SPM;
        clear SPM
        
        FSLstring = ['fslmerge -t ' groupSPM.swd '/temp_alphasimdata'];
        for scan = 1:singlesubSPM.nscan
            stringlength = strfind(singlesubSPM.xY.P(scan,:),'.img');
            FSLstring = [FSLstring ' ' singlesubSPM.xY.P(scan,1:(stringlength+3))];
        end
        eval(['!' FSLstring]);
        eval(['!fslchfiletype NIFTI ' groupSPM.swd '/temp_alphasimdata']);
        
    else %if this group-level analysis is not a second-level analysis
        disp('Saving scan from random subject...')
        eval(['!fslchfiletype NIFTI ' groupSPM.xY.P{randsubject}(1:end-2) ' ' groupSPM.swd '/temp_alphasimdata']);
        secondlevel = 0;
    end

    disp('Converting scans into afni format...')

    t = evalc(['!module load afni; 3dcopy ' groupSPM.swd '/temp_alphasimdata.nii ' groupSPM.swd '/temp_alphasimdata_afni']);
    
    disp('Calculating smoothness...')
    if secondlevel == 1
        t = evalc(['!module load afni; 3dFWHMx -detrend 20 ' groupSPM.swd '/temp_alphasimdata_afni+orig']);
    else
        t = evalc(['!module load afni; 3dFWHMx ' groupSPM.swd '/temp_alphasimdata_afni+orig']);
    end
    xsmooth(sub) = str2num(t(end-25:end-19));
    ysmooth(sub) = str2num(t(end-16:end-10));
    zsmooth(sub) = str2num(t(end-7:end-1));

    disp(['X smoothness: ' num2str(xsmooth(sub)) '; Y smoothness: ' num2str(ysmooth(sub)) '; Z smoothness: ' num2str(zsmooth(sub))])
    disp([]);
    
    fprintf(fid,'%s\n\r',['Subject ' num2str(sub) ' smoothness: X smoothness: ' num2str(xsmooth(sub)) '; Y smoothness: ' num2str(ysmooth(sub)) '; Z smoothness: ' num2str(zsmooth(sub))]);
    
    delete([groupSPM.swd '/temp_alphasimdata*']);
    
end

xsmoothMean = num2str(mean(xsmooth));
ysmoothMean = num2str(mean(ysmooth));
zsmoothMean = num2str(mean(zsmooth));

disp('Running alphasim...')

resultsstring = evalc(['!module load afni; 3dClustSim -mask ' groupSPM.swd '/temp_alphasimmask_afni+orig -fwhmxyz ' xsmoothMean ' ' ysmoothMean ' ' zsmoothMean]);
disp(resultsstring)


fprintf(fid,'%s',resultsstring);
fclose(fid);

delete([groupSPM.swd '/temp_alphasim*']);

    
