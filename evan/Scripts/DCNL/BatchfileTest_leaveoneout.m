load Batch.mat

scans1 = matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1;
scans2 = matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2;
matlabbatch{1}.spm.stats.factorial_design.dir{1} = '/data/projects/ASD/Data/Data_Analysis/EmoDotProbe/EmoDotProbe_Group_Results/Con_vs_ASD_Angry_vs_Neut_leaveoneout/';

filename = '/data/projects/ASD/Data/Data_Analysis/EmoDotProbe/EmoDotProbe_Group_Results/Con_vs_ASD_Angry_vs_Neut/testresults.txt';

delete(filename);
fid = fopen(filename,'at'); %open the output file for writing
fprintf(fid,'%s\n\r\','Test Results:'); %write the output file header
fclose(fid);
dlmwrite(filename,' ','-append');


for subject = 1:length(scans1)
    
    scanstoadd = [scans1(1:subject-1); scans1(subject+1:end)];
    matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1 = scanstoadd;
    
    try rmdir(matlabbatch{1}.spm.stats.factorial_design.dir{1},'s'); catch; end
    mkdir(matlabbatch{1}.spm.stats.factorial_design.dir{1});

    save([matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat'], 'matlabbatch');
    
    try
        spm_jobman('run',[matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat']);
        dlmwrite(filename,['Worked without ' scans1{subject}],'-append','delimiter','');
    catch
        dlmwrite(filename,['Failed without ' scans1{subject}],'-append','delimiter','');
    end
    close all

end

matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1 = scans1;

for subject = 1:length(scans2)
    
    scanstoadd = [scans2(1:subject-1); scans2(subject+1:end)];
    matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2 = scanstoadd;
    
    try rmdir(matlabbatch{1}.spm.stats.factorial_design.dir{1},'s'); catch; end
    mkdir(matlabbatch{1}.spm.stats.factorial_design.dir{1});

    save([matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat'], 'matlabbatch');
    
    try
        spm_jobman('run',[matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat']);
        dlmwrite(filename,['Worked without ' scans2{subject}],'-append','delimiter','');
    catch
        dlmwrite(filename,['Failed without ' scans2{subject}],'-append','delimiter','');
    end
    close all

end