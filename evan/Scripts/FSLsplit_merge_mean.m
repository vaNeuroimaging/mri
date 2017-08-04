subjectlist = '/data/cn4/evan/Task_parcellation/HCPSubjects/SubjectList3.txt';

subjects = textread(subjectlist,'%s');

betalist = '/data/cn4/evan/Task_parcellation/HCPSubjects/MainTaskContrasts_torun.txt';

[ign betanums ign ign ign ign] = textread(betalist,'%s%s%s%s%s%s');

datafolder = '/data/cn4/evan/Task_parcellation/HCPSubjects/';



thisdir = pwd;

cd(datafolder)

%allmergestr = 'fslmerge -t AllSubAllbetas.nii.gz';

for s = 1:length(subjects)
    
%    system(['fslsplit ' subjects{s} '_mergeFile.nii.gz ' subjects{s} '_mergeFile']);

 %   allmergestr = [allmergestr ' ' subjects{s} '_mergeFile.nii.gz'];
    
end

%Finalmergestr = 'fslmerge -t MeanSub_allbetas.nii.gz';



% for beta = 2:length(betanums)
%     
%     betamergestring = ['fslmerge -t Allsubs_beta' num2str(betanums{beta}) '.nii.gz'];
%     
%     for s = 1:length(subjects)
%         
%         betamergestring = [betamergestring ' ' subjects{s} '_mergeFile'  sprintf('%04d',str2num(betanums{beta})-1) '.nii.gz'];
%         
%         allmergestr = [allmergestr ' ' subjects{s} '_mergeFile'  sprintf('%04d',str2num(betanums{beta})-1) '.nii.gz'];
%         
%     end
%     
%     system(betamergestring)
%     
%     system(['fslmaths Allsubs_beta' num2str(betanums{beta}) '.nii.gz -ing 1 -Tmean Mean_beta' num2str(betanums{beta}) '.nii.gz'])
%     
%     Finalmergestr = [Finalmergestr ' Mean_beta' num2str(betanums{beta}) '.nii.gz'];
%     
% end

%system(allmergestr)

%system(Finalmergestr)

%system('fslchfiletype NIFTI MeanSub_allbetas.nii.gz');

cd(thisdir)