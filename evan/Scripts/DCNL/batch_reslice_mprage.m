subjects = {'415'};

rootdir='/fmri/data3/Evan/Gene-Rest-Nback/Data/';

for subnum = 1:length(subjects)

    subject = subjects{subnum};
    datadir = [rootdir subject '/Struct/FIRST_Hip_seg/']; 
%    copyfile([datadir 'MRPAGE.nii'],'temp.nii');
    
%     rmdir(datadir,'s');
%     mkdir(datadir);
%     
%     copyfile('temp.nii',[datadir 'MPRAGE.nii']);
    
    reslice_nii([datadir 'MPRAGE.nii'],[datadir 'MPRAGE_reslice_step2.nii'],[1 1 1], [1], [], [1], [], []);
 
    %eval(['!flirt -dof 6 -in ' datadir 'MPRAGE.nii -ref ' datadir 'MPRAGE_reslice_step2.nii -out ' datadir 'MPRAGE_reslice_step3']);
    
    
end  

