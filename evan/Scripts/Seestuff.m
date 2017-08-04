subjects = {'vc34096' 'vc33457' 'vc34125' 'vc34126' 'vc34128' 'vc34140' 'vc34141' 'vc34198' 'vc34199' 'vc34200' 'vc34201' 'vc34220' 'vc33378' 'vc35175' 'vc34252' 'vc33775_2' 'vc35469' 'vc34306' 'vc34307' 'vc34308' 'vc34330' 'vc34331' 'vc34401' 'vc34402' 'vc34403' 'vc34404' 'vc33769' 'vc34408'};

for subject = 1:length(subjects)
try ls(['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' subjects{subject} '/7112b_fs_LR/Native/']);
%    disp([subjects{subject} ' exists'])
catch
end
    
end