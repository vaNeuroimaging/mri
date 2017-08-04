function make_slicetiming_files_jdp(subjects)

slicetime_info_file = '@/home/data/scripts/Resources/slice_sec.txt';
TR = 2.2;
frameskip = 5;


if iscell(subjects)

elseif exist(subjects,'file')==2
    
    subjects = textread(subjects,'%s');
    
elseif ischar(subjects)
    
    subjects = {subjects};
    
end

for s = 1:length(subjects)
    
    folder = ['/home/data/subjects/' subjects{s} '/preproc_JDP/'];
    restfiles = dir([folder '*est*.nii.gz']);
    for f = 1:length(restfiles)
        fileprefix = [folder restfiles(f).name(1:end-7)];
        dlmwrite([fileprefix '.SLICETIME'],slicetime_info_file,'delimiter','');
        dlmwrite([fileprefix '.TR'],TR,'delimiter','');
        dlmwrite([fileprefix '.SKIPTRS'],frameskip,'delimiter','');
    end
end