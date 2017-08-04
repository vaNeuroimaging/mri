maxworkers = 20;

MAVsubs = dir('/home/data/subjects/MAV*'); MAVsubs = struct2cell(MAVsubs); MAVsubs = MAVsubs(1,:)';
GOBAsubs = dir('/home/data/subjects/Y01*'); GOBAsubs = struct2cell(GOBAsubs); GOBAsubs = GOBAsubs(1,:)';
subjects = [MAVsubs;GOBAsubs];

subcount = length(subjects);

nworkers = min([subcount maxworkers]);

processingpool = parpool(nworkers);

parfor subnum = 1:subcount
    subject = subjects{subnum};
    
    
    %raw
    rawdirs = dir(['/home/data/subjects/' subject '/raw/']);
    for i = 1:length(rawdirs)
        T1files = dir(['/home/data/subjects/' subject '/raw/' rawdirs(i).name '/*T1*.nii.gz']);
        T2files = dir(['/home/data/subjects/' subject '/raw/' rawdirs(i).name '/*T2*.nii.gz']);
        T1files = struct2cell(T1files); T1files = T1files(1,:)';
        T2files = struct2cell(T2files); T2files = T2files(1,:)';
        files = [T1files;T2files];
        
        for j = 1:length(files)
            files{j} = ['/home/data/subjects/' subject '/raw/' rawdirs(i).name '/' files{j}];
            [failed, message] = system(['/home/data/scripts/Processing_Pipeline/mri_deface-v1.22-Linux64 ' files{j} ' /home/data/scripts/Processing_Pipeline/talairach_mixed_with_skull.gca /home/data/scripts/Processing_Pipeline/face.gca ' files{j}]);
            if failed
                disp('Defacing FAILED!')
                disp(message)
            else
                disp(['defaced ' files{j}])
            end
        end
        
    end
    
    
    %preprocessed
    T1files = dir(['/home/data/subjects/' subject '/preprocessed/T1*.nii.gz']);
    T2files = dir(['/home/data/subjects/' subject '/preprocessed/T2*.nii.gz']);
    T1files = struct2cell(T1files); T1files = T1files(1,:)';
    T2files = struct2cell(T2files); T2files = T2files(1,:)';
    files = [T1files;T2files];
        
    for j = 1:length(files)
        files{j} = ['/home/data/subjects/' subject '/preprocessed/' files{j}];
        [failed, message] = system(['/home/data/scripts/Processing_Pipeline/mri_deface-v1.22-Linux64 ' files{j} ' /home/data/scripts/Processing_Pipeline/talairach_mixed_with_skull.gca /home/data/scripts/Processing_Pipeline/face.gca ' files{j}]);
        if failed
            disp('Defacing FAILED!')
            disp(message)
        else
           disp(['defaced ' files{j}])
        end
    end
    
    
    %freesurfer
    freesurferfolder = ['/home/data/subjects/' subject '/freesurfer/'];
    files = {[freesurferfolder 'T1_avg_biascorr.nii.gz'];[freesurferfolder '/mri/T1.mgz'];[freesurferfolder '/mri/orig.mgz'];[freesurferfolder '/mri/orig_nu.mgz'];[freesurferfolder '/mri/orig/001.mgz']};
    for j = 1:length(files)
        if exist(files{j},'file')
            [failed, message] = system(['/home/data/scripts/Processing_Pipeline/mri_deface-v1.22-Linux64 ' files{j} ' /home/data/scripts/Processing_Pipeline/talairach_mixed_with_skull.gca /home/data/scripts/Processing_Pipeline/face.gca ' files{j}]);
            if failed
                disp('Defacing FAILED!')
                disp(message)
            else
                disp(['defaced ' files{j}])
            end
        end
    end
    
    
    %fs_LR
    fs_LRfolder = ['/home/data/subjects/' subject '/fs_LR/'];
    files = {[fs_LRfolder '/MNI/T1_avg_biascorr.nii.gz'];[fs_LRfolder '/NativeVol/T1_avg_biascorr.nii.gz'];[fs_LRfolder '/MNI/T1_avg_biascorr_MNI.nii.gz'];[fs_LRfolder '/NativeVol/T1_avg_biascorr_MNI.nii.gz']};
    for j = 1:length(files)
        if exist(files{j},'file')
            [failed, message] = system(['/home/data/scripts/Processing_Pipeline/mri_deface-v1.22-Linux64 ' files{j} ' /home/data/scripts/Processing_Pipeline/talairach_mixed_with_skull.gca /home/data/scripts/Processing_Pipeline/face.gca ' files{j}]);
            if failed
                disp('Defacing FAILED!')
                disp(message)
            else
                disp(['defaced ' files{j}])
            end
        end
    end
end