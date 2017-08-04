subjects = {'vc34268'};
tosubjects = {'vc35469'};
%
for subject = 1:length(subjects)
    
    sourcedir = ['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' subjects{subject} '/7112b_fs_LR/'];
    targetdir = ['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' tosubjects{subject} '/7112b_fs_LR/'];
    
    mkdir(targetdir);
    
    files = dir([sourcedir '*.*']);
    
    for file = 3:length(files)
        index = strfind(files(file).name,subjects{subject});
        
        if isempty(index)
            copyfile([sourcedir files(file).name],[targetdir files(file).name]);
        else
            copyfile([sourcedir files(file).name],[targetdir tosubjects{subject} files(file).name(length(subjects{subject})+1:end)])
        end
    end
    
    
    sourcedir = ['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' subjects{subject} '/7112b_fs_LR/fsaverage/'];
    targetdir = ['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' tosubjects{subject} '/7112b_fs_LR/fsaverage/'];
    
    mkdir(targetdir);
    
    files = dir([sourcedir '*.*']);
    
    for file = 3:length(files)
        index = strfind(files(file).name,subjects{subject});
        
        if isempty(index)
            copyfile([sourcedir files(file).name],[targetdir files(file).name]);
        else
            copyfile([sourcedir files(file).name],[targetdir tosubjects{subject} files(file).name(length(subjects{subject})+1:end)])
        end
    end
    
    
    sourcedir = ['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' subjects{subject} '/7112b_fs_LR/fsaverage_LR32k/'];
    targetdir = ['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' tosubjects{subject} '/7112b_fs_LR/fsaverage_LR32k/'];
    
    mkdir(targetdir);
    
    files = dir([sourcedir '*.*']);
    
    for file = 3:length(files)
        index = strfind(files(file).name,subjects{subject});
        
        if isempty(index)
            copyfile([sourcedir files(file).name],[targetdir files(file).name]);
        else
            copyfile([sourcedir files(file).name],[targetdir tosubjects{subject} files(file).name(length(subjects{subject})+1:end)])
        end
    end
        
    
    sourcedir = ['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' subjects{subject} '/7112b_fs_LR/Native/'];
    targetdir = ['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' tosubjects{subject} '/7112b_fs_LR/Native/'];
    
    mkdir(targetdir);
    
    files = dir([sourcedir '*.*']);
    
    for file = 3:length(files)
        index = strfind(files(file).name,subjects{subject});
        
        if isempty(index)
            copyfile([sourcedir files(file).name],[targetdir files(file).name]);
        else
            copyfile([sourcedir files(file).name],[targetdir tosubjects{subject} files(file).name(length(subjects{subject})+1:end)])
        end
    end
    
    
    

end