function copyfiles_frommedia_DART(targetstr)

if ~exist('targetstr','var')
    targetstr = 'MAV';
end

targetdirroot = '/home/data/subjects/';

overwrite = 0;

media_root = '/media/HERMANOSBAC/NP985/';

possible_prefixes = {'MAV','vc3'};

deface_file_strings = {'prage','T2'};

one_session_notinname = true;

origdir = pwd;

sourcedirs = cell(1,0);
medianames = dir([media_root '/*']);
for medianum = 3:length(medianames)
    if strcmp(medianames(medianum).name(1:length(targetstr)),targetstr)% && exist([media_root '/' medianames(medianum).name '/DICOM'],'dir')
        thissourcedir = [media_root '/' medianames(medianum).name '/'];
        sourcedirs(end+1) = {thissourcedir};
    else
        subdirs = dir([media_root '/' medianames(medianum).name '/']);
        for subdirnum = 1:length(subdirs);
            if (length(subdirs(subdirnum).name) >= length(targetstr)) && (strcmp(subdirs(subdirnum).name(1:length(targetstr)),targetstr))% && exist([media_root '/' medianames(medianum).name '/' subdirs(subdirnum).name '/DICOM'],'dir'))
                thissourcedir = [media_root '/' medianames(medianum).name '/' subdirs(subdirnum).name '/'];
                sourcedirs(end+1) = {thissourcedir};
            end
        end
    end
end

for sourcedirnum = 1:length(sourcedirs)
    conversionworked = 0;
    this_sourcedir = sourcedirs{sourcedirnum};
    
    tokens = tokenize(this_sourcedir,'/');
    for i = length(tokens):-1:1
        if ~isempty(tokens{i})
            sessname = tokens{i};
            break
        end
    end
    
    for i = 1:length(possible_prefixes)
        if strcmp(sessname(1:length(possible_prefixes{i})),possible_prefixes{i})
            prefix = possible_prefixes{i};
            break
        end
    end
    
    thissub = sessname;
    sessnum = '1';
    
    thistargetdir = [targetdirroot thissub '/raw/' thissub '-' sessnum '/'];
    
    %disp(['Copying from ' this_sourcedir ' to ' thistargetdir ', converting, and defacing data'])
    if exist(thistargetdir,'dir') && (length(dir(thistargetdir)) > 2)
        %disp([thistargetdir 'DICOM/ already exists!'])
        if overwrite
            disp(['Converting data from ' this_sourcedir ' to ' thistargetdir ', overwriting previous data'])
            [~,~,~] = rmdir(thistargetdir,'s');
            mkdir(thistargetdir);
            
            [failed, message] = system(['dcm2nii -o ' thistargetdir ' ' this_sourcedir]);
            if failed
                disp('Dicom conversion FAILED!')
                disp(message)
            else
                conversionworked = 1;
            end
                
        end
    else
        disp(['Converting data from ' this_sourcedir ' to ' thistargetdir ])
        mkdir(thistargetdir);
        

        [failed, message] = system(['dcm2nii -o ' thistargetdir ' ' this_sourcedir]);
        if failed
            disp('Dicom conversion FAILED!')
            disp(message)
        else
            conversionworked = 1;
        end
    end
    
%     if conversionworked
%         
%         for defacestr = 1:length(deface_file_strings)
%             deface_files = dir([thistargetdir '/*' deface_file_strings{defacestr} '*.nii.gz']);
%             for filenum = 1:length(deface_files)
%                 cd(thistargetdir)
%                 this_file = [thistargetdir '/' deface_files(filenum).name];
%                 [failed, message] = system(['/home/data/scripts/Processing_Pipeline/mri_deface-v1.22-Linux64 ' this_file ' /home/data/scripts/Processing_Pipeline/talairach_mixed_with_skull.gca /home/data/scripts/Processing_Pipeline/face.gca ' this_file]);
%                 if failed
%                     disp('Defacing FAILED!')
%                     disp(message)
%                 end
%             end
%         end
%         
%         
%     end
    
    [failed, message] = system(['chmod ga+r ' thistargetdir '/*']);
    
end

cd(origdir)
