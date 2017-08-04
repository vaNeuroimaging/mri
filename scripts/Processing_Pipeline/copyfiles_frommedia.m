function copyfiles_frommedia(targetstrings,usedate)
%copyfiles_frommedia(targetstr,usedate)

% if ~exist('targetstrings','var')
%     targetstrings = 'MAV';
% end

targetdirroot = '/home/data/subjects/';

overwrite = 0;

media_root = '/media/';

possible_prefixes = {'MAV','MBG','ROBI','Movie'};

deface_file_strings = {'T1','T2'};

one_session_notinname = false;

origdir = pwd;


if ~iscell(targetstrings)
    targetstrings = {targetstrings};
end

for stringnum = 1:length(targetstrings)
    targetstr = targetstrings{stringnum};
    
    sourcedirs = cell(1,0);
    medianames = dir([media_root '/*']);
    for medianum = 3:length(medianames)
        if (length(medianames) >= length(targetstr)) && strcmp(medianames(medianum).name(1:length(targetstr)),targetstr) && exist([media_root '/' medianames(medianum).name '/DICOM'],'dir')
            thissourcedir = [media_root '/' medianames(medianum).name '/'];
            sourcedirs(end+1) = {thissourcedir};
        else
            subdirs = dir([media_root '/' medianames(medianum).name '/']);
            for subdirnum = 1:length(subdirs);
                thissourcedir = [media_root '/' medianames(medianum).name '/' subdirs(subdirnum).name '/'];
                if exist('usedate','var') && (length(subdirs(subdirnum).name)>2) &&  ~strcmp(subdirs(subdirnum).name(1),'$') && (datetime(datevec(targetstr)) < datetime(datevec(subdirs(subdirnum).date)))
                    sourcedirs(end+1) = {thissourcedir};
                
                elseif ~exist('usedate','var') && (length(subdirs(subdirnum).name) >= length(targetstr)) && (strcmp(subdirs(subdirnum).name(1:length(targetstr)),targetstr) && exist([media_root '/' medianames(medianum).name '/' subdirs(subdirnum).name '/DICOM'],'dir'))
                    
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
        
        thissub = sessname(1:length(prefix)+3);
        if strcmp(sessname(length(prefix)+4),'-')
            sessname(length(prefix)+4) = [];
        end
        if one_session_notinname
            sessnum = 1;
        else
            sessnum = num2str(str2num(sessname(length(prefix)+4 : end)));
        end
        thistargetdir = [targetdirroot thissub '/raw/' thissub '-' sessnum '/'];
        
        %disp(['Copying from ' this_sourcedir ' to ' thistargetdir ', converting, and defacing data'])
        if exist(thistargetdir,'dir') && (length(dir(thistargetdir)) > 2)
            %disp([thistargetdir 'DICOM/ already exists!'])
            if overwrite
                disp(['Converting data from ' this_sourcedir ' to ' thistargetdir ', overwriting previous data, and defacing high-res data'])
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
            disp(['Converting data from ' this_sourcedir ' to ' thistargetdir ' and defacing high-res data'])
            mkdir(thistargetdir);
            
            
            [failed, message] = system(['dcm2nii -o ' thistargetdir ' ' this_sourcedir]);
            if failed
                disp('Dicom conversion FAILED!')
                disp(message)
            else
                conversionworked = 1;
            end
        end
        
        if conversionworked
            
            for defacestr = 1:length(deface_file_strings)
                deface_files = dir([thistargetdir '/*' deface_file_strings{defacestr} '*.nii.gz']);
                for filenum = 1:length(deface_files)
                    cd(thistargetdir)
                    this_file = [thistargetdir '/' deface_files(filenum).name];
                    [failed, message] = system(['/home/data/scripts/Processing_Pipeline/mri_deface-v1.22-Linux64 ' this_file ' /home/data/scripts/Processing_Pipeline/talairach_mixed_with_skull.gca /home/data/scripts/Processing_Pipeline/face.gca ' this_file]);
                    if failed
                        disp('Defacing FAILED!')
                        disp(message)
                    end
                end
            end
            
            
        end
        
        [failed, message] = system(['chmod ga+r ' thistargetdir '/*']);
        
    end
end

cd(origdir)
