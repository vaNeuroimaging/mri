function check_motion_firstrun(subject,runs)

subdirroot = '/home/data/subjects/';

functional_sequences = {'restingstate'};

FDthresh = .2;
brainradius = 50;
contiguoustimepoints = 5;

if ~exist('runnum')
    runs = 1;
end

if strcmp(runs,'all')
    runs = 1:5;
end

for seq = 1:length(functional_sequences)
    
    FD = [];
    
    for runnum = runs
        
        subject_rawdir = [subdirroot subject '/raw/' subject '-' num2str(runnum) '/'];
        
        if ~exist(subject_rawdir,'dir') || length(dir(subject_rawdir)) < 2
            disp([subject_rawdir ' doesn''t exist; looking for the session on inserted media'])
            
            failed = copyfiles_frommedia_firstsess(subject,subject_rawdir);
            if failed
                error('copying or conversion failure')
            end
            
        end
        
        
        thisseq_files = dir([subject_rawdir '*' functional_sequences{seq} '*.nii.gz']);
        disp(['Found ' num2str(length(thisseq_files)) ' ' functional_sequences{seq} ' runs in ' subject_rawdir '. Getting motion estimates.'])
        
        for f = 1:length(thisseq_files);
            
            data = [subject_rawdir thisseq_files(f).name];
            [fail, result] = system(['mcflirt -in ' data ' -refvol 0 -plots']);
            if fail
                error(result)
            end
            
            paramsfile = [data(1:end-7) '_mcf.par'];
            
            thisrun_params = load(paramsfile);
            thisrun_rot = thisrun_params(:,1:3);
            thisrun_rot_mm = thisrun_rot * brainradius;
            thisrun_params(:,1:3) = thisrun_rot_mm;
            thisrun_params_delta = [zeros(1,6) ; [thisrun_params(2:end,:) - thisrun_params(1:end-1,:)]];
            %thisrun_FD = [sum(abs(thisrun_params_delta),2)];
            thisrun_params_delta_sq = thisrun_params_delta .^2;
            thisrun_FD = sum([sqrt(sum(thisrun_params_delta_sq(:,1:3),2)) sqrt(sum(thisrun_params_delta_sq(:,4:6),2))],2);
            
            FD = [FD ; thisrun_FD];
            
            delete([data(1:end-7) '_mcf*'])
            
        end
        
    end
    
    
    
    tmask = FD <= FDthresh;
    
    
    chunkedges = [1; diff([tmask;0])];
    chunkstarts = find(chunkedges==1); chunkstops = (find(chunkedges==-1) -1);
    chunksizes = zeros(length(chunkstarts),1);
    chunkID = zeros(size(tmask));
    for chunknum = 1:length(chunkstarts)
        chunkID(chunkstarts(chunknum) : chunkstops(chunknum)) = chunknum;
        chunksizes(chunknum) = nnz(chunkID==chunknum);
        if chunksizes(chunknum) < contiguoustimepoints
            tmask(chunkID==chunknum) = 0;
        end
    end
    
    
    pct_retained = nnz(tmask) / numel(tmask);
    disp(['Subject ' subject ', ' functional_sequences{seq} ': ' num2str(pct_retained*100) '% data retained from ' num2str(numel(tmask)) ' timepoints in session ' num2str(runs)])
    
    %dlmwrite([subject_rawdir functional_sequences{seq} '_FD.txt'],FD,'\t')
    %dlmwrite([subject_rawdir functional_sequences{seq} '_all_tmask.txt'],tmask,'\t')
    
    indices = 1:length(FD);
    figure; plot(indices,FD,'-r',indices,repmat(FDthresh,length(FD),1),'--k')
    %export_fig(gcf,[subject_rawdir functional_sequences{seq} '_FD.pdf'])
    
    
end



end

function failed = copyfiles_frommedia_firstsess(subject,thistargetdir)

media_root = '/media/';


breakout = 0;
medianames = dir([media_root '/*']);
for medianum = 3:length(medianames)
    if strcmp(medianames(medianum).name(1:length(subject)),subject) && (strcmp(medianames(medianum).name(end-1:end),'01') || strcmp(medianames(medianum).name(end-1:end),'-1')) && exist([media_root '/' medianames(medianum).name '/DICOM'],'dir')
        thissourcedir = [media_root '/' medianames(medianum).name '/'];
        break
    else
        subdirs = dir([media_root '/' medianames(medianum).name '/']);
        for subdirnum = 3:length(subdirs);
            if strcmp(medianames(medianum).name(1:length(subject)),subject) && (strcmp(medianames(medianum).name(end-1:end),'01') || strcmp(medianames(medianum).name(end-1:end),'-1')) && exist([media_root '/' medianames(medianum).name '/DICOM'],'dir')
                thissourcedir = [media_root '/' medianames(medianum).name '/' subdirs(subdirnum).name '/'];
                breakout = 1;
                break
            end
        end
        if breakout
            break
        end
    end
end


disp(['Copying from ' thissourcedir ' to ' thistargetdir ' and converting data'])

mkdir(thistargetdir);

[failed, message] = system(['cp -r ' thissourcedir '/* ' thistargetdir ]);
if failed
    disp('Copy FAILED!')
    disp(message)
end
[failed, message] = system(['dcm2nii ' thistargetdir ]);
if failed
    disp('Dicom conversion FAILED!')
    disp(message)
end


end