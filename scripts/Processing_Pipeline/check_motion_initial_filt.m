function check_motion_initial_filt(subject,runs)

subdirroot = '/home/data/subjects/';

functional_sequences = {'RSFC'};

brainradius = 50;
contiguoustimepoints = 5;
lowpassfilt = .08; %Hz
highpassfilt = .009; %Hz
minpctframes = .1;

lowpassmotionfilt = .025; %Hz

FDthresh = .04;

TR = 3.0;

if ~exist('runs')
    runs = 1;
end

if strcmp(runs,'all')
    runs = 1:5;
end

warning off

lopassmotioncutoff=lowpassmotionfilt/(0.5/TR); 
[butta_motion, buttb_motion]=butter(1,lopassmotioncutoff,'low');

for seq = 1:length(functional_sequences)
    
    FD = [];
    
    for runnum = runs
        
        subject_rawdir = [subdirroot subject '/raw/' subject '-' num2str(runnum) '/'];
        
        if exist(subject_rawdir,'dir')
            
            
            
            thisseq_files = dir([subject_rawdir '*' functional_sequences{seq} '*.nii.gz']);
            disp(['Found ' num2str(length(thisseq_files)) ' ' functional_sequences{seq} ' runs in ' subject_rawdir '. Getting motion estimates.'])
            
            for f = 1:length(thisseq_files);
                
                data = [subject_rawdir thisseq_files(f).name];
                
                paramsfile = [data(1:end-7) '_mcf.par'];
                
                if ~exist(paramsfile,'file')
                    [fail, result] = system(['mcflirt -in ' data ' -refvol 0 -plots']);
                    if fail
                        error(result)
                    end
                end
                delete([paramsfile(1:end-4) '.nii.gz'])
                
                thisrun_params = load(paramsfile);
                
                thisrun_params_unfilt = load(paramsfile);
                thisrun_rot = thisrun_params_unfilt(:,1:3);
                thisrun_rot_mm = thisrun_rot * brainradius;
                thisrun_params_unfilt(:,1:3) = thisrun_rot_mm;
                %figure;
                %plot(thisrun_params_unfilt)
                
                thisrun_params = filtfilt(butta_motion,buttb_motion,double(thisrun_params_unfilt));
                
                thisrun_params_delta = [zeros(1,6) ; diff(thisrun_params)];
                thisrun_FD = [sum(abs(thisrun_params_delta),2)];
                
%                 thisrun_rot = thisrun_params(:,1:3);
%                 thisrun_rot_mm = thisrun_rot * brainradius;
%                 thisrun_params(:,1:3) = thisrun_rot_mm;
%                 thisrun_params_delta = [zeros(1,6) ; diff(thisrun_params)];
%                 thisrun_FD = [sum(abs(thisrun_params_delta),2)];
                %thisrun_params_delta_sq = thisrun_params_delta .^2;
                %thisrun_FD = sum([sqrt(sum(thisrun_params_delta_sq(:,1:3),2)) sqrt(sum(thisrun_params_delta_sq(:,4:6),2))],2);
                
                FD = [FD ; thisrun_FD];
                
                %delete([data(1:end-7) '_mcf.nii.gz'])
                
            end
            
        end
    end
    
    
    
    tmask = (FD <= FDthresh);
    
    
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
    
    plotted_tmask = single(tmask);
    plotted_tmask(tmask==0) = NaN;
    
    indices = 1:length(FD);
    figure; plot(indices,FD,'-r',indices,repmat(FDthresh,length(FD),1),'--k',indices,(plotted_tmask-1),'.g')
    %export_fig(gcf,[subject_rawdir functional_sequences{seq} '_FD.pdf'])
    
    
end



