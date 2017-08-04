function check_motion(subject,runs)

subdirroot = '/home/data/subjects/';

functional_sequences = {'restingstate'};

FDthresh = .2;
brainradius = 50;
contiguoustimepoints = 5;

if ~exist('runs')
    runs = 1:1000;
end


warning off

for seq = 1:length(functional_sequences)
    
    FD = [];
    
    %for runnum = runs
    runcounter = 0;
        
        subject_rawdir = [subdirroot subject '/preprocessed/'];
        
        if exist(subject_rawdir,'dir')
            
           for run = runs
            
            this_file = [subject_rawdir functional_sequences{seq} '_' num2str(run) '.nii.gz'];
            if exist(this_file)
            %disp(['Getting motion estimates for ' this_file])
            
            runcounter = runcounter+1;
                
                paramsfile = [this_file(1:end-7) '_mcf.par'];
                
                if ~exist(paramsfile,'file')
                    [fail, result] = system(['mcflirt -in ' this_file ' -refvol 0 -plots']);
                    if fail
                        error(result)
                    end
                end
                
                thisrun_params = load(paramsfile);
                thisrun_rot = thisrun_params(:,1:3);
                thisrun_rot_mm = thisrun_rot * brainradius;
                thisrun_params(:,1:3) = thisrun_rot_mm;
                thisrun_params_delta = [zeros(1,6) ; diff(thisrun_params)];
                thisrun_FD = [sum(abs(thisrun_params_delta),2)];
                %thisrun_params_delta_sq = thisrun_params_delta .^2;
                %thisrun_FD = sum([sqrt(sum(thisrun_params_delta_sq(:,1:3),2)) sqrt(sum(thisrun_params_delta_sq(:,4:6),2))],2);
                
                FD = [FD ; thisrun_FD];
                
                %delete([data(1:end-7) '_mcf.nii.gz'])
                
            else
                break
            end
           end
            
        end
    %end
    
    if runcounter > 0
    
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
    disp(['Subject ' subject ', ' functional_sequences{seq} ': ' num2str(pct_retained*100) '% data retained from ' num2str(numel(tmask)) ' timepoints from ' num2str(runcounter) ' runs'])
    
    %dlmwrite([subject_rawdir functional_sequences{seq} '_FD.txt'],FD,'\t')
    %dlmwrite([subject_rawdir functional_sequences{seq} '_all_tmask.txt'],tmask,'\t')
    
    plotted_tmask = single(tmask);
    plotted_tmask(tmask==0) = NaN;
    
    indices = 1:length(FD);
    %figure; plot(indices,FD,'-r',indices,repmat(FDthresh,length(FD),1),'--k',indices,(plotted_tmask-1),'.g')
    %export_fig(gcf,[subject_rawdir functional_sequences{seq} '_FD.pdf'])
    end
    
    
end



