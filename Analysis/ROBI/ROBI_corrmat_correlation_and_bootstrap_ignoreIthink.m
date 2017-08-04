cd /home/data/Analysis/ROBI/
subjectnames = {'ROBI003'};%{'Movie002','Movie003','Movie005','Movie006'};
bootstrap_iterations = 500;
for s = 1:length(subjectnames)
    
    subjectname = subjectnames{s};
    
    datafile = ['/home/data/subjects/' subjectname '/cifti/cifti_timeseries_normalwall/RSFC_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii'];
    data = ft_read_cifti_mod(datafile);
    
    datastruct = data; datastruct.data = [];
    datastruct.dimord = 'pos_pos';
    
    runs_sessionsfile = ['/home/data/subjects/' subjectname '/fc_processed/RSFC_runs_sessions.txt'];
    
    runs_sessions = load(runs_sessionsfile);
    sessions = runs_sessions(:,2);
    clear run_sessions
    
    tmask = load(['/home/data/subjects/' subjectname '/fc_processed/RSFC_all_tmask.txt']);
    sessions = sessions(logical(tmask));
    
    
    
    
    
    RSFC = ft_read_cifti_mod(['/home/data/subjects/' subjectname '/connectome/RSFC_corr.dconn.nii']);
    Movie = ft_read_cifti_mod(['/home/data/subjects/' subjectname '/connectome/Movie_corr.dconn.nii']);
    
    out = zeros(length(RSFC.data),1);
    for i = 1:length(RSFC.data)
        out(i) = paircorr_mod(RSFC.data(:,i),Movie.data(:,i));
    end
    
    outstruct = RSFC;
    clear Movie RSFC
    outstruct.data = out;
    outstruct.dimord = 'pos_time';
    ft_write_cifti_mod([subjectname '_RSFC_v_Movie_correlation'],outstruct)
    
    prevstring = [];
    
    %bootstrap
    RSFC_tc = ft_read_cifti_mod(['/home/data/subjects/' subjectname '/cifti/cifti_timeseries_normalwall/RSFC_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii']);
    [RSFCruns,~] = textread(['/home/data/subjects/' subjectname '/fc_processed/RSFC_runs_sessions.txt'],'%f%f');
    tmask = load(['/home/data/subjects/' subjectname '/fc_processed/RSFC_all_tmask.txt']);
    RSFCruns = RSFCruns(logical(tmask));
    RSFCrunIDs = unique(RSFCruns);
    
    Movie_tc = ft_read_cifti_mod(['/home/data/subjects/' subjectname '/cifti/cifti_timeseries_normalwall/Movie_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii']);
    [Movieruns,~] = textread(['/home/data/subjects/' subjectname '/fc_processed/Movie_runs_sessions.txt'],'%f%f');
    tmask = load(['/home/data/subjects/' subjectname '/fc_processed/Movie_all_tmask.txt']);
    Movieruns = Movieruns(logical(tmask)) + max(RSFCruns);
    MovierunIDs = unique(Movieruns);
    
    bootstrap_out = zeros(size(RSFC_tc.data,1),bootstrap_iterations);
    outstruct = RSFC_tc; outstruct.data = [];
    
    both = [RSFC_tc.data Movie_tc.data];
    bothrunIDvec = [RSFCruns; Movieruns];
    MovieRestIDvec = [ones(size(RSFCrunIDs)); ones(size(MovierunIDs)).*2];
    runIDs = [RSFCrunIDs; MovierunIDs];
    
    %MovieRestIDvec = [ones(1,size(RSFC_tc.data,2)) ones(1,size(Movie_tc.data,2)).*2];
    clear RSFC_tc Movie_tc
    for i = 1:bootstrap_iterations
        
        %randIDvec = MovieRestIDvec(randperm(length(MovieRestIDvec)));
                
        randMovieRestorder = MovieRestIDvec(randperm(length(MovieRestIDvec)));
        
        randIDvec = zeros(size(bothrunIDvec));
        for j = 1:length(randMovieRestorder)
            randIDvec(bothrunIDvec==runIDs(j)) = randMovieRestorder(j);
        end
        
        RSFC = single(paircorr_mod(both(:,randIDvec==1)')); RSFC(isnan(RSFC)) = 0; RSFC = FisherTransform(RSFC);
        Movie = single(paircorr_mod(both(:,randIDvec==2)')); Movie(isnan(Movie)) = 0; Movie = FisherTransform(Movie);
        for j = 1:length(RSFC)
            string = [subjectname ': iteration ' num2str(i) ', vertex ' num2str(j)];
            fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
            prevstring = string;
            bootstrap_out(j,i) = paircorr_mod(double(RSFC(:,j)),double(Movie(:,j)));
        end
        clear Movie RSFC
    end
    
    outstruct.data = bootstrap_out;
    ft_write_cifti_mod([subjectname '_RSFC_v_Movie_correlation_bootstrapsv2'],outstruct)
    
    bootstrap_percentile = sum(bootstrap_out > repmat(out,1,bootstrap_iterations),2);
    outstruct.data = bootstrap_percentile;
    ft_write_cifti_mod([subjectname '_RSFC_v_Movie_correlation_bootstrappedv2_percentile'],outstruct)
    
    disp(' ')
    
end