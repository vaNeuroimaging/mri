outfolder = '/home/data/Analysis/ROBI/';
subjectnames = {'ROBI003'};%{'Movie002','Movie003','Movie005','Movie006'};
bootstrap_iterations = 100;
for s = 1:length(subjectnames)
    
    subjectname = subjectnames{s};
    
    RSFC = ft_read_cifti_mod(['/home/data/subjects/' subjectname '/cifti/cifti_timeseries_normalwall/RSFC_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii']);
    RSFC_tc = RSFC.data; RSFC.data = []; outstruct = RSFC; clear RSFC
    [runs,sessions] = textread(['/home/data/subjects/' subjectname '/fc_processed/RSFC_runs_sessions.txt'],'%f%f');
    tmask = logical(load(['/home/data/subjects/' subjectname '/fc_processed/RSFC_all_tmask.txt']));
    runs = runs(tmask); sessions = sessions(tmask);
    
    Pre = single(paircorr_mod(RSFC_tc(:,sessions<100)'));
    Pre(isnan(Pre)) = 0;
    Pre = FisherTransform(Pre);
    
    Post = single(paircorr_mod(RSFC_tc(:,sessions>100)'));
    Post(isnan(Post)) = 0;
    Post = FisherTransform(Post);
    
    
    correlation_out = zeros(size(RSFC_tc,1),1);
    difference_out = zeros(size(RSFC_tc,1),1);
    for i = 1:size(RSFC_tc,1)
        correlation_out(i) = paircorr_mod(double(Pre(:,i)),double(Post(:,i)));
        difference_out(i) = mean(abs(Post(:,i) - Pre(:,i)));
    end
    
    clear Pre Post
    outstruct.data = correlation_out;
    outstruct.dimord = 'pos_time';
    ft_write_cifti_mod([outfolder '/' subjectname '_Pre_v_Post_correlation'],outstruct)
    
    outstruct.data = difference_out;
    ft_write_cifti_mod([outfolder '/' subjectname '_Pre_v_Post_difference'],outstruct)
    
    prevstring = [];
    
    
    runIDs = unique(runs);
    PrePostIDs_byrun = zeros(length(runIDs),1);
    for i = 1:length(runIDs)
        PrePostIDs_byrun(i) = single(mode(sessions(runs==runIDs(i))) > 100) +1;
    end
    
    
    correlation_bootstrap_out = zeros(size(RSFC_tc,1),bootstrap_iterations);
    difference_bootstrap_out = zeros(size(RSFC_tc,1),bootstrap_iterations);
    outstruct.data = [];
    
    
    for i = 1:bootstrap_iterations
        
        randPrePostIDs_byrun = PrePostIDs_byrun(randperm(length(PrePostIDs_byrun)));
        
        randPrePostIDs_bytime = zeros(size(runs));
        for j = 1:length(runIDs)
            randPrePostIDs_bytime(runs==runIDs(j)) = randPrePostIDs_byrun(j);
        end
        
        
        Pre = single(paircorr_mod(RSFC_tc(:,randPrePostIDs_bytime==1)'));
        Pre(isnan(Pre)) = 0;
        Pre = FisherTransform(Pre);
        
        Post = single(paircorr_mod(RSFC_tc(:,randPrePostIDs_bytime==2)'));
        Post(isnan(Post)) = 0;
        Post = FisherTransform(Post);
        
        for j = 1:length(Pre)
            string = [subjectname ': iteration ' num2str(i) ', vertex ' num2str(j)];
            fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
            prevstring = string;
            correlation_bootstrap_out(j,i) = paircorr_mod(double(Pre(:,j)),double(Post(:,j)));
            difference_bootstrap_out(j,i) = mean(abs(Post(:,i) - Pre(:,i)));
        end
        clear Pre Post
    end
    
    outstruct.data = correlation_bootstrap_out;
    ft_write_cifti_mod([outfolder '/' subjectname '_Pre_v_Post_correlation_bootstraps'],outstruct)
    
    bootstrap_percentile = sum(correlation_bootstrap_out > repmat(correlation_out,1,bootstrap_iterations),2);
    outstruct.data = bootstrap_percentile;
    ft_write_cifti_mod([outfolder '/' subjectname '_Pre_v_Post_correlation_bootstrap_percentile'],outstruct)
    
    outstruct.data = difference_bootstrap_out;
    ft_write_cifti_mod([outfolder '/' subjectname '_Pre_v_Post_difference_bootstraps'],outstruct)
    
    bootstrap_percentile = sum(difference_bootstrap_out > repmat(difference_out,1,bootstrap_iterations),2);
    outstruct.data = bootstrap_percentile;
    ft_write_cifti_mod([outfolder '/' subjectname '_Pre_v_Post_difference_bootstrap_percentile'],outstruct)
    
    disp(' ')
    
end