outfolder = '/home/data/Analysis/ROBI/';
subjectnames = {'ROBI003'};%{'Movie002','Movie003','Movie005','Movie006'};
bootstrap_iterations = 1000;
for s = 1:length(subjectnames)
    
    subjectname = subjectnames{s};
    parcelassigns = load(['/home/data/subjects/' subjectname '/parcellation/parcel_infomap/rawassn_minsize5_regularized_manual.txt']);
    
    parcelsfile = ['/home/data/subjects/' subjectname '/parcellation/RSFC_parcels_edgethresh_0.5.dtseries.nii'];
    parcels = ft_read_cifti_mod(parcelsfile);
    parcelIDs = unique(parcels.data); parcelIDs(parcelIDs<1) = [];
    datafile = ['/home/data/subjects/' subjectname '/cifti/cifti_timeseries_normalwall/RSFC_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii'];
    data = ft_read_cifti_mod(datafile);
    RSFC_tc = zeros(length(parcelIDs),size(data.data,2));
    for i = 1:length(parcelIDs);
        RSFC_tc(i,:) = mean(data.data(parcels.data==parcelIDs(i),:),1);
    end
    clear data;
    outstruct = parcels;
    
    %RSFC = ft_read_cifti_mod(['/home/data/subjects/' subjectname '/cifti/parcel_timecourses/RSFC_parcels_edgethresh_0.5.ptseries.nii']);
    %RSFC_tc = RSFC.data; RSFC.data = []; outstruct = RSFC; clear RSFC
    [runs,sessions] = textread(['/home/data/subjects/' subjectname '/fc_processed/RSFC_runs_sessions.txt'],'%f%f');
    tmask = logical(load(['/home/data/subjects/' subjectname '/fc_processed/RSFC_all_tmask.txt']));
    runs = runs(tmask); sessions = sessions(tmask);
    
    Pre = single(paircorr_mod(RSFC_tc(:,sessions<100)'));
    Pre(isnan(Pre)) = 0;
    Pre = FisherTransform(Pre);
    parcel_correlmat_figmaker(Pre,parcelassigns,[-.5 .5],[subjectname '_Pre-Treatment'])
    export_fig(gca,[subjectname '_Pre-Treatment_parcel_correlmat.pdf'])
    
    Post = single(paircorr_mod(RSFC_tc(:,sessions>100)'));
    Post(isnan(Post)) = 0;
    Post = FisherTransform(Post);
    parcel_correlmat_figmaker(Pre,parcelassigns,[-.5 .5],[subjectname '_Post-Treatment'])
    export_fig(gca,[subjectname '_Post-Treatment_parcel_correlmat.pdf'])
    
    parcel_correlmat_figmaker(Post-Pre,parcelassigns,[-.3 .3],[subjectname '_Post_minus_Pre-Treatment'])
    export_fig(gca,[subjectname '_Post_minus_Pre-Treatment_parcel_correlmat.pdf'])
    
    
    out = zeros(size(RSFC_tc,1),1);
    for i = 1:size(RSFC_tc,1)
        out(i) = paircorr_mod(double(Pre(:,i)),double(Post(:,i)));
    end
    
    clear Pre Post
    outstruct.data = out;
    %outstruct.dimord = 'pos_time';
    ft_write_cifti_mod([outfolder '/' subjectname '_Pre_v_Post_parcelcorrelation'],outstruct)
    
    prevstring = [];
    
    
    runIDs = unique(runs);
    PrePostIDs_byrun = zeros(length(runIDs),1);
    for i = 1:length(runIDs)
        PrePostIDs_byrun(i) = single(mode(sessions(runs==runIDs(i))) > 100) +1;
    end
    
    
    bootstrap_out = zeros(size(RSFC_tc,1),bootstrap_iterations);
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
            bootstrap_out(j,i) = paircorr_mod(double(Pre(:,j)),double(Post(:,j)));
        end
        clear Movie RSFC
    end
    
    outstruct.data = bootstrap_out;
    ft_write_cifti_mod([outfolder '/' subjectname '_Pre_v_Post_parcelcorrelation_bootstraps'],outstruct)
    
    bootstrap_percentile = sum(bootstrap_out > repmat(out,1,bootstrap_iterations),2);
    outstruct.data = bootstrap_percentile;
    ft_write_cifti_mod([outfolder '/' subjectname '_Pre_v_Post_parcelcorrelation_bootstrap_percentile'],outstruct)
    
    disp(' ')
    
end