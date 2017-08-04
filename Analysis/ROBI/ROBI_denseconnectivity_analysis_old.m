subjects = {'ROBI008','ROBI002'};
outfolder = '/home/data/Analysis/ROBI/';

for s = 1:length(subjects)
    
    datafile = ['/home/data/subjects/' subjects{s} '/cifti/cifti_timeseries_normalwall/RSFC_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii'];
    data = ft_read_cifti_mod(datafile);
    
    datastruct = data; datastruct.data = [];
    datastruct.dimord = 'pos_pos';
    
    runs_sessionsfile = ['/home/data/subjects/' subjects{s} '/fc_processed/RSFC_runs_sessions.txt'];
    
    runs_sessions = load(runs_sessionsfile);
    sessions = runs_sessions(:,2);
    clear run_sessions
    
    tmask = load(['/home/data/subjects/' subjects{s} '/fc_processed/RSFC_all_tmask.txt']);
    sessions = sessions(logical(tmask));
    
    
    corr_pre = paircorr_mod(data.data(:,sessions<100)');
    corr_pre(isnan(corr_pre)) = 0;
    corr_pre = FisherTransform(corr_pre);
    
    datastruct.data = corr_pre;
    ft_write_cifti_mod([outfolder '/' subjects{s} '_pre'],datastruct);
    datastruct.data = [];
    
    clear corr_pre
    
    
    
    corr_post = paircorr_mod(data.data(:,sessions>100)');
    corr_post(isnan(corr_post)) = 0;
    corr_post = FisherTransform(corr_post);
    
    clear data
    
    datastruct.data = corr_post;
    ft_write_cifti_mod([outfolder '/' subjects{s} '_post'],datastruct);
    datastruct.data = [];
    
    corr_pre = ft_read_cifti_mod([outfolder '/' subjects{s} '_pre.dconn.nii']);
    
    corr_postminpre = corr_post - corr_pre.data;
    clear corr_pre corr_post
    
    datastruct.data = corr_postminpre;
    ft_write_cifti_mod([outfolder '/' subjects{s} '_postminpre'],datastruct);
    datastruct.data = [];
    
    postminpre_mean = mean(abs(corr_postminpre),2);
    datastruct.data = postminpre_mean;
    datastruct.dimord = 'pos_time';
    ft_write_cifti_mod([outfolder '/' subjects{s} '_postminpre'],datastruct);
    
    clear corr_postminpre
    
end