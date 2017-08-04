subjects = '/home/data/subjects/processing_list_062017.txt';
subjects = '/home/data/subjects/processing_list_080217.txt';
sessions_totest = [101 102 103];
network_pairs_totest = {[1 1],[3 3],[5 5],[10 10],[1 3],[1 5],[1 10]};
network_pairs_totest_names = {'DMN','FP','DA','SMh','DMNvFP','DMNvDA','DMNvSMh'};
outfolder = '/home/data/Analysis/ROBI/';
communities = load('/home/data/atlases/Group_parcellation/Parcel_Communities.txt');



if iscell(subjects)
elseif exist(subjects,'file')
    subjects = textread(subjects,'%s');
elseif ischar(subjects)
    subjects = {subjects};
end

for s = 1%:length(subjects)
    
    communities = ft_read_cifti_mod(['/home/data/subjects/' subjects{s} '/infomap/RSFC/rawassn_minsize400_regularized_recolored.dscalar.nii']);
    communities = communities.data;
    
    datafile = ['/home/data/subjects/' subjects{s} '/cifti/cifti_timeseries_normalwall/RSFC_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii'];
    data = ft_read_cifti_mod(datafile);
    
    runs_sessionsfile = ['/home/data/subjects/' subjects{s} '/fc_processed/RSFC_runs_sessions.txt'];
    
    runs_sessions = load(runs_sessionsfile);
    sessions = runs_sessions(:,2);
    clear run_sessions
    
    tmask = load(['/home/data/subjects/' subjects{s} '/fc_processed/RSFC_all_tmask.txt']);
    sessions = sessions(logical(tmask));
    
    
    corr_pre = paircorr_mod(data.data(:,sessions<100)');
    corr_pre(isnan(corr_pre)) = 0;
    corr_pre = FisherTransform(corr_pre);
    
    legendnames = {'pre'};
    
    network_pairs_fc = zeros(length(network_pairs_totest),length(sessions_totest)+1);
    
    for p = 1:length(network_pairs_totest)
        fcvals = corr_pre(communities==network_pairs_totest{p}(1),communities==network_pairs_totest{p}(2));
        network_pairs_fc(p,1) = mean(fcvals(:));
    end
    clear corr_pre
    
    possible_barcolors = {[0 1 0],[0 0 1],[1 1 0],[1 0 1],[0 1 1]};
    barcolors = {[1 0 0]};
    for sessnum = 1:length(sessions_totest)
        
        corr_post = paircorr_mod(data.data(:,sessions==sessions_totest(sessnum))');
        corr_post(isnan(corr_post)) = 0;
        corr_post = FisherTransform(corr_post);
        
        for p = 1:length(network_pairs_totest)
            fcvals = corr_post(communities==network_pairs_totest{p}(1),communities==network_pairs_totest{p}(2));
            network_pairs_fc(p,sessnum+1) = mean(fcvals(:));
        end
        
        legendnames(end+1) = {['Post' num2str(sessions_totest(sessnum))]};
        barcolors(end+1) = possible_barcolors(sessnum);
        clear corr_post
        
    end
    h = figure('Color','white','position',[ 949         474        1415         721],'DefaultAxesFontSize',20);
    b=bar(network_pairs_fc);
    for barnum = 1:length(barcolors)
        b(barnum).FaceColor = barcolors{barnum};
    end
    
    legend(legendnames)
    a=gca;
    a.XTickLabel = network_pairs_totest_names;
    title(subjects{s})
    export_fig(gca,[subjects{s} '_infomap_pairwise.pdf']);
    
end