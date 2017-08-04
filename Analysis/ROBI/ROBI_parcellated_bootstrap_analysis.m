subjects = {'ROBI003','ROBI004'};
outfolder = '/home/data/Analysis/ROBI/';
numbootstraps = 100000;

for s = 1%:length(subjects)
    
    parcelsfile = ['/home/data/subjects/' subjects{s} '/parcellation/RSFC_parcels_edgethresh_0.5.dtseries.nii']; %'/home/data/atlases/Group_parcellation/Parcels_LR.dtseries.nii';
    
    parcelslabelfile = ['/home/data/subjects/' subjects{s} '/parcellation/RSFC_parcels_edgethresh_0.5.dlabel.nii']; %'/home/data/atlases/Group_parcellation/Parcels_LR.dlabel.nii';
    
    if exist([outfolder '/' subjects{s} '_pre.dconn.nii'])
        system(['wb_command -cifti-parcellate ' outfolder '/' subjects{s} '_pre.dconn.nii ' parcelslabelfile ' COLUMN ' outfolder '/' subjects{s} '_pre_subparcels.pdconn.nii'])
    end
    if exist([outfolder '/' subjects{s} '_post.dconn.nii'])
        system(['wb_command -cifti-parcellate ' outfolder '/' subjects{s} '_post.dconn.nii ' parcelslabelfile ' COLUMN ' outfolder '/' subjects{s} '_post_subparcels.pdconn.nii'])
    end
%     if exist([outfolder '/' subjects{s} '_postminpre.dconn.nii'])
%         system(['wb_command -cifti-parcellate ' outfolder '/' subjects{s} '_postminpre.dconn.nii ' parcelslabelfile ' COLUMN ' outfolder '/' subjects{s} '_postminpre_subparcels.pdconn.nii'])
%     end
    
    
    prevstring = [];
    string = [subjects{s} ': real data'];
    fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
    prevstring = string;
    
    datafile = ['/home/data/subjects/' subjects{s} '/cifti/cifti_timeseries_normalwall/RSFC_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii'];
    data = ft_read_cifti_mod(datafile);
    
    tmask = load(['/home/data/subjects/' subjects{s} '/fc_processed/RSFC_all_tmask.txt']);
    
    runs_sessionsfile = ['/home/data/subjects/' subjects{s} '/fc_processed/RSFC_runs_sessions.txt'];
    
    runs_sessions = load(runs_sessionsfile);
    runs_sessions = runs_sessions(logical(tmask),:);
    
    runs = runs_sessions(:,1);
    sessions = runs_sessions(:,2);
    runIDs = unique(runs);
    clear run_sessions
    
    
    parcels = ft_read_cifti_mod(parcelsfile);
    parcelIDs = unique(parcels.data); parcelIDs(parcelIDs<1) = [];
    out = parcels;
    
    parceltcs = zeros(size(data.data,2),length(parcelIDs));
    for IDnum = 1:length(parcelIDs)
        parceltcs(:,IDnum) = mean(data.data(parcels.data==parcelIDs(IDnum),:),1);
    end
    
    
    clear data
    
    
    parcelcorrmaps_pre = paircorr_mod(parceltcs((sessions < 100),:));
    parcelcorrmaps_pre(isnan(parcelcorrmaps_pre)) = 0;  parcelcorrmaps_pre = FisherTransform(parcelcorrmaps_pre);
    
    parcelcorrmaps_post = paircorr_mod(parceltcs((sessions > 100),:));
    parcelcorrmaps_post(isnan(parcelcorrmaps_post)) = 0;  parcelcorrmaps_post = FisherTransform(parcelcorrmaps_post);
    
    parcelcorrmaps_preminpost = parcelcorrmaps_pre - parcelcorrmaps_post;
    
    parcelcorrmaps_preminpost_avg = mean(abs(parcelcorrmaps_preminpost),2);
    
    for IDnum = 1:length(parcelIDs)
        out.data(parcels.data==parcelIDs(IDnum)) = parcelcorrmaps_preminpost_avg(IDnum);
    end
    ft_write_cifti_mod([outfolder '/' subjects{s} '_postminpre_bysubparcel'],out)
    
    
    
    
    %bootstrap
    parcelcorrmaps_preminpost_avg_permute = zeros(length(parcelcorrmaps_preminpost_avg),numbootstraps);
    for i = 1:numbootstraps
        
        string = [subjects{s} ': permuted data, iteration ' num2str(i)];
        fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
        prevstring = string;
        
        %randsessions = sessions(randperm(length(sessions)));
        done = 0;
        while done == 0
            runorder = randperm(length(runIDs));
            randsessions = [];
            for runnum = 1:length(runorder)
                randsessions = [randsessions ; sessions(runs==runIDs(runorder(runnum)))];
            end
            if ~all((randsessions<100) == (sessions < 100))
                done = 1;
            end
        end
        
        parcelcorrmaps_pre = paircorr_mod(parceltcs((randsessions < 100),:));
        parcelcorrmaps_pre(isnan(parcelcorrmaps_pre)) = 0;  parcelcorrmaps_pre = FisherTransform(parcelcorrmaps_pre);
        
        parcelcorrmaps_post = paircorr_mod(parceltcs((randsessions > 100),:));
        parcelcorrmaps_post(isnan(parcelcorrmaps_post)) = 0;  parcelcorrmaps_post = FisherTransform(parcelcorrmaps_post);
        
        parcelcorrmaps_preminpost = parcelcorrmaps_pre - parcelcorrmaps_post;
        
        parcelcorrmaps_preminpost_avg_permute(:,i) = mean(abs(parcelcorrmaps_preminpost),2);
    end
    
    parcelcorrmaps_preminpost_pctile_ofpermuted = sum((repmat(parcelcorrmaps_preminpost_avg,1,numbootstraps) > parcelcorrmaps_preminpost_avg_permute),2) ./ numbootstraps;
    for IDnum = 1:length(parcelIDs)
        out.data(parcels.data==parcelIDs(IDnum)) = parcelcorrmaps_preminpost_pctile_ofpermuted(IDnum);
    end
    ft_write_cifti_mod([outfolder '/' subjects{s} '_postminpre_bysubparcel_percentileofpermuted'],out)
    
    parcelcorrmaps_preminpost_pctile_ofpermuted_sig = (1 - parcelcorrmaps_preminpost_pctile_ofpermuted) * length(parcelIDs);
    for IDnum = 1:length(parcelIDs)
        out.data(parcels.data==parcelIDs(IDnum)) = parcelcorrmaps_preminpost_pctile_ofpermuted_sig(IDnum);
    end
    ft_write_cifti_mod([outfolder '/' subjects{s} '_postminpre_bysubparcel_correctedsignifianceVSpermuted'],out)
    
    
    disp(' ') 
    
    
    
end