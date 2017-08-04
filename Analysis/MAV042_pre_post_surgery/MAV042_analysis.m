[runs,sessions] = textread('/home/data/subjects/MAV042/fc_processed/RSFC_runs_sessions.txt');
tmask = load('/home/data/subjects/MAV042/fc_processed/RSFC_all_tmask.txt');
presessioninds = (sessions <= 4); 
presessioninds_tmasked = presessioninds(logical(tmask));
postsessioninds = (sessions > 4);
postsessioninds_tmasked = postsessioninds(logical(tmask));

data = ft_read_cifti_mod('/home/data/subjects/MAV042/cifti/cifti_timeseries_normalwall/RSFC_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii');

datapre = data.data(:,presessioninds_tmasked);
datapost = data.data(:,postsessioninds_tmasked);
data.data = [];


parcels = ft_read_cifti_mod('/home/data/subjects/MAV042/parcellation/RSFC_parcels_edgethresh_0.5.dtseries.nii');
parcelIDs = unique(parcels.data); parcelIDs(parcelIDs<1) = [];
assignments = load('/home/data/subjects/MAV042/parcellation/parcel_infomap/rawassn_minsize5_regularized_recoloredv4.txt');
assignments(assignments > 17) = 18;

pretcs = zeros(size(datapre,2),length(parcelIDs));
posttcs = zeros(size(datapost,2),length(parcelIDs));

for IDnum = 1:length(parcelIDs)
    pretcs(:,IDnum) = mean(datapre(parcels.data==parcelIDs(IDnum),:),1);
    posttcs(:,IDnum) = mean(datapost(parcels.data==parcelIDs(IDnum),:),1);
end

premat = paircorr_mod(pretcs);
premat(isnan(premat)) = 0;
premat = FisherTransform(premat);


postmat = paircorr_mod(posttcs);
postmat(isnan(postmat)) = 0;
postmat = FisherTransform(postmat);


parcel_correlmat_figmaker_alt(premat,assignments,[-.6 .6],'Pre-Surgery')
parcel_correlmat_figmaker_alt(postmat,assignments,[-.6 .6],'Post-Surgery')
parcel_correlmat_figmaker_alt(postmat - premat,assignments,[-.25 .25],'Post-Surgery Changes')


