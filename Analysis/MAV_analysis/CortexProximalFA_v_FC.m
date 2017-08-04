[subjects,TBIstatus,age] = textread('MAV_subs_TBIstatus_age.txt','%s%s%s');

parcels = ft_read_cifti_mod('/home/data/atlases/Group_parcellation/Parcels_LR.dtseries.nii');
parcelIDs = unique(parcels.data); parcelIDs(parcelIDs<1) = [];

for dist = 1%1:7

    corrmats = zeros(length(parcelIDs),length(parcelIDs), length(subjects));
    FAvals = zeros(length(parcelIDs), length(subjects));
    
for s = 1:length(subjects)
    
    corrmats(:,:,s) = smartload(['/home/data/subjects/' subjects{s} '/connectome/RSFC_Parcels_LR_corr.mat']);
     
    
    FA = ft_read_cifti_mod(['/home/data/subjects/' subjects{s} '/DTI/FA_near_cortex_LR.dscalar.nii']);
    
    for IDnum = 1:length(parcelIDs)
        FAvals(IDnum,s) = mean(FA.data(parcels.data==parcelIDs(IDnum),dist));
    end
end

rsqmat = zeros(length(parcelIDs),length(parcelIDs));

for IDnum = 1:length(parcelIDs)
    for IDnum2 = (IDnum + 1):length(parcelIDs)
        
        [~,~,~,~,STATS] = regress(squeeze(corrmats(IDnum,IDnum2,:)),[FAvals(IDnum,:)' FAvals(IDnum2,:)' ones(length(subjects),1)]);
        rsqmat(IDnum,IDnum2) = STATS(3);
    end
end

parcel_correlmat_figmaker((rsqmat + rsqmat') < .01,'/home/data/atlases/Group_parcellation/Parcel_Communities.txt',[-.4 .4],['cpFA dist' num2str(dist+1) ' vs FC'])

end