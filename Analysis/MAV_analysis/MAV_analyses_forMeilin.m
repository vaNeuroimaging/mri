cd /home/data/Analysis/MAV_analysis

subjects = {'MAV006','MAV008'};
for s = 1:length(subjects)
    alldata(:,:,s) = smartload(['/home/data/subjects/' subjects{s} '/connectome/RSFC_Parcels_LR_corr.mat']);
end


%%

parcel_communities = load('/home/data/atlases/Group_parcellation/Parcel_Communities.txt');

within_inds = false(333);
across_inds_pos = false(333);
across_inds_neg = false(333);
IDs = unique(parcel_communities); IDs(IDs<1) = [];

mean_connect_bynetwork = diag(2*ones(100,1));

for IDnum = 1:length(IDs)
    ID = IDs(IDnum);
    within_inds(parcel_communities==ID,parcel_communities==ID) = true;
    for IDnum2 = (IDnum+1) : length(IDs)
        ID2 = IDs(IDnum2);
        thiscombo_inds = false(333);
        thiscombo_inds(parcel_communities==ID,parcel_communities==ID2) = true;
        temp = mean(alldata,3);
        if mean(temp(thiscombo_inds)) > 0
            mean_connect_bynetwork(ID,ID2) = 1;
            across_inds_pos(thiscombo_inds) = true;
        elseif mean(temp(thiscombo_inds)) < 0
            across_inds_neg(thiscombo_inds) = true;
            mean_connect_bynetwork(ID,ID2) = -1;
        end
    end
end
within_inds(diag(true(333,1),0)) = false;
within_inds(tril(true(333),-1)) = false;
across_inds = ~within_inds;
across_inds(diag(true(333,1),0)) = false;

indsmat = zeros(333);
indsmat(within_inds) = 1;
indsmat(across_inds_pos) = 5;
indsmat(across_inds_neg) = -5;

%%

for i = 1:size(alldata,3)
    temp = alldata(:,:,i);
    allwithin_subs(i) = mean(temp(within_inds));
    allacrosspos_subs(i) = mean(temp(across_inds_pos));
    allacrossneg_subs(i) = mean(temp(across_inds_neg));
end

