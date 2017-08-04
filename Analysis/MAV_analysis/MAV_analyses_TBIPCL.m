
load('/home/meilin/MAV/TBI_PCL_data.mat')

%% Comparing TBI yes/no PCL scores (t-test)

TBIyesall_indices = find(TBI_all==1);
TBInoall_indices = find(TBI_all==0);

TBIyesall_PCL = PCL_all(TBIyesall_indices);
TBInoall_PCL = PCL_all(TBInoall_indices);

[h, p] = ttest2(TBIyesall_PCL, TBInoall_PCL)

%% TBI num and PCL scores (r)

[r, p] = corr(PCL_all, TBI_numall)

%% TBI sev and PCL scores (ANOVA)

TBIsev_indices = find(TBI_sevall>0);
TBI_sevall_overzero = TBI_sevall(TBIsev_indices);
TBI_sevall_PCL = PCL_all(TBIsev_indices);

[p,tbl,stats] = anova1(TBI_sevall_PCL, TBI_sevall_overzero,'off')

%% Between/Within 

cd /home/data/Analysis/MAV_analysis

subject = {'MAV006','MAV008','MAV009','MAV014','MAV015','MAV016','MAV024','MAV029','MAV033','MAV035','MAV036','MAV040','MAV041','MAV042','MAV043','MAV044','MAV045','MAV047','MAV053'};

for s = 1:length(subject)
    alldata(:,:,s) = smartload(['/home/data/subjects/' subject{s} '/connectome/RSFC_Parcels_LR_corr.mat']);
end

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

%%

[within_r,within_p] = corr(allwithin_subs(:),PCL(:))
[pos_r,pos_p] = corr(allacrosspos_subs(:), PCL(:))
[neg_r,neg_p] = corr(allacrossneg_subs(:), PCL(:))

plot(allacrosspos_subs(:), PCL(:),'k.');


%% FA

tracts = load_untouch_nii_2D('/usr/local/fsl/data/atlases/JHU/JHU-ICBM-labels-2mm.nii.gz');
for s = 1:length(subject)
    data = load_untouch_nii_2D(['/home/data/subjects/' subject{s} '/DTI/DTI_avg_ec_FA_MNI.nii.gz']);
    subFA = data.img(logical(tracts.img));
    mean_subFA(s) = mean(subFA);
end

[PCL_r,PCL_p] = corr(mean_subFA(:),PCL(:));

TBIyes_indices = find(TBI==1);
TBIno_indices = find(TBI==0);

TBIyes_FA = mean_subFA(TBIyes_indices);
TBIno_FA = mean_subFA(TBIno_indices);

[h, p] = ttest2(TBIyes_FA, TBIno_FA)

%% TBI --> FA

[FAwithin_r,FAwithin_p] = corr(allwithin_subs(:),mean_subFA(:))
[FApos_r,FApos_p] = corr(allacrosspos_subs(:), mean_subFA(:))
[FAneg_r,FAneg_p] = corr(allacrossneg_subs(:), mean_subFA(:))

plot(allacrosspos_subs(:), mean_subFA(:),'k.')





