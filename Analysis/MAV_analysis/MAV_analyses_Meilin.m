cd /home/data/Analysis/MAV_analysis

subjects = {'MAV006','MAV008','MAV014','MAV015','MAV016','MAV024','MAV029','MAV033','MAV035','MAV036','MAV040','MAV041','MAV042','MAV043','MAV044','MAV045','MAV047','MAV053'};
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

%%

FAC1_t = FAC1b.';
FAC2_t = FAC2b.';
FAC3_t = FAC3b.';
FAC4_t = FAC4b.';

[within_r1,within_p1] = corr(allwithin_subs(:),FAC1_t(:));
[within_r2,within_p2] = corr(allwithin_subs(:),FAC2_t(:));
[within_r3,within_p3] = corr(allwithin_subs(:),FAC3_t(:));
[within_r4,within_p4] = corr(allwithin_subs(:),FAC4_t(:));

[pos_r1,pos_p1] = corr(allacrosspos_subs(:), FAC1_t(:));
[pos_r2,pos_p2] = corr(allacrosspos_subs(:), FAC2_t(:));
[pos_r3,pos_p3] = corr(allacrosspos_subs(:), FAC3_t(:));
[pos_r4,pos_p4] = corr(allacrosspos_subs(:), FAC4_t(:));

[neg_r1,neg_p1] = corr(allacrossneg_subs(:), FAC1_t(:));
[neg_r2,neg_p2] = corr(allacrossneg_subs(:), FAC2_t(:));
[neg_r3,neg_p3] = corr(allacrossneg_subs(:), FAC3_t(:));
[neg_r4,neg_p4] = corr(allacrossneg_subs(:), FAC4_t(:));

withincorrcoeff = [within_r1, within_r2, within_r3, within_r4; within_p1, within_p2, within_p3, within_p4]
poscorrceoff = [pos_r1, pos_r2, pos_r3, pos_r4; pos_p1, pos_p2, pos_p3, pos_p4]
negcorrcoeff = [neg_r1, neg_r2, neg_r3, neg_r4; neg_p1, neg_p2, neg_p3, neg_p4]

%% Controlling for TBI yes/no

TBI_t = TBI.';

[within_TBIr1,within_TBIp1] = partialcorr(allwithin_subs(:),FAC1_t(:), TBI_t(:));
[within_TBIr2,within_TBIp2] = partialcorr(allwithin_subs(:),FAC2_t(:), TBI_t(:));
[within_TBIr3,within_TBIp3] = partialcorr(allwithin_subs(:),FAC3_t(:), TBI_t(:));
[within_TBIr4,within_TBIp4] = partialcorr(allwithin_subs(:),FAC4_t(:), TBI_t(:));

[pos_TBIr1,pos_TBIp1] = partialcorr(allacrosspos_subs(:), FAC1_t(:), TBI_t(:));
[pos_TBIr2,pos_TBIp2] = partialcorr(allacrosspos_subs(:), FAC2_t(:), TBI_t(:));
[pos_TBIr3,pos_TBIp3] = partialcorr(allacrosspos_subs(:), FAC3_t(:), TBI_t(:));
[pos_TBIr4,pos_TBIp4] = partialcorr(allacrosspos_subs(:), FAC4_t(:), TBI_t(:));

[neg_TBIr1,neg_TBIp1] = partialcorr(allacrossneg_subs(:), FAC1_t(:), TBI_t(:));
[neg_TBIr2,neg_TBIp2] = partialcorr(allacrossneg_subs(:), FAC2_t(:), TBI_t(:));
[neg_TBIr3,neg_TBIp3] = partialcorr(allacrossneg_subs(:), FAC3_t(:), TBI_t(:));
[neg_TBIr4,neg_TBIp4] = partialcorr(allacrossneg_subs(:), FAC4_t(:), TBI_t(:));

TBIwithincorrcoeff = [within_TBIr1, within_TBIr2, within_TBIr3, within_TBIr4; within_TBIp1, within_TBIp2, within_TBIp3, within_TBIp4]
TBIposcorrcoeff = [pos_TBIr1, pos_TBIr2, pos_TBIr3, pos_TBIr4; pos_TBIp1, pos_TBIp2, pos_TBIp3, pos_TBIp4]
TBInegcorrcoeff = [neg_TBIr1, neg_TBIr2, neg_TBIr3, neg_TBIr4; neg_TBIp1, neg_TBIp2, neg_TBIp3, neg_TBIp4]

%% Controlling for severity of worst TBI

TBI_sev_t = TBI_sev.';

[within_TBIsevr1,within_TBIsevp1] = partialcorr(allwithin_subs(:),FAC1_t(:), TBI_sev_t(:));
[within_TBIsevr2,within_TBIsevp2] = partialcorr(allwithin_subs(:),FAC2_t(:), TBI_sev_t(:));
[within_TBIsevr3,within_TBIsevp3] = partialcorr(allwithin_subs(:),FAC3_t(:), TBI_sev_t(:));
[within_TBIsevr4,within_TBIsevp4] = partialcorr(allwithin_subs(:),FAC4_t(:), TBI_sev_t(:));

[pos_TBIsevr1,pos_TBIsevp1] = partialcorr(allacrosspos_subs(:), FAC1_t(:), TBI_sev_t(:));
[pos_TBIsevr2,pos_TBIsevp2] = partialcorr(allacrosspos_subs(:), FAC2_t(:), TBI_sev_t(:));
[pos_TBIsevr3,pos_TBIsevp3] = partialcorr(allacrosspos_subs(:), FAC3_t(:), TBI_sev_t(:));
[pos_TBIsevr4,pos_TBIsevp4] = partialcorr(allacrosspos_subs(:), FAC4_t(:), TBI_sev_t(:));

[neg_TBIsevr1,neg_TBIsevp1] = partialcorr(allacrossneg_subs(:), FAC1_t(:), TBI_sev_t(:));
[neg_TBIsevr2,neg_TBIsevp2] = partialcorr(allacrossneg_subs(:), FAC2_t(:), TBI_sev_t(:));
[neg_TBIsevr3,neg_TBIsevp3] = partialcorr(allacrossneg_subs(:), FAC3_t(:), TBI_sev_t(:));
[neg_TBIsevr4,neg_TBIsevp4] = partialcorr(allacrossneg_subs(:), FAC4_t(:), TBI_sev_t(:));

TBIsevwithincorrcoeff = [within_TBIsevr1, within_TBIsevr2, within_TBIsevr3, within_TBIsevr4; within_TBIsevp1, within_TBIsevp2, within_TBIsevp3, within_TBIsevp4]
TBIsevposcorrcoeff = [pos_TBIsevr1, pos_TBIsevr2, pos_TBIsevr3, pos_TBIsevr4; pos_TBIsevp1, pos_TBIsevp2, pos_TBIsevp3, pos_TBIsevp4]
TBIsevnegcorrcoeff = [neg_TBIsevr1, neg_TBIsevr2, neg_TBIsevr3, neg_TBIsevr4; neg_TBIsevp1, neg_TBIsevp2, neg_TBIsevp3, neg_TBIsevp4]
 
%% Controlling for number of TBIs

TBI_num_t = TBI_num.';

[within_TBInumr1,within_TBInump1] = partialcorr(allwithin_subs(:),FAC1_t(:), TBI_num_t(:));
[within_TBInumr2,within_TBInump2] = partialcorr(allwithin_subs(:),FAC2_t(:), TBI_num_t(:));
[within_TBInumr3,within_TBInump3] = partialcorr(allwithin_subs(:),FAC3_t(:), TBI_num_t(:));
[within_TBInumr4,within_TBInump4] = partialcorr(allwithin_subs(:),FAC4_t(:), TBI_num_t(:));

[pos_TBInumr1,pos_TBInump1] = partialcorr(allacrosspos_subs(:), FAC1_t(:), TBI_num_t(:));
[pos_TBInumr2,pos_TBInump2] = partialcorr(allacrosspos_subs(:), FAC2_t(:), TBI_num_t(:));
[pos_TBInumr3,pos_TBInump3] = partialcorr(allacrosspos_subs(:), FAC3_t(:), TBI_num_t(:));
[pos_TBInumr4,pos_TBInump4] = partialcorr(allacrosspos_subs(:), FAC4_t(:), TBI_num_t(:));

[neg_TBInumr1,neg_TBInump1] = partialcorr(allacrossneg_subs(:), FAC1_t(:), TBI_num_t(:));
[neg_TBInumr2,neg_TBInump2] = partialcorr(allacrossneg_subs(:), FAC2_t(:), TBI_num_t(:));
[neg_TBInumr3,neg_TBInump3] = partialcorr(allacrossneg_subs(:), FAC3_t(:), TBI_num_t(:));
[neg_TBInumr4,neg_TBInump4] = partialcorr(allacrossneg_subs(:), FAC4_t(:), TBI_num_t(:));

TBInumwithincorrcoeff = [within_TBInumr1, within_TBInumr2, within_TBInumr3, within_TBInumr4; within_TBInump1, within_TBInump2, within_TBInump3, within_TBInump4]
TBInumposcorrcoeff = [pos_TBInumr1, pos_TBInumr2, pos_TBInumr3, pos_TBInumr4; pos_TBInump1, pos_TBInump2, pos_TBInump3, pos_TBInump4]
TBInumnegcorrcoeff = [neg_TBInumr1, neg_TBInumr2, neg_TBInumr3, neg_TBInumr4; neg_TBInump1, neg_TBInump2, neg_TBInump3, neg_TBInump4]

 
