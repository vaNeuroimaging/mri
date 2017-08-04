subjects = {'MAV006','MAV008', 'MAV009','MAV014','MAV015','MAV016','MAV024','MAV029','MAV033','MAV035','MAV036','MAV040','MAV041','MAV042','MAV043','MAV044','MAV045','MAV047','MAV053'};
tracts = load_untouch_nii_2D('/usr/local/fsl/data/atlases/JHU/JHU-ICBM-labels-2mm.nii.gz');

for s = 1:length(subjects)
    data = load_untouch_nii_2D(['/home/data/subjects/' subjects{s} '/DTI/DTI_avg_ec_FA_MNI.nii.gz']);
    subFA = data.img(logical(tracts.img));
    mean_subFA(s) = mean(subFA);
end

%% Behav Components

FAC1_t = FAC1.';
FAC2_t = FAC2.';
FAC3_t = FAC3.';
FAC4_t = FAC4.';

[FAC1_r,FAC1_p] = corr(mean_subFA(:),FAC1_t(:));
[FAC2_r,FAC2_p] = corr(mean_subFA(:),FAC2_t(:));
[FAC3_r,FAC3_p] = corr(mean_subFA(:),FAC3_t(:));
[FAC4_r,FAC4_p] = corr(mean_subFA(:),FAC4_t(:));    

FA_corrcomponents = [FAC1_r, FAC2_r, FAC3_r, FAC4_r; FAC1_p, FAC2_p, FAC3_p, FAC4_p]
    
%% Controlling for TBI

TBI_t = TBI.';

[FAC1_TBI_r,FAC1_TBI_p] = partialcorr(mean_subFA(:),FAC1_t(:), TBI_t(:));
[FAC2_TBI_r,FAC2_TBI_p] = partialcorr(mean_subFA(:),FAC2_t(:), TBI_t(:));
[FAC3_TBI_r,FAC3_TBI_p] = partialcorr(mean_subFA(:),FAC3_t(:), TBI_t(:));
[FAC4_TBI_r,FAC4_TBI_p] = partialcorr(mean_subFA(:),FAC4_t(:), TBI_t(:));

FA_corrcomponents_TBI = [FAC1_TBI_r, FAC2_TBI_r, FAC3_TBI_r, FAC4_TBI_r; FAC1_TBI_p, FAC2_TBI_p, FAC3_TBI_p, FAC4_TBI_p]

%% Controlling for severity of worst TBI

TBI_sev_t = TBI_sev.';

[FAC1_sev_r,FAC1_sev_p] = partialcorr(mean_subFA(:),FAC1_t(:), TBI_sev_t(:));
[FAC2_sev_r,FAC2_sev_p] = partialcorr(mean_subFA(:),FAC2_t(:), TBI_sev_t(:));
[FAC3_sev_r,FAC3_sev_p] = partialcorr(mean_subFA(:),FAC3_t(:), TBI_sev_t(:));
[FAC4_sev_r,FAC4_sev_p] = partialcorr(mean_subFA(:),FAC4_t(:), TBI_sev_t(:));

FA_corrcomponents_TBIsev = [FAC1_sev_r, FAC2_sev_r, FAC3_sev_r, FAC4_sev_r; FAC1_sev_p, FAC2_sev_p, FAC3_sev_p, FAC4_sev_p]

%% Controlling for number of TBIs

TBI_num_t = TBI_num.';

[FAC1_num_r,FAC1_num_p] = partialcorr(mean_subFA(:),FAC1_t(:), TBI_num_t(:));
[FAC2_num_r,FAC2_num_p] = partialcorr(mean_subFA(:),FAC2_t(:), TBI_num_t(:));
[FAC3_num_r,FAC3_num_p] = partialcorr(mean_subFA(:),FAC3_t(:), TBI_num_t(:));
[FAC4_num_r,FAC4_num_p] = partialcorr(mean_subFA(:),FAC4_t(:), TBI_num_t(:));

FA_corrcomponents_TBInum = [FAC1_num_r, FAC2_num_r, FAC3_num_r, FAC4_num_r; FAC1_num_p, FAC2_num_p, FAC3_num_p, FAC4_num_p]

%% TBI t-test

TBIyes_indices = find(TBI_t==1);
TBIno_indices = find(TBI_t==0);

TBIyes_FA = mean_subFA(TBIyes_indices);
TBIno_FA = mean_subFA(TBIno_indices);

[h, p] = ttest2(TBIyes_FA, TBIno_FA)

%% TBI sev ANOVA

ANOVAp = anova1(mean_subFA, TBI_sev)


%% All subjects with TBI and FA data

allsubjects = {'MAV004','MAV006','MAV008','MAV009','MAV014','MAV015','MAV016','MAV019','MAV024','MAV029','MAV033','MAV035','MAV036','MAV040','MAV041','MAV042','MAV043','MAV044','MAV045','MAV047','MAV053'};
alltracts = load_untouch_nii_2D('/usr/local/fsl/data/atlases/JHU/JHU-ICBM-labels-2mm.nii.gz');

for s = 1:length(allsubjects)
    alldata = load_untouch_nii_2D(['/home/data/subjects/' allsubjects{s} '/DTI/DTI_avg_ec_FA_MNI.nii.gz']);
    allsubFA = alldata.img(logical(alltracts.img));
    allmean_subFA(s) = mean(allsubFA);
end

%% TBI t-test

allTBI_t = allTBI.';

allTBIyes_indices = find(allTBI_t==1);
allTBIno_indices = find(allTBI_t==0);

allTBIyes_FA = allmean_subFA(TBIyes_indices);
allTBIno_FA = allmean_subFA(TBIno_indices);

[ah, ap] = ttest2(allTBIyes_FA, allTBIno_FA)

%% TBI sev ANOVA

allANOVAp = anova1(allmean_subFA, allTBI_sev)


        
        




