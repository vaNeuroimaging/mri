subjects = {'MBG001','MBG002','MBG004','MBG005','MBG006','MBG007','MBG009','MBG010','MBG011','MBG012','MBG013','MBG014','MBG015','MBG016','MBG017','MBG018','MBG019','MBG020','MBG021','MBG022','MBG023','MBG024','MBG025','MBG026','MBG027','MBG028','MBG029','MBG030'};
[NUM,TXT,RAW]=xlsread('MBG_values_May2017.xlsx');

tracts = load_untouch_nii_2D('/home/data/Analysis/MBG_FA/JHU-ICBM-labels-2mm_selected.nii.gz');
tractIDs = unique(tracts.img); tractIDs(tractIDs<1) = [];
tractFAvals = zeros(length(subjects),length(tractIDs));

for s = 1:length(subjects)
    
    FAdata = load_untouch_nii_2D(['/home/data/subjects/' subjects{s} '/DTI/DTI_avg_ec_FA_MNI.nii.gz']);
    
    WMmask = load_untouch_nii_2D(['/home/data/subjects/' subjects{s} '/DTI/aparc+aseg_cerebralwm_mask_222.nii.gz']);
    
    for IDnum = 1:length(tractIDs)
        tractFAvals(s,IDnum) = mean(FAdata.img((tracts.img==tractIDs(IDnum)) & logical(WMmask.img)));
    end
    
    meanFA_intracts(s) = mean(FAdata.img(logical(WMmask.img) & logical(tracts.img)));
    
    allFA(:,s) = FAdata.img;
    
    subind_inexcel = strcmp(TXT,subjects{s});
    MBGvalues(s,1) = NUM(subind_inexcel);
    
end

%%

for t = 1:length(tractIDs)
    nonan_inds = ~isnan(tractFAvals(:,t));
    [tracts_corr(t) tracts_p(t)] = corr(tractFAvals(nonan_inds,t), MBGvalues(nonan_inds));
    
end

%%

MBGvalues_t = MBGvalues.';

[r1, p1] = corr(meanFA_intracts(:), MBGvalues_t(:));

%%

voxelwise = paircorr_mod(allFA',MBGvalues) .* single(logical(tracts.img));
FAdata.img = voxelwise;
save_untouch_nii_2D(FAdata,'FA_vs_MBG.nii.gz')


