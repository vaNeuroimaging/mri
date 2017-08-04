subjects = {'MBG001','MBG002','MBG004','MBG005','MBG006','MBG007','MBG009','MBG010','MBG011','MBG012','MBG013','MBG014','MBG015','MBG016','MBG017','MBG018','MBG019','MBG020','MBG021','MBG022','MBG023','MBG024','MBG025','MBG026','MBG027','MBG028','MBG029','MBG030'};
[NUM,TXT,RAW]=xlsread('../MBG_values_May2017.xlsx');

%[subjects, serum] = textread('/home/data/Analysis/MBG_FA/MBG_serum_levels.txt','%s%f');
for s = 1:length(subjects)
    
    [runs,sessions] = textread(['/home/data/subjects/' subjects{s} '/DTI/DTIruns_sessions.txt'],'%s%f');
    sess1inds = find(sessions==1);
    for sess = 1:length(sess1inds)
        sessdata = load_untouch_nii_2D(['/home/data/subjects/' subjects{s} '/DTI/DTI_' num2str(sess1inds(sess)) '_ec_FA_MNI.nii.gz']);
        subdata(:,sess) = sessdata.img;
    end
    sessdata.img = mean(subdata,2);
    clear subdata
    save_untouch_nii_2D(sessdata,[subjects{s} '_FA.nii.gz'])
    system(['fslmaths /home/data/subjects/' subjects{s} '/DTI/T1_biascorr_bet_MNI_222.nii.gz -thr 100 -bin -mul ' subjects{s} '_FA.nii.gz ' subjects{s} '_FA.nii.gz'])
    
    subind_inexcel = strcmp(TXT,subjects{s});
    MBGvalues(s,1) = NUM(subind_inexcel);
end

system('tbss_1_preproc *_FA.nii.gz')
system('tbss_2_reg -T')
system('tbss_3_postreg -S')
system('tbss_4_prestats 0.2')

%%

all_FA_skel = load_untouch_nii_2D('stats/all_FA_skeletonised.nii.gz');
skelinds = any(all_FA_skel.img,2);
MBGcorr = paircorr_mod(all_FA_skel.img(skelinds,:)',MBGvalues);
MBGcorr(isnan(MBGcorr)) = 0;
all_FA_skel.img = zeros(size(all_FA_skel.img,1),1);
all_FA_skel.img(skelinds) = MBGcorr;
all_FA_skel.hdr.dime.dim(5) = 1;
save_untouch_nii_2D(all_FA_skel,'MBG_vs_skeletonizedFA.nii.gz')


all_FA_skel.img(skelinds) = -MBGcorr;
save_untouch_nii_2D(all_FA_skel,'MBG_vs_skeletonizedFA_neg.nii.gz')

system('cluster -i MBG_vs_skeletonizedFA_neg.nii.gz -t .5 --minextent=20 --othresh=MBG_vs_skeletonizedFA_neg_rthreshp5_clusterthresh20.nii.gz')
system('cluster -i MBG_vs_skeletonizedFA.nii.gz -t .5 --minextent=20 --othresh=MBG_vs_skeletonizedFA_rthreshp5_clusterthresh20.nii.gz')
system('fslmaths MBG_vs_skeletonizedFA_neg_rthreshp5_clusterthresh20.nii.gz -mul -1 -add MBG_vs_skeletonizedFA_rthreshp5_clusterthresh20.nii.gz MBG_vs_skeletonizedFA_posneg_rthreshp5_clusterthresh20.nii.gz')
