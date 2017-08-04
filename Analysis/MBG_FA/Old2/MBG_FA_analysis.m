[subjects, urine] = textread('MBG_urine_levels.txt','%s%f');
for s = 1:length(subjects)
    
    [runs,sessions] = textread(['/home/data/subjects/' subjects{s} '/DTI/DTIruns_sessions.txt'],'%s%f');
    sess1inds = find(sessions==1);
    for sess = 1:length(sess1inds)
        sessdata = load_untouch_nii_2D(['/home/data/subjects/' subjects{s} '/DTI/DTI_' num2str(sess1inds(sess)) '_ec_FA_MNI.nii.gz']);
        subdata(:,sess) = sessdata.img;
    end
    
alldata(:,s) = mean(subdata,2);
clear subdata
end

tracts = load_untouch_nii_2D('/usr/local/fsl/data/atlases/JHU/JHU-ICBM-tracts-maxprob-thr25-2mm.nii.gz');

%% Urine

urine_subinds = ~isnan(urine);
urinecorr = corr(alldata(:,urine_subinds)',urine(urine_subinds),'type','Spearman');
urinecorr(isnan(urinecorr)) = 0;
sessdata.img = urinecorr;
save_untouch_nii_2D(sessdata,'Urine_MBG_vs_FA.nii.gz')

urinecorr = urinecorr .* logical(tracts.img);
sessdata.img = urinecorr;
save_untouch_nii_2D(sessdata,'Urine_MBG_vs_FA_intracts_thr25.nii.gz')

sessdata.img = -urinecorr;
save_untouch_nii_2D(sessdata,'Urine_MBG_vs_FA_intracts_thr25_neg.nii.gz')

system('cluster -i Urine_MBG_vs_FA_intracts_thr25_neg.nii.gz -t .4 --minextent=20 --othresh=Urine_MBG_vs_FA_intracts_thr25_neg_rthreshp4_clusterthresh20.nii.gz')
system('cluster -i Urine_MBG_vs_FA_intracts_thr25.nii.gz -t .4 --minextent=20 --othresh=Urine_MBG_vs_FA_intracts_thr25_rthreshp4_clusterthresh20.nii.gz')
system('fslmaths Urine_MBG_vs_FA_intracts_thr25_neg_rthreshp4_clusterthresh20.nii.gz -mul -1 -add Urine_MBG_vs_FA_intracts_thr25_rthreshp4_clusterthresh20.nii.gz Urine_MBG_vs_FA_intracts_thr25_posneg_rthreshp4_clusterthresh20.nii.gz')

tractIDs = unique(tracts.img); tractIDs(tractIDs<1) = [];
out_corr = sessdata;
out_ps = sessdata;
for ID = tractIDs(:)'
    FAs_intract = mean(alldata(tracts.img==ID,urine_subinds),1)';
    [FA_vMBG_corr,FA_vMBG_ps] = corr(FAs_intract,urine(urine_subinds),'type','Spearman');
    out_corr.img(tracts.img==ID) = FA_vMBG_corr;
    out_ps.img(tracts.img==ID) = FA_vMBG_ps;
end
save_untouch_nii_2D(out_corr,'Urine_MBG_vs_FA_tractavg_corr.nii.gz')
save_untouch_nii_2D(out_ps,'Urine_MBG_vs_FA_tractavg_ps.nii.gz')
out_ps.img = out_ps.img .* numel(tractIDs);
save_untouch_nii_2D(out_ps,'Urine_MBG_vs_FA_tractavg_ps_corrected.nii.gz')


% %% Serum
% 
% serum_subinds = ~isnan(serum);
% serumcorr = paircorr_mod(alldata(:,serum_subinds)',serum(serum_subinds));
% serumcorr(isnan(serumcorr)) = 0;
% subdata.img = serumcorr;
% save_untouch_nii_2D(subdata,'Serum_MBG_vs_FA.nii.gz')
% 
% serumcorr = serumcorr .* logical(tracts.img);
% subdata.img = serumcorr;
% save_untouch_nii_2D(subdata,'Serum_MBG_vs_FA_intracts_thr25.nii.gz')
% 
% subdata.img = -serumcorr;
% save_untouch_nii_2D(subdata,'Serum_MBG_vs_FA_intracts_thr25_neg.nii.gz')
% 
% system('cluster -i Serum_MBG_vs_FA_intracts_thr25.nii.gz -t .5 --minextent=20 --othresh=Serum_MBG_vs_FA_intracts_thr25_rthreshp5_clusterthresh20.nii.gz')
% system('cluster -i Serum_MBG_vs_FA_intracts_thr25_neg.nii.gz -t .5 --minextent=20 --othresh=Serum_MBG_vs_FA_intracts_thr25_neg_rthreshp5_clusterthresh20.nii.gz')
% system('fslmaths Serum_MBG_vs_FA_intracts_thr25_neg_rthreshp5_clusterthresh20.nii.gz -mul -1 -add Serum_MBG_vs_FA_intracts_thr25_rthreshp5_clusterthresh20.nii.gz Serum_MBG_vs_FA_intracts_thr25_posneg_rthreshp5_clusterthresh20.nii.gz')


