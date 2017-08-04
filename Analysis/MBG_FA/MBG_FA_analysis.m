subjects = {'MBG001','MBG002','MBG004','MBG005','MBG006','MBG007','MBG009','MBG010','MBG011','MBG012','MBG013','MBG014','MBG015','MBG016','MBG017','MBG018','MBG019','MBG020','MBG021','MBG022','MBG023','MBG024','MBG025','MBG026','MBG027','MBG028','MBG029','MBG030'};
[NUM,TXT,RAW]=xlsread('MBG_values_May2017.xlsx');

for s = 1:length(subjects)
    
    [runs,sessions] = textread(['/home/data/subjects/' subjects{s} '/DTI/DTIruns_sessions.txt'],'%s%f');
    sess1inds = find(sessions==1);
    for sess = 1:length(sess1inds)
        sessdata = load_untouch_nii_2D(['/home/data/subjects/' subjects{s} '/DTI/DTI_' num2str(sess1inds(sess)) '_ec_FA_MNI.nii.gz']);
        subdata(:,sess) = sessdata.img;
    end
    
    alldata(:,s) = mean(subdata,2);
    clear subdata

    subind_inexcel = strcmp(TXT,subjects{s});
    MBG(s,1) = NUM(subind_inexcel);
end

tracts = load_untouch_nii_2D('JHU-ICBM-labels-2mm_selected.nii.gz');

%% Serum

MBG_subinds = ~isnan(MBG);
MBGcorr = corr(alldata(:,MBG_subinds)',MBG(MBG_subinds),'type','Spearman');
MBGcorr(isnan(MBGcorr)) = 0;
sessdata.img = MBGcorr;
save_untouch_nii_2D(sessdata,'MBG_vs_FA.nii.gz')

MBGcorr = MBGcorr .* logical(tracts.img);
sessdata.img = MBGcorr;
save_untouch_nii_2D(sessdata,'MBG_vs_FA_intracts_thr25.nii.gz')

sessdata.img = -MBGcorr;
save_untouch_nii_2D(sessdata,'MBG_vs_FA_intracts_thr25_neg.nii.gz')

system('cluster -i MBG_vs_FA_intracts_thr25_neg.nii.gz -t .4 --minextent=20 --othresh=MBG_vs_FA_intracts_thr25_neg_rthreshp4_clusterthresh20.nii.gz --oindex=MBG_vs_FA_intracts_thr25_neg_rthreshp4_clusterthresh20_sep')
system('cluster -i MBG_vs_FA_intracts_thr25.nii.gz -t .4 --minextent=20 --othresh=MBG_vs_FA_intracts_thr25_rthreshp4_clusterthresh20.nii.gz --oindex=MBG_vs_FA_intracts_thr25_rthreshp4_clusterthresh20_sep')
system('fslmaths MBG_vs_FA_intracts_thr25_neg_rthreshp4_clusterthresh20.nii.gz -mul -1 -add MBG_vs_FA_intracts_thr25_rthreshp4_clusterthresh20.nii.gz MBG_vs_FA_intracts_thr25_posneg_rthreshp4_clusterthresh20.nii.gz')
%%
clusterinds_img = load_untouch_nii_2D('MBG_vs_FA_intracts_thr25_neg_rthreshp4_clusterthresh20_sep.nii.gz');
clusterinds = unique(clusterinds_img.img); clusterinds(clusterinds<1) = [];
for clusternum = 1:length(clusterinds)
    clustermeans = mean(alldata(clusterinds_img.img==clusterinds(clusternum),:),1);
    figure;
    plot(MBG(MBG_subinds),clustermeans(MBG_subinds),'k.','MarkerSize',20);
    hold on
    thisfit = fit(MBG(MBG_subinds),clustermeans(MBG_subinds)','poly1');
    plot(thisfit,'r-');
    xlim([min(MBG(MBG_subinds))-2 max(MBG(MBG_subinds))+2])
    ylim([max([0 min(clustermeans(MBG_subinds))-.01]),max(clustermeans(MBG_subinds))+.01])
    xlabel('MBG')
    ylabel('FA')
    legend off
    export_fig(gcf,['Cluster' num2str(clusternum) '_scatter.pdf'])
    close
end

prevclustcount = clusternum;
clusterinds_img = load_untouch_nii_2D('MBG_vs_FA_intracts_thr25_rthreshp4_clusterthresh20_sep.nii.gz');
clusterinds = unique(clusterinds_img.img); clusterinds(clusterinds<1) = [];
for clusternum = 1:length(clusterinds)
    clustermeans = mean(alldata(clusterinds_img.img==clusterinds(clusternum),:),1);
    figure;
    plot(MBG(MBG_subinds),clustermeans(MBG_subinds),'k.','MarkerSize',20);
    hold on
    thisfit = fit(MBG(MBG_subinds),clustermeans(MBG_subinds)','poly1');
    plot(thisfit,'r-');
    xlim([min(MBG(MBG_subinds))-2 max(MBG(MBG_subinds))+2])
    ylim([max([0 min(clustermeans(MBG_subinds))-.01]),max(clustermeans(MBG_subinds))+.01])
    xlabel('MBG')
    ylabel('FA')
    legend off
    export_fig(gcf,['Cluster' num2str(clusternum+prevclustcount) '_scatter.pdf'])
end

%%


tractIDs = unique(tracts.img); tractIDs(tractIDs<1) = [];
out_corr = sessdata;
out_ps = sessdata;
for ID = tractIDs(:)'
    FAs_intract = mean(alldata(tracts.img==ID,MBG_subinds),1)';
    [FA_vMBG_corr,FA_vMBG_ps] = corr(FAs_intract,MBG(MBG_subinds),'type','Spearman');
    out_corr.img(tracts.img==ID) = FA_vMBG_corr;
    out_ps.img(tracts.img==ID) = FA_vMBG_ps;
end
save_untouch_nii_2D(out_corr,'MBG_vs_FA_tractavg_corr.nii.gz')
save_untouch_nii_2D(out_ps,'MBG_vs_FA_tractavg_ps.nii.gz')
out_ps.img = out_ps.img .* numel(tractIDs);
save_untouch_nii_2D(out_ps,'MBG_vs_FA_tractavg_ps_corrected.nii.gz')


% %% Serum
% 
% MBG_subinds = ~isnan(MBG);
% MBGcorr = paircorr_mod(alldata(:,MBG_subinds)',MBG(MBG_subinds));
% MBGcorr(isnan(MBGcorr)) = 0;
% subdata.img = MBGcorr;
% save_untouch_nii_2D(subdata,'MBG_vs_FA.nii.gz')
% 
% MBGcorr = MBGcorr .* logical(tracts.img);
% subdata.img = MBGcorr;
% save_untouch_nii_2D(subdata,'MBG_vs_FA_intracts_thr25.nii.gz')
% 
% subdata.img = -MBGcorr;
% save_untouch_nii_2D(subdata,'MBG_vs_FA_intracts_thr25_neg.nii.gz')
% 
% system('cluster -i MBG_vs_FA_intracts_thr25.nii.gz -t .5 --minextent=20 --othresh=MBG_vs_FA_intracts_thr25_rthreshp5_clusterthresh20.nii.gz')
% system('cluster -i MBG_vs_FA_intracts_thr25_neg.nii.gz -t .5 --minextent=20 --othresh=MBG_vs_FA_intracts_thr25_neg_rthreshp5_clusterthresh20.nii.gz')
% system('fslmaths MBG_vs_FA_intracts_thr25_neg_rthreshp5_clusterthresh20.nii.gz -mul -1 -add MBG_vs_FA_intracts_thr25_rthreshp5_clusterthresh20.nii.gz MBG_vs_FA_intracts_thr25_posneg_rthreshp5_clusterthresh20.nii.gz')


