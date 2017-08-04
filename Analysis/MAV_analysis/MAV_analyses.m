cd /home/data/Analysis/MAV_analysis

initsubjects = textread('MAV_list.txt','%s');
initsubjects(strcmp(initsubjects,'MAV019')) = [];
%initsubjects(strcmp(initsubjects,'MAV066')) = [];
%initsubjects(strcmp(initsubjects,'MAV065')) = [];
initsubjects(strcmp(initsubjects,'MAV067')) = [];

%[subjects, tbivec, ages,BAI,BDI,AUDIT,PCL,PSIQ] = textread('MAV_data.txt','%s%n%n%n%n%n%n%n');

%[~,~,data]=xlsread('MAV_AssessmentData_cleaned.xlsx');
[~,~,data]=xlsread('MAV_AssessmentData_cleanedv2.xlsx');
data(strcmp(data,'')) = {NaN};
data(strcmp(data,'#DIV/0!')) = {NaN};

subjects = cell(0,1);
subs_withdata_indvec = false(size(data,1),1);
for s = 1:size(data,1)
    subs_withdata_indvec(s) = any(strcmp(data{s,1},initsubjects));
    if subs_withdata_indvec(s)
        subjects(end+1,1) = {data{s,1}};
    end
end


tbicol_inexcel = strcmp('Vasterling_Severity',data(1,:));
%tbivec = cell2mat(data(subs_withdata_indvec,tbicol_inexcel));

tbinumcol_inexcel = strcmp('Vasterling_Number',data(1,:));
%tbicountvec = cell2mat(data(subs_withdata_indvec,tbinumcol_inexcel));



PCLcol_inexcel = strcmp('PCL',data(1,:));
DSMBcol_inexcel = find(strcmp('DSM_B',data(1,:)));

FCEScol_inexcel = find(strcmp('FCES_34_Total',data(1,:)));

tbivec = zeros(length(subjects),1);
tbinumvec = zeros(length(subjects),1);
PCLvec = zeros(length(subjects),1);
FCESvec = zeros(length(subjects),1);
DSMBvec = zeros(length(subjects),1);
DSMCvec = zeros(length(subjects),1);
DSMDvec = zeros(length(subjects),1);
DSMEvec = zeros(length(subjects),1);

for s = 1:length(subjects)
    
    excel_subind = strcmp(subjects{s},data(:,1));
    
    if any(excel_subind)
        tbivec(s) = data{excel_subind,tbicol_inexcel};
        tbinumvec(s) = data{excel_subind,tbinumcol_inexcel};
        PCLvec(s) = data{excel_subind,PCLcol_inexcel};
        DSMBvec(s) = data{excel_subind,DSMBcol_inexcel};
        DSMCvec(s) = data{excel_subind,DSMBcol_inexcel+1};
        DSMDvec(s) = data{excel_subind,DSMBcol_inexcel+2};
        DSMEvec(s) = data{excel_subind,DSMBcol_inexcel+3};
        FCESvec(s) = data{excel_subind,FCEScol_inexcel};
        
    else
        tbivec(s) = NaN;
        tbinumvec(s) = NaN;
        PCLvec(s) = NaN;
        DSMBvec(s) = NaN;
        FCESvec(s) = NaN;
    end

    alldata(:,:,s) = smartload(['/home/data/subjects/' subjects{s} '/connectome/RSFC_Parcels_LR_corr.mat']);
end

tbi0 = alldata(:,:,tbivec==0);
tbi1 = alldata(:,:,tbivec==1);
tbi2 = alldata(:,:,tbivec==2);
tbi12 = alldata(:,:,tbivec>=1);

PCL_subinds = ~isnan(PCLvec);

%%

groupparcelfile = '/home/data/atlases/Group_parcellation/Parcels_LR.dlabel.nii';
groupparcelcommunitiesfile = '/home/data/atlases/Group_parcellation/Parcel_Communities.dtseries.nii';
controlavg = smartload('/home/data/atlases/Group_parcellation/Parcels_LR_avgcorr_120.mat');
parcels = ft_read_cifti_mod(groupparcelfile);
parcelIDs = unique(parcels.data); parcelIDs(parcelIDs<1) = [];
communities = ft_read_cifti_mod(groupparcelcommunitiesfile);
parcels.data = parcels.data(communities.brainstructure(1:min(length(communities.brainstructure),length(parcels.brainstructure))) > 0);
parcel_communities = zeros(length(parcelIDs),1);
for IDnum = 1:length(parcelIDs)
    parcel_communities(IDnum,1) = mode(communities.data(parcels.data==parcelIDs(IDnum)));
end



% parcel_correlmat_figmaker(mean(tbi0,3),parcel_communities,[-.5 .5],'No TBI');
% export_fig(gca,['No_TBI_corr.pdf'])
% parcel_correlmat_figmaker(mean(tbi1,3),parcel_communities,[-.5 .5],'Mild TBI');
% export_fig(gca,['Mild_TBI_corr.pdf'])
% parcel_correlmat_figmaker(mean(tbi2,3),parcel_communities,[-.5 .5],'Moderate TBI');
% export_fig(gca,['Moderate_TBI_corr.pdf'])
% parcel_correlmat_figmaker(mean(tbi0,3)-mean(tbi1,3),parcel_communities,[-.4 .4],'No TBI minus Mild TBI');
% export_fig(gca,['NoTBI_min_MildTBIcorr.pdf'])
% parcel_correlmat_figmaker(mean(tbi12,3),parcel_communities,[-.5 .5],'Any TBI');
% export_fig(gca,['AnyTBIcorr.pdf'])
% parcel_correlmat_figmaker(mean(tbi0,3) - mean(tbi12,3),parcel_communities,[-.4 .4],'No TBI minus Any TBI');
% export_fig(gca,['NoTBI_min_AnyTBIcorr.pdf'])


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
        if mean(controlavg(thiscombo_inds)) > 0
            mean_connect_bynetwork(ID,ID2) = 1;
            across_inds_pos(thiscombo_inds) = true;
        elseif mean(controlavg(thiscombo_inds)) < 0
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
indsmat(within_inds) = 3;
indsmat(across_inds_pos) = 2;
indsmat(across_inds_neg) = -2;

%parcel_correlmat_figmaker_v2(indsmat,parcel_communities,[-3 3])



for i = 1:size(alldata,3)
    temp = alldata(:,:,i);
    allwithin_subs(i) = mean(temp(within_inds));
    allacrosspos_subs(i) = mean(temp(across_inds_pos));
    allacrossneg_subs(i) = mean(temp(across_inds_neg));
end


%%  FA
%tracts = load_untouch_nii_2D('/usr/local/fsl/data/atlases/JHU/JHU-ICBM-labels-2mm.nii.gz');
tracts = load_untouch_nii_2D('JHU-ICBM-labels-2mm_selected.nii.gz');
%tracts = load_untouch_nii_2D('JHU-ICBM-labels-2mm-ero.nii.gz');
tractIDs = unique(tracts.img); tractIDs(tractIDs<1) = [];
tractFAvals = zeros(length(subjects),length(tractIDs));

for s = 1:length(subjects)
    
    
    FAdata = load_untouch_nii_2D(['/home/data/subjects/' subjects{s} '/DTI/DTI_avg_ec_FA_MNI.nii.gz']);
    
    WMmask = load_untouch_nii_2D(['/home/data/subjects/' subjects{s} '/DTI/aparc+aseg_cerebralwm_mask_222.nii.gz']);
    
    for IDnum = 1:length(tractIDs)
        tractFAvals(s,IDnum) = mean(FAdata.img((tracts.img==tractIDs(IDnum)) & logical(WMmask.img)));
    end
    
    meanFA_intracts(s) = mean(FAdata.img(logical(WMmask.img) & logical(tracts.img)));
    
    FAdata.img(~logical(WMmask.img)) = NaN;
    
    %allFA(:,s) = FAdata.img;
    
    
end

% %% modularity
% addpath /home/data/scripts/Resources/BrainConnectivityToolbox/
% modularities = zeros(size(alldata,3),1);
% for i = 1:size(alldata,3)
%     temp = alldata(:,:,i);
%     %[matrix,~,~] = matrix_thresholder_Evan(temp,0.05,'kden');
%     %[Q,~,~,~] = calc_modularity_TL(parcel_communities,logical(matrix));
%     [~,Q] = modularity_louvain_und_sign(temp);
%     modularities(i) = Q;
% end
    
    %% TBI severity vs PCL

[P,ANOVATAB,STATS] = anova1(PCLvec(PCL_subinds)',tbivec(PCL_subinds));
ANOVATAB

figure('Color','white','position',[1982 478 1352 804]);
bh = boxplot(PCLvec(PCL_subinds),tbivec(PCL_subinds),'colors','k','symbol','k','labels',{'No TBI','Mild TBI','Moderate TBI'});
set(bh,'linewidth',3)
set(bh,'MarkerEdgeColor',[0 0 0])
set(gca,'FontSize',40,'FontWeight','bold','LineWidth',2)
set(gca,'YTick',[0:20:80])
box off
export_fig(gca,['PCL_byTBI.pdf'])

[H,P,CI,STATS] = ttest2(PCLvec(tbivec(PCL_subinds)==1)' ,PCLvec(tbivec(PCL_subinds)==0)')
[H,P,CI,STATS] = ttest2(PCLvec(tbivec(PCL_subinds)==2)' , PCLvec(tbivec(PCL_subinds)==1)')
[H,P,CI,STATS] = ttest2(PCLvec(tbivec(PCL_subinds)==2)' , PCLvec(tbivec(PCL_subinds)==0)')

% 
% figure;
% plotSpread([allwithin_subs(tbivec==0)';allwithin_subs(tbivec==1)';allwithin_subs(tbivec==2)'],'categoryIdx',[zeros(nnz(tbivec==0),1);ones(nnz(tbivec==1),1);ones(nnz(tbivec==2),1)*2],'categoryColors',{'b','r','g'},'MarkerSize',40)
% set(gcf,'Color',[1 1 1])
% set(gca,'FontSize',15)
% export_fig(gca,['MeanWithinNetworkConnectivity.pdf'])

%% TBI number vs PCL

test_subinds = (tbivec > 0) & (~isnan(PCLvec));%true(length(tbivec),1);%


[r,p] = corr(tbinumvec(test_subinds),PCLvec(test_subinds))%,'type','Spearman')
prettyscatterplot(tbinumvec(test_subinds),PCLvec(test_subinds),1)
export_fig(gca,'TBInum_vs_PCL.pdf')
%figure;plot(tbinumvec(test_subinds),PCLvec(test_subinds),'k.','MarkerSize',20)
%xlim([0 max(tbinumvec(test_subinds))+1])
[r,p] = partialcorr([PCLvec(PCL_subinds) tbivec(PCL_subinds) tbinumvec(PCL_subinds)])

%% TBI severity vs RSFC

[P,ANOVATAB,STATS] = anova1(allwithin_subs',tbivec);
ANOVATAB
 [H,P,CI,STATS] = ttest2([allwithin_subs(tbivec==1)' ; allwithin_subs(tbivec==2)'],allwithin_subs(tbivec==0)');
 STATS
 P
[H,P,CI,STATS] = ttest2(allwithin_subs(tbivec==1)',allwithin_subs(tbivec==2)');
% figure;
% plotSpread([allwithin_subs(tbivec==0)';allwithin_subs(tbivec==1)';allwithin_subs(tbivec==2)'],'categoryIdx',[zeros(nnz(tbivec==0),1);ones(nnz(tbivec==1),1);ones(nnz(tbivec==2),1)*2],'categoryColors',{'b','r','g'},'MarkerSize',40)
% set(gcf,'Color',[1 1 1])
% set(gca,'FontSize',15)
% export_fig(gca,['MeanWithinNetworkConnectivity.pdf'])

figure('Color','white','position',[1982 478 1352 804]);
bh = boxplot(allwithin_subs,tbivec,'colors','k','symbol','k','labels',{'No TBI','Mild TBI','Moderate TBI'});
set(bh,'linewidth',3)
set(bh,'MarkerEdgeColor',[0 0 0])
set(gca,'FontSize',40,'FontWeight','bold','LineWidth',2)
box off
 export_fig(gca,['MeanWithinNetworkConnectivity_byTBI.pdf'])
% 
% 
% 
% 
% [P,ANOVATAB,STATS] = anova1(allacrosspos_subs',tbivec);
% ANOVATAB
%  [H,P,CI,STATS] = ttest2([allacrosspos_subs(tbivec==1)' ; allacrosspos_subs(tbivec==2)'],allacrosspos_subs(tbivec==0)');
%  STATS
%  P
% 
% % figure;
% % plotSpread([allacrosspos_subs(tbivec==0)';allacrosspos_subs(tbivec==1)';allacrosspos_subs(tbivec==2)'],'categoryIdx',[zeros(nnz(tbivec==0),1);ones(nnz(tbivec==1),1);ones(nnz(tbivec==2),1)*2],'categoryColors',{'b','r','g'},'MarkerSize',40)
% % set(gcf,'Color',[1 1 1])
% % set(gca,'FontSize',15)
% % export_fig(gca,['MeanPositiveAcrossNetworkConnectivity.pdf'])
% 
% figure('Color','white','position',[1982 478 1352 804]);
% bh = boxplot(allacrosspos_subs,tbivec,'colors','k','symbol','k','labels',{'No TBI','Mild TBI','Moderate TBI'});
% set(bh,'linewidth',3)
% set(bh,'MarkerEdgeColor',[0 0 0])
% set(gca,'FontSize',20,'FontWeight','bold','LineWidth',2)
% export_fig(gca,['MeanPositiveAcrossNetworkConnectivity_byTBI.pdf'])
% 
% 
% 
% 
% [P,ANOVATAB,STATS] = anova1(allacrossneg_subs',tbivec);
% ANOVATAB
%  [H,P,CI,STATS] = ttest2([allacrossneg_subs(tbivec==1)' ; allacrossneg_subs(tbivec==2)'],allacrossneg_subs(tbivec==0)');
%  STATS
%  P
% 
% % figure;
% % plotSpread([allacrossneg_subs(tbivec==0)';allacrossneg_subs(tbivec==1)';allacrossneg_subs(tbivec==2)'],'categoryIdx',[zeros(nnz(tbivec==0),1);ones(nnz(tbivec==1),1);ones(nnz(tbivec==2),1)*2],'categoryColors',{'b','r','g'},'MarkerSize',40)
% % set(gcf,'Color',[1 1 1])
% % set(gca,'FontSize',15)
% % export_fig(gca,['MeanNegativeAcrossNetworkConnectivity.pdf'])
% 
% figure('Color','white','position',[1982 478 1352 804]);
% bh = boxplot(allacrossneg_subs,tbivec,'colors','k','symbol','k','labels',{'No TBI','Mild TBI','Moderate TBI'});
% set(bh,'linewidth',3)
% set(bh,'MarkerEdgeColor',[0 0 0])
% set(gca,'FontSize',20,'FontWeight','bold','LineWidth',2)
% export_fig(gca,['MeanNegativeAcrossNetworkConnectivity_byTBI.pdf'])


%% TBI number vs RSFC

test_subinds = tbivec > 0;%true(size(tbivec));%

[r_within,p_within] = corr(tbinumvec(test_subinds),allwithin_subs(test_subinds)')%,'type','Spearman')
prettyscatterplot(tbinumvec(test_subinds),allwithin_subs(test_subinds),1)
export_fig(gca,'TBInum_vs_WithinRSFC.pdf')

% figure;plot(tbinumvec(test_subinds),allwithin_subs(test_subinds),'k.','MarkerSize',20)
% xlim([0 max(tbinumvec(test_subinds))+1])
% 
% [r_acrosspos,p_acrosspos] = corr(tbinumvec(test_subinds),allacrosspos_subs(test_subinds)')%,'type','Spearman')
% figure;plot(tbinumvec(test_subinds),allacrosspos_subs(test_subinds),'k.','MarkerSize',20)
% xlim([0 max(tbinumvec(test_subinds))+1])
% 
% [r_acrossneg,p_acrossneg] = corr(tbinumvec(test_subinds),allacrossneg_subs(test_subinds)')%,'type','Spearman')
% figure;plot(tbinumvec(test_subinds),allacrossneg_subs(test_subinds),'k.','MarkerSize',20)
% xlim([0 max(tbinumvec(test_subinds))+1])


[r,p] = partialcorr(tbinumvec,allwithin_subs',tbivec)
% [r,p] = partialcorr(tbinumvec,allacrosspos_subs',tbivec)
% [r,p] = partialcorr(tbinumvec,allacrossneg_subs',tbivec)


%% RSFC vs PCL



[r_within,p_within] = corr(PCLvec(PCL_subinds),allwithin_subs(PCL_subinds)')%,'type','Spearman')
prettyscatterplot(allwithin_subs(PCL_subinds),PCLvec(PCL_subinds),1)
 export_fig(gca,'WithinFC_vs_PCL.pdf')
% 
% [r_acrosspos,p_acrosspos] = corr(PCLvec(PCL_subinds),allacrosspos_subs(PCL_subinds)')%,'type','Spearman')
% figure('Color','white','position',[1982 478 1352 804]);
% plot(allacrosspos_subs(PCL_subinds),PCLvec(PCL_subinds),'k.','MarkerSize',40)
% set(gca,'FontSize',20,'FontWeight','bold','LineWidth',2)
% hold on
% limits = [.02 .065];
% xlim(limits)
% ylim([0 max(PCLvec(PCL_subinds))+1])
% FO = fit(allacrosspos_subs(PCL_subinds)', PCLvec(PCL_subinds), 'poly1');
% fitted = [limits(1) allacrosspos_subs(PCL_subinds) limits(2)] .* FO.p1 + FO.p2;
% [sorted,sorti] = sort([limits(1) allacrosspos_subs(PCL_subinds) limits(2)]);
% plot(sorted,fitted(sorti),'r-','LineWidth',3)
% export_fig(gca,'AcrossPosFC_vs_PCL.pdf')
% % figure;plot(PCLvec(PCL_subinds),allacrosspos_subs(PCL_subinds),'k.','MarkerSize',20)
% % xlim([0 max(PCLvec(PCL_subinds))+1])
% 
% [r_acrossneg,p_acrossneg] = corr(PCLvec(PCL_subinds),allacrossneg_subs(PCL_subinds)')%,'type','Spearman')
% figure('Color','white','position',[1982 478 1352 804]);
% plot(allacrossneg_subs(PCL_subinds),PCLvec(PCL_subinds),'k.','MarkerSize',40)
% set(gca,'FontSize',20,'FontWeight','bold','LineWidth',2)
% hold on
% limits = [-.09 -.03];
% xlim(limits)
% ylim([0 max(PCLvec(PCL_subinds))+1])
% FO = fit(allacrossneg_subs(PCL_subinds)', PCLvec(PCL_subinds), 'poly1');
% fitted = [limits(1) allacrossneg_subs(PCL_subinds) limits(2)] .* FO.p1 + FO.p2;
% [sorted,sorti] = sort([limits(1) allacrossneg_subs(PCL_subinds) limits(2)]);
% plot(sorted,fitted(sorti),'r-','LineWidth',3)
% export_fig(gca,'AcrossNegFC_vs_PCL.pdf')
% % figure;plot(PCLvec(PCL_subinds),allacrossneg_subs(PCL_subinds),'k.','MarkerSize',20)
% % xlim([0 max(PCLvec(PCL_subinds))+1])




%% TBI severity vs PCL control for RSFC

[~,T,~,~]=anovan(PCLvec(PCL_subinds)',{tbivec(PCL_subinds),allwithin_subs(PCL_subinds)'},'continuous',2)

%[~,T,~,~]=anovan(PCLvec(PCL_subinds)',{tbivec(PCL_subinds),allacrosspos_subs(PCL_subinds)'},'continuous',2)

%[~,T,~,~]=anovan(PCLvec(PCL_subinds)',{tbivec(PCL_subinds),allacrossneg_subs(PCL_subinds)'},'continuous',2)



%% All measures vs RSFC



% for i = 2:size(data,2)
%     subsdata = cell2mat(data(subs_withdata_indvec,i));
%     hasdata = ~isnan(subsdata);
%     [rs_within(i-1),ps_within(i-1)] = corr(allwithin_subs(hasdata)',subsdata(hasdata),'type','Spearman');
%     [rs_acrosspos(i-1),ps_acrosspos(i-1)] = corr(allacrosspos_subs(hasdata)',subsdata(hasdata),'type','Spearman');
%     [rs_acrossneg(i-1),ps_acrossneg(i-1)] = corr(allacrossneg_subs(hasdata)',subsdata(hasdata),'type','Spearman');
% end




 %% Surf FA and thickness
% 
% for s = 1:length(subjects)
%     
%     FAsurf = ft_read_cifti_mod(['/home/data/subjects/' subjects{s} '/DTI/FA_near_cortex_LR.dscalar.nii']);
%     
%     Ldepth = gifti(['/home/data/subjects/' subjects{s} '/fs_LR/MNI/fsaverage_LR32k/' subjects{s} '.L.sulc.32k_fs_LR.shape.gii']);
%     Rdepth = gifti(['/home/data/subjects/' subjects{s} '/fs_LR/MNI/fsaverage_LR32k/' subjects{s} '.R.sulc.32k_fs_LR.shape.gii']);
%     
%     LRdepth = [Ldepth.cdata(FAsurf.brainstructure(1:32492)>0) ; Rdepth.cdata(FAsurf.brainstructure(32493:64984)>0)];
%     
%     meanFA_onsurf(s) = mean(FAsurf.data(LRdepth>0,2));
%     
%     allFA_onsurf(:,s) = FAsurf.data(:,2);
% %     
% %     
% %     
%     Lthick = gifti(['/home/data/subjects/' subjects{s} '/fs_LR/MNI/fsaverage_LR32k/' subjects{s} '.L.thickness.32k_fs_LR.shape.gii']);
%     Rthick = gifti(['/home/data/subjects/' subjects{s} '/fs_LR/MNI/fsaverage_LR32k/' subjects{s} '.R.thickness.32k_fs_LR.shape.gii']);
%     
%     LRthickness = [Lthick.cdata(FAsurf.brainstructure(1:32492)>0) ; Rthick.cdata(FAsurf.brainstructure(32493:64984)>0)];
%     
%     meancorticalthickness(s) = mean(LRthickness);
%     
%     allcorticalthickness(:,s) = LRthickness;
%     
% end

%% TBI severity vs FA

[P,ANOVATAB,STATS] = anova1(meanFA_intracts',tbivec);
ANOVATAB
figure('Color','white','position',[1982 478 1352 804]);
bh = boxplot(meanFA_intracts,tbivec,'colors','k','symbol','k','labels',{'No TBI','Mild TBI','Moderate TBI'});
set(bh,'linewidth',3)
set(bh,'MarkerEdgeColor',[0 0 0])
set(gca,'FontSize',40,'FontWeight','bold','LineWidth',2)
ylim([.46 .580001])
box off
export_fig(gca,['FA_byTBI.pdf'])

figure;
plotSpread([meanFA_intracts(tbivec==0)';meanFA_intracts(tbivec==1)';meanFA_intracts(tbivec==2)'],'categoryIdx',[zeros(nnz(tbivec==0),1);ones(nnz(tbivec==1),1);ones(nnz(tbivec==2),1)*2],'categoryColors',{'b','r','g'},'MarkerSize',40)
set(gcf,'Color',[1 1 1])
set(gca,'FontSize',15)


%% TBI number vs FA

%  [r_within,p_within] = corr(tbinumvec,meanFA_intracts')%,'type','Spearman')
%  prettyscatterplot(tbinumvec,meanFA_intracts,1)
%  export_fig(gca,'TBInum_vs_FA.pdf')
% figure;plot(tbinumvec,meanFA_intracts,'k.','MarkerSize',20)
% xlim([0 max(tbinumvec)+1])


test_subinds = tbivec > 0;

[r_within,p_within] = corr(tbinumvec(test_subinds),meanFA_intracts(test_subinds)')%,'type','Spearman')
prettyscatterplot(tbinumvec(test_subinds),meanFA_intracts(test_subinds),1)
 export_fig(gca,'TBInum_vs_FA.pdf')






%% TBI severity vs FA in each tract

pvals = zeros(size(tractFAvals,2),1);
Fs = zeros(size(tractFAvals,2),1);
for i = 1:size(tractFAvals,2)
    nonaninds = ~isnan(tractFAvals(:,i));
    [P,ANOVATAB,STATS] = anova1(tractFAvals(nonaninds,i),tbivec(nonaninds),'off');
    pvals(i) = P;
    Fs(i) = ANOVATAB{2,5};
end


%[~,Ps,~,STATS] = ttest2(tractFAvals(tbivec==0,:),tractFAvals(tbivec>0,:));

%% FA vc FC

[r_within,p_within] = corr(meanFA_intracts',allwithin_subs')%,'type','Spearman')
figure;plot(meanFA_intracts,allwithin_subs,'k.','MarkerSize',20)
%xlim([.2 max(meanFA_intracts)+.1])

% [r_acrosspos,p_acrosspos] = corr(meanFA_intracts',allacrosspos_subs')%,'type','Spearman')
% figure;plot(meanFA_intracts,allacrosspos_subs,'k.','MarkerSize',20)
% xlim([.2 max(meanFA_intracts)+.1])
% 
% [r_acrossneg,p_acrossneg] = corr(meanFA_intracts',allacrossneg_subs')%,'type','Spearman')
% figure;plot(meanFA_intracts,allacrossneg_subs,'k.','MarkerSize',20)
% xlim([.2 max(meanFA_intracts)+.1])



%% FA vs PCL



[r_FA,p_FA] = corr(PCLvec(PCL_subinds),meanFA_intracts(PCL_subinds)')%,'type','Spearman')
prettyscatterplot(meanFA_intracts(PCL_subinds),PCLvec(PCL_subinds),1)
%xlim([.46 .580000001])
export_fig(gca,'FA_vs_PCL.pdf')
% figure('Color','white','position',[1982 478 1352 804]);
% plot(meanFA_intracts(PCL_subinds),PCLvec(PCL_subinds),'k.','MarkerSize',40)
% set(gca,'FontSize',20,'FontWeight','bold','LineWidth',2)
% hold on
% xlim([.45 .6])
% ylim([0 max(PCLvec(PCL_subinds))+1])
% FO = fit(meanFA_intracts(PCL_subinds)', PCLvec(PCL_subinds), 'poly1');
% fitted = [.45 meanFA_intracts(PCL_subinds) .6] .* FO.p1 + FO.p2;
% plot([.45 meanFA_intracts(PCL_subinds) .6],fitted,'r-','LineWidth',3)


%% TBI vs PCL control for FA

[~,T,~,~]=anovan(PCLvec(PCL_subinds)',{tbivec(PCL_subinds),meanFA_intracts(PCL_subinds)'},'continuous',2)


%% FA/FC vs PCL

[~,T,~,~]=anovan(PCLvec(PCL_subinds)',{meanFA_intracts(PCL_subinds)',allwithin_subs(PCL_subinds)'},'continuous',[1 2])

[r,p] = partialcorr([PCLvec(PCL_subinds) meanFA_intracts(PCL_subinds)' allwithin_subs(PCL_subinds)'])

[~,~,PCL_residafterFA] = regress(PCLvec(PCL_subinds),[meanFA_intracts(PCL_subinds)' ones(nnz(PCL_subinds),1)]);
prettyscatterplot(allwithin_subs(PCL_subinds),PCL_residafterFA,1)
export_fig(gca,'RSFC_vs_PCL_controlFA.pdf')

[~,~,PCL_residafterRSFC] = regress(PCLvec(PCL_subinds),[allwithin_subs(PCL_subinds)' ones(nnz(PCL_subinds),1)]);
prettyscatterplot(meanFA_intracts(PCL_subinds),PCL_residafterRSFC,1)
%xlim([.46 .58000001])
export_fig(gca,'FA_vs_PCL_controlRSFC.pdf')

% [~,T,~,~]=anovan(PCLvec(PCL_subinds)',{meanFA_intracts(PCL_subinds)',allacrosspos_subs(PCL_subinds)'},'continuous',[1 2])
% 
% [~,T,~,~]=anovan(PCLvec(PCL_subinds)',{meanFA_intracts(PCL_subinds)',allacrossneg_subs(PCL_subinds)'},'continuous',[1 2])

%% TBI vs PCL control for FA and FC

[~,T,~,~]=anovan(PCLvec(PCL_subinds)',{tbivec(PCL_subinds),meanFA_intracts(PCL_subinds)',allwithin_subs(PCL_subinds)'},'continuous',[2 3])

% [~,T,~,~]=anovan(PCLvec(PCL_subinds)',{tbivec(PCL_subinds),meanFA_intracts(PCL_subinds)',allacrosspos_subs(PCL_subinds)'},'continuous',[2 3])
% 
% [~,T,~,~]=anovan(PCLvec(PCL_subinds)',{tbivec(PCL_subinds),meanFA_intracts(PCL_subinds)',allacrossneg_subs(PCL_subinds)'},'continuous',[2 3])

%% FA/FC vs PCL control for TBI number and severity

[r,p] = partialcorr([PCLvec(PCL_subinds) meanFA_intracts(PCL_subinds)' allwithin_subs(PCL_subinds)' tbivec(PCL_subinds) tbinumvec(PCL_subinds)])
[~,T,~,~]=anovan(PCLvec(PCL_subinds)',{tbivec(PCL_subinds),tbinumvec(PCL_subinds),meanFA_intracts(PCL_subinds)',allwithin_subs(PCL_subinds)'},'continuous',[2 3 4])



%%
% 
% [~,~,~,STATS] = ttest2(allFA(:,tbivec==0)',allFA(:,tbivec>0)');
% STATS.tstat(~logical(tracts.img)) = 0;
% FAdata.img = STATS.tstat';
% save_untouch_nii_2D(FAdata,'FA_ttest_byTBI_T.nii.gz')
% 
% % [~,~,~,STATS] = ttest2(allFA_onsurf(:,tbivec==0)',allFA_onsurf(:,tbivec>0)');
% % FAsurf.data = STATS.tstat';
% % FAsurf.dimord = 'pos_time';
% ft_write_cifti_mod('FAsurf_ttest_byTBI_T',FAsurf)

%  [~,~,~,STATS] = ttest2(allcorticalthickness(:,tbivec==0)',allcorticalthickness(:,tbivec>0)');
%  FAsurf.data = STATS.tstat';
%  FAsurf.dimord = 'pos_time';
%  ft_write_cifti_mod('CorticalThickness_ttest_byTBI_T',FAsurf)










    
