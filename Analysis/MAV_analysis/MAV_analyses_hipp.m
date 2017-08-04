cd /home/data/Analysis/MAV_analysis

initsubjects = textread('MAV_list.txt','%s');
%initsubjects(strcmp(initsubjects,'MAV009')) = [];
initsubjects(strcmp(initsubjects,'MAV019')) = [];
initsubjects(strcmp(initsubjects,'MAV066')) = [];
initsubjects(strcmp(initsubjects,'MAV065')) = [];

%[subjects, tbivec, ages,BAI,BDI,AUDIT,PCL,PSIQ] = textread('MAV_data.txt','%s%n%n%n%n%n%n%n');

[~,~,data]=xlsread('MAV_AssessmentData_cleaned.xlsx');
%[~,~,data]=xlsread('MAV_AssessmentData_cleanedv2.xlsx');
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
DSMBcol_inexcel = strcmp('DSM_B',data(1,:));

tbivec = zeros(length(subjects),1);
tbinumvec = zeros(length(subjects),1);
PCLvec = zeros(length(subjects),1);
DSMBvec = zeros(length(subjects),1);

for s = 1:length(subjects)
    
    excel_subind = strcmp(subjects{s},data(:,1));
    
    if any(excel_subind)
        tbivec(s) = data{excel_subind,tbicol_inexcel};
        tbinumvec(s) = data{excel_subind,tbinumcol_inexcel};
        PCLvec(s) = data{excel_subind,PCLcol_inexcel};
        DSMBvec(s) = data{excel_subind,DSMBcol_inexcel};
        
    else
        tbivec(s) = NaN;
        tbinumvec(s) = NaN;
        PCLvec(s) = NaN;
        DSMBvec(s) = NaN;
    end

   
end

PCL_subinds = ~isnan(PCLvec);

%%





groupparcelfile = '/home/data/atlases/Group_parcellation/Parcels_LR.dlabel.nii';
groupparcelcommunitiesfile = '/home/data/atlases/Group_parcellation/Parcel_Communities.dtseries.nii';
parcels = ft_read_cifti_mod(groupparcelfile);
vmPFCinds = logical((parcels.data==117) + (parcels.data==279));

hipp_vmPFC_fc = zeros(length(subjects),1);

for s = 1:length(subjects)
    
    data = ft_read_cifti_mod(['/home/data/subjects/' subjects{s} '/cifti/cifti_timeseries_normalwall/RSFC_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii']);
    structinds= data.brainstructure(data.brainstructure > 0);
    hipptc = mean(data.data((structinds==11) | (structinds==12) , :),1)';
    vmPFCtc = mean(data.data(vmPFCinds,:),1)';
    hipp_vmPFC_fc(s) = FisherTransform(paircorr_mod(hipptc,vmPFCtc));
end



%% TBI severity vs RSFC

[P,ANOVATAB,STATS] = anova1(hipp_vmPFC_fc,tbivec);
ANOVATAB
 [H,P,CI,STATS] = ttest2([hipp_vmPFC_fc(tbivec==1) ; hipp_vmPFC_fc(tbivec==2)],hipp_vmPFC_fc(tbivec==0));
 STATS
 P


figure('Color','white','position',[1982 478 1352 804]);
bh = boxplot(hipp_vmPFC_fc,tbivec,'colors','k','symbol','k','labels',{'No TBI','Mild TBI','Moderate TBI'});
set(bh,'linewidth',3)
set(bh,'MarkerEdgeColor',[0 0 0])
set(gca,'FontSize',40,'FontWeight','bold','LineWidth',2)
box off
%export_fig(gca,['MeanWithinNetworkConnectivity_byTBI.pdf'])



%% TBI number vs RSFC

test_subinds = tbivec > 0;%true(size(tbivec));%

[r_within,p_within] = corr(tbinumvec(test_subinds),hipp_vmPFC_fc(test_subinds))%,'type','Spearman')
prettyscatterplot(tbinumvec(test_subinds),hipp_vmPFC_fc(test_subinds),1)
%export_fig(gca,'TBInum_vs_WithinRSFC.pdf')



%[r,p] = partialcorr(tbinumvec,hipp_vmPFC_fc,tbivec)


%% RSFC vs PCL



[r_within,p_within] = corr(PCLvec(PCL_subinds),hipp_vmPFC_fc(PCL_subinds))%,'type','Spearman')
prettyscatterplot(hipp_vmPFC_fc(PCL_subinds),PCLvec(PCL_subinds),1)



% %% TBI severity vs PCL
% 
% [P,ANOVATAB,STATS] = anova1(PCLvec(PCL_subinds)',tbivec(PCL_subinds));
% ANOVATAB
% 
% figure('Color','white','position',[1982 478 1352 804]);
% bh = boxplot(PCLvec(PCL_subinds),tbivec(PCL_subinds),'colors','k','symbol','k','labels',{'No TBI','Mild TBI','Moderate TBI'});
% set(bh,'linewidth',3)
% set(bh,'MarkerEdgeColor',[0 0 0])
% set(gca,'FontSize',40,'FontWeight','bold','LineWidth',2)
% set(gca,'YTick',[0:20:80])
% box off
% export_fig(gca,['PCL_byTBI.pdf'])
% 
% [H,P,CI,STATS] = ttest2(PCLvec(tbivec(PCL_subinds)==1)' ,PCLvec(tbivec(PCL_subinds)==0)')
% [H,P,CI,STATS] = ttest2(PCLvec(tbivec(PCL_subinds)==2)' , PCLvec(tbivec(PCL_subinds)==1)')
% 
% 


%% TBI severity vs PCL control for RSFC

[~,T,~,~]=anovan(PCLvec(PCL_subinds)',{tbivec(PCL_subinds),hipp_vmPFC_fc(PCL_subinds)},'continuous',2)



%%  FA
%tracts = load_untouch_nii_2D('/usr/local/fsl/data/atlases/JHU/JHU-ICBM-labels-2mm.nii.gz');
tracts = load_untouch_nii_2D('JHU-ICBM-labels-2mm_selected.nii.gz');

uncinateFA = zeros(length(subjects),1);
for s = 1:length(subjects)
    
    
    FAdata = load_untouch_nii_2D(['/home/data/subjects/' subjects{s} '/DTI/DTI_avg_ec_FA_MNI.nii.gz']);
    
    WMmask = load_untouch_nii_2D(['/home/data/subjects/' subjects{s} '/DTI/aparc+aseg_cerebralwm_mask_222.nii.gz']);
    
    uncinateFA(s) = mean(FAdata.img(logical(WMmask.img) & ((tracts.img==45) | (tracts.img==46))));
    
    
end


%% TBI severity vs FA

[P,ANOVATAB,STATS] = anova1(uncinateFA,tbivec);
ANOVATAB
figure('Color','white','position',[1982 478 1352 804]);
bh = boxplot(uncinateFA,tbivec,'colors','k','symbol','k','labels',{'No TBI','Mild TBI','Moderate TBI'});
set(bh,'linewidth',3)
set(bh,'MarkerEdgeColor',[0 0 0])
set(gca,'FontSize',40,'FontWeight','bold','LineWidth',2)
%ylim([.46 .580001])
box off
%export_fig(gca,['FA_byTBI.pdf'])
% 
% figure;
% plotSpread([meanFA_intracts(tbivec==0)';meanFA_intracts(tbivec==1)';meanFA_intracts(tbivec==2)'],'categoryIdx',[zeros(nnz(tbivec==0),1);ones(nnz(tbivec==1),1);ones(nnz(tbivec==2),1)*2],'categoryColors',{'b','r','g'},'MarkerSize',40)
% set(gcf,'Color',[1 1 1])
% set(gca,'FontSize',15)


%% TBI number vs FA

 [r_within,p_within] = corr(tbinumvec,uncinateFA)%,'type','Spearman')
 %prettyscatterplot(tbinumvec,meanFA_intracts,1)
 %export_fig(gca,'TBInum_vs_FA.pdf')
% figure;plot(tbinumvec,meanFA_intracts,'k.','MarkerSize',20)
% xlim([0 max(tbinumvec)+1])


% test_subinds = tbivec > 0;
% 
% [r_within,p_within] = corr(tbinumvec(test_subinds),meanFA_intracts(test_subinds)')%,'type','Spearman')
% prettyscatterplot(tbinumvec(test_subinds),meanFA_intracts(test_subinds),1)
%  export_fig(gca,'TBInum_vs_FA.pdf')


% %% TBI number vs PCL
% 
% test_subinds = (tbivec > 0) & (~isnan(PCLvec));%true(length(tbivec),1);%
% 
% 
% [r,p] = corr(tbinumvec(test_subinds),PCLvec(test_subinds))%,'type','Spearman')
% prettyscatterplot(tbinumvec(test_subinds),PCLvec(test_subinds),1)
% export_fig(gca,'TBInum_vs_PCL.pdf')
% %figure;plot(tbinumvec(test_subinds),PCLvec(test_subinds),'k.','MarkerSize',20)
% %xlim([0 max(tbinumvec(test_subinds))+1])



%% FA vc FC

[r_within,p_within] = corr(uncinateFA,hipp_vmPFC_fc)%,'type','Spearman')
figure;plot(uncinateFA,hipp_vmPFC_fc,'k.','MarkerSize',20)



%% FA vs PCL

PCL_subinds = ~isnan(PCLvec);

[r_FA,p_FA] = corr(PCLvec(PCL_subinds),uncinateFA(PCL_subinds))%,'type','Spearman')
prettyscatterplot(uncinateFA(PCL_subinds),PCLvec(PCL_subinds),1)
%xlim([.46 .580000001])
%export_fig(gca,'FA_vs_PCL.pdf')
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

[~,T,~,~]=anovan(PCLvec(PCL_subinds)',{tbivec(PCL_subinds),uncinateFA(PCL_subinds)'},'continuous',2)


%% FA/FC vs PCL

[~,T,~,~]=anovan(PCLvec(PCL_subinds)',{uncinateFA(PCL_subinds),hipp_vmPFC_fc(PCL_subinds)},'continuous',[1 2])
% 
% [r,p] = partialcorr([PCLvec(PCL_subinds) meanFA_intracts(PCL_subinds)' allwithin_subs(PCL_subinds)'])
% 
% [~,~,PCL_residafterFA] = regress(PCLvec(PCL_subinds),[meanFA_intracts(PCL_subinds)' ones(nnz(PCL_subinds),1)]);
% prettyscatterplot(allwithin_subs(PCL_subinds),PCL_residafterFA,1)
% export_fig(gca,'RSFC_vs_PCL_controlFA.pdf')
% 
% [~,~,PCL_residafterRSFC] = regress(PCLvec(PCL_subinds),[allwithin_subs(PCL_subinds)' ones(nnz(PCL_subinds),1)]);
% prettyscatterplot(meanFA_intracts(PCL_subinds),PCL_residafterRSFC,1)
% xlim([.46 .58000001])
% export_fig(gca,'FA_vs_PCL_controlRSFC.pdf')

% [~,T,~,~]=anovan(PCLvec(PCL_subinds)',{meanFA_intracts(PCL_subinds)',allacrosspos_subs(PCL_subinds)'},'continuous',[1 2])
% 
% [~,T,~,~]=anovan(PCLvec(PCL_subinds)',{meanFA_intracts(PCL_subinds)',allacrossneg_subs(PCL_subinds)'},'continuous',[1 2])

%% TBI vs PCL control for FA and FC

[~,T,~,~]=anovan(PCLvec(PCL_subinds)',{tbivec(PCL_subinds),uncinateFA(PCL_subinds),hipp_vmPFC_fc(PCL_subinds)},'continuous',[2 3])


%% FA/FC vs PCL control for TBI number and severity

[r,p] = partialcorr([PCLvec(PCL_subinds) uncinateFA(PCL_subinds) hipp_vmPFC_fc(PCL_subinds) tbivec(PCL_subinds) tbinumvec(PCL_subinds)])











    
