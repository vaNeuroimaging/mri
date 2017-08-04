cd /home/data/Analysis/MAV_analysis
subjects = textread('MAV_list.txt','%s');
xdist = 30;

%%
communities_120 = load('/home/data/atlases/Group_parcellation/Parcel_Communities_17.txt');
communityIDs = unique(communities_120); communityIDs(communityIDs<1) = [];

corrmat = smartload('/home/data/atlases/Group_parcellation/Parcels_LR_avgcorr_120.mat');

community_relationships = zeros(length(communityIDs));

for IDnum = 1:length(communityIDs)
    for IDnum2 = IDnum : length(communityIDs)
        if IDnum==IDnum2
            community_relationships(IDnum,IDnum2) = 2;
        else
            meancorr = mean(mean(corrmat(communities_120==communityIDs(IDnum),communities_120==communityIDs(IDnum2))));
            community_relationships(IDnum,IDnum2) = sign(meancorr);
            community_relationships(IDnum2,IDnum) = sign(meancorr);
        end
    end
end




%%


%distances = smartload('/home/data/scripts/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Conte69.LR.32k_fs_LR.distances_surfgeo_voleuc_normalwall_xhemlarge_standardsubcort_uint8.mat');

for s = 1:length(subjects)
    disp(subjects{s})
    
    networks = ft_read_cifti_mod(['/home/data/subjects/' subjects{s} '/template_matching/RSFC_templatematch_jaccard_kden0.05.dtseries.nii']);
    networks = networks.data;
    
    community_relationships_full = zeros(size(distances,1),size(distances,2),'int8');
    for IDnum = 1:length(communityIDs)
        for IDnum2 = IDnum : length(communityIDs)
            community_relationships_full(networks==IDnum,networks==IDnum2) = community_relationships(IDnum2,IDnum);
            community_relationships_full(networks==IDnum2,networks==IDnum) = community_relationships(IDnum2,IDnum);
        end
    end
    community_relationships_full(distances<xdist) = 0;
    community_relationships_full = community_relationships_full .* triu(ones(size(community_relationships_full),'int8'),1);
    
    
    data = ft_read_cifti_mod(['/home/data/subjects/' subjects{s} '/cifti/cifti_timeseries_normalwall/RSFC_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii']);
    datacorr = paircorr_mod(data.data');
    clear data
    datacorr(isnan(datacorr)) = 0;
    datacorr = FisherTransform(datacorr);
    
    
    allwithin_subs(s) = mean(datacorr(community_relationships_full==2));
    allacrosspos_subs(s) = mean(datacorr(community_relationships_full==1));
    allacrossneg_subs(s) = mean(datacorr(community_relationships_full==-1));
    
    clear datacorr community_relationships_full
    
end
    
    
%%


[~,~,data]=xlsread('MAV_AssessmentData_cleaned.xlsx');

tbicol_inexcel = strcmp('Vasterling_Severity',data(1,:));
tbinumcol_inexcel = strcmp('Vasterling_Number',data(1,:));

tbivec = zeros(length(subjects),1);
tbinumvec = zeros(length(subjects),1);

for s = 1:length(subjects)
    
    excel_subind = strcmp(subjects{s},data(:,1));
    
    if any(excel_subind)
        tbivec(s) = data{excel_subind,tbicol_inexcel};
        tbinumvec(s) = data{excel_subind,tbinumcol_inexcel};
        
    else
        tbivec(s) = NaN;
        tbinumvec(s) = NaN;
    end

end

    
    
   %%

[P,ANOVATAB,STATS] = anova1(allwithin_subs',tbivec);
ANOVATAB
 [H,P,CI,STATS] = ttest2([allwithin_subs(tbivec==1)' ; allwithin_subs(tbivec==2)'],allwithin_subs(tbivec==0)');
 STATS
 P

figure;
plotSpread([allwithin_subs(tbivec==0)';allwithin_subs(tbivec==1)';allwithin_subs(tbivec==2)'],'categoryIdx',[zeros(nnz(tbivec==0),1);ones(nnz(tbivec==1),1);ones(nnz(tbivec==2),1)*2],'categoryColors',{'b','r','g'},'MarkerSize',40)
set(gcf,'Color',[1 1 1])
set(gca,'FontSize',15)
export_fig(gca,['MeanWithinNetworkConnectivity_subnetworks.pdf'])




[P,ANOVATAB,STATS] = anova1(allacrosspos_subs',tbivec);
ANOVATAB
 [H,P,CI,STATS] = ttest2([allacrosspos_subs(tbivec==1)' ; allacrosspos_subs(tbivec==2)'],allacrosspos_subs(tbivec==0)');
 STATS
 P

figure;
plotSpread([allacrosspos_subs(tbivec==0)';allacrosspos_subs(tbivec==1)';allacrosspos_subs(tbivec==2)'],'categoryIdx',[zeros(nnz(tbivec==0),1);ones(nnz(tbivec==1),1);ones(nnz(tbivec==2),1)*2],'categoryColors',{'b','r','g'},'MarkerSize',40)
set(gcf,'Color',[1 1 1])
set(gca,'FontSize',15)
export_fig(gca,['MeanPositiveAcrossNetworkConnectivity_subnetworks.pdf'])




[P,ANOVATAB,STATS] = anova1(allacrossneg_subs',tbivec);
ANOVATAB
 [H,P,CI,STATS] = ttest2([allacrossneg_subs(tbivec==1)' ; allacrossneg_subs(tbivec==2)'],allacrossneg_subs(tbivec==0)');
 STATS
 P

figure;
plotSpread([allacrossneg_subs(tbivec==0)';allacrossneg_subs(tbivec==1)';allacrossneg_subs(tbivec==2)'],'categoryIdx',[zeros(nnz(tbivec==0),1);ones(nnz(tbivec==1),1);ones(nnz(tbivec==2),1)*2],'categoryColors',{'b','r','g'},'MarkerSize',40)
set(gcf,'Color',[1 1 1])
set(gca,'FontSize',15)
export_fig(gca,['MeanNegativeAcrossNetworkConnectivity_subnetworks.pdf'])
      


% %%
% 
% [H,P,CI,STATS] = ttest2(allwithin_subs(tbivec>0)',allwithin_subs(tbivec==0)')
% 
% figure;
% plotSpread([allwithin_subs(tbivec==0)';allwithin_subs(tbivec==1)';allwithin_subs(tbivec==2)'],'categoryIdx',[zeros(4,1);ones(9,1);ones(3,1)*2],'categoryColors',{'b','r','g'},'MarkerSize',40)
% set(gcf,'Color',[1 1 1])
% set(gca,'FontSize',15)
% export_fig(gca,['MeanWithinNetworkConnectivity_subnetworks.pdf'])
% 
% 
% 
% 
% [H,P,CI,STATS] = ttest2(allacrosspos_subs(tbivec>0)',allacrosspos_subs(tbivec==0)')
% 
% figure;
% plotSpread([allacrosspos_subs(tbivec==0)';allacrosspos_subs(tbivec==1)';allacrosspos_subs(tbivec==2)'],'categoryIdx',[zeros(4,1);ones(9,1);ones(3,1)*2],'categoryColors',{'b','r','g'},'MarkerSize',40)
% set(gcf,'Color',[1 1 1])
% set(gca,'FontSize',15)
% export_fig(gca,['MeanPositiveAcrossNetworkConnectivity_subnetworks.pdf'])
% 
% 
% 
% 
% 
% [H,P,CI,STATS] = ttest2(allacrossneg_subs(tbivec>0)',allacrossneg_subs(tbivec==0)')
% 
% figure;
% plotSpread([allacrossneg_subs(tbivec==0)';allacrossneg_subs(tbivec==1)';allacrossneg_subs(tbivec==2)'],'categoryIdx',[zeros(4,1);ones(9,1);ones(3,1)*2],'categoryColors',{'b','r','g'},'MarkerSize',40)
% set(gcf,'Color',[1 1 1])
% set(gca,'FontSize',15)
% export_fig(gca,['MeanNegativeAcrossNetworkConnectivity_subnetworks.pdf'])