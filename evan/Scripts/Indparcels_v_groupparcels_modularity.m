kdenthresholds = [.005:.005:.1];
basedir = '/data/cn4/evan/RestingState/Ind_variability/Subjects/';
datalist = '/data/cn4/evan/RestingState/Ind_variability/Subjects/LFRS_DATALIST.txt';
[subjects subdata] = textread(datalist,'%s%s');

for s = 1:length(subjects)
    
   ind_rawclrs = load([basedir subjects{s} '/Parcels_target240_noratio_infomap/watershed_bothhem_Tk0005to01in0005_S1to1_xd20_BI_INFMAP/rawassn.txt']);
   load([basedir subjects{s} '/Parcels_target240_noratio_infomap/corrmat.mat']);
   
   for i = 1:size(ind_rawclrs,2)
       indQ(s,i) = M_calc_modularity(ind_rawclrs(:,i),all_water_corrmat); 
   end
   
   %ind_corrmat(:,:,s) = all_water_corrmat;
   
   group_rawclrs = load([basedir subjects{s} '/Parcels_group_target240_infomap/watershed_bothhem_Tk0005to01in0005_S1to1_xd20_BI_INFMAP/rawassn.txt']);
   load([basedir subjects{s} '/Parcels_group_target240_infomap/corrmat.mat']);
   
   for i = 1:size(group_rawclrs,2)
       groupQ(s,i) = M_calc_modularity(group_rawclrs(:,i),all_water_corrmat); 
   end
   
   %group_corrmat(:,:,s) = all_water_corrmat;
   
   plot(kdenthresholds,indQ(s,:),'r',kdenthresholds,groupQ(s,:),'b')
   pause
    
   
end

for i = 1:size(ind_rawclrs,2)
     [H,P,CI,STATS] = ttest(indQ(:,i),groupQ(:,i));
     disp(['T = ' num2str(STATS.tstat) ', p = ' num2str(P)])
end

% imagesc(mean(ind_corrmat,3))
% figure
% imagesc(mean(group_corrmat,3))
% figure
plot(kdenthresholds,mean(indQ,1),'r',kdenthresholds,mean(groupQ,1),'b')