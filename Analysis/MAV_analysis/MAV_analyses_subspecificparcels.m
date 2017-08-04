cd /home/data/Analysis/MAV_analysis


subjects = textread('MAV_list.txt','%s');

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



[subjects, tbivec, ages] = textread('MAV_subs_TBIstatus_age.txt','%s%n%n');
for s = 1:length(subjects)
    disp(s)
    corrmat = smartload(['/home/data/subjects/' subjects{s} '/parcellation/parcel_infomap/corrmat.mat']);
    alldata{s} = corrmat;
    parcel_communities{s} = load(['/home/data/subjects/' subjects{s} '/parcellation/parcel_infomap/rawassn_minsize5_regularized_recoloredv4.txt']);
    
    
end

networkmat = zeros(17,17,nnz(tbivec==0));
for s = 1:find(tbivec==0)
    for i = 1:17
        for j = i+1 : 17
            
            thiscomparison_data = alldata{s}(parcel_communities{s}==i,parcel_communities{s}==j);
            if ~isempty(thiscomparison_data)
                networkmat(i,j,s) = mean(thiscomparison_data(:));
            end
        end
    end
end

avg_networkmat = mean(networkmat,3);
acrossindspos = avg_networkmat>0;
acrossindsneg = avg_networkmat<0;


for s = 1:length(subjects)
    within = [];
    acrosspos = [];
    acrossneg = [];
    for i = 1:17
        for j = i : 17
            thiscomparison_data = alldata{s}(parcel_communities{s}==i,parcel_communities{s}==j);
            if i==j
                thiscomparison_data_inds = ~diag(true(size(thiscomparison_data,1),1),0);
                within = [within ; thiscomparison_data(thiscomparison_data_inds)];
            elseif acrossindspos(i,j)
                acrosspos = [acrosspos; thiscomparison_data(:)];
            elseif acrossindsneg(i,j)
                acrossneg = [acrossneg; thiscomparison_data(:)];
            end
        end
    end
    
    allwithin_subs(s) = mean(within);
    allacrosspos_subs(s) = mean(acrosspos);
    allacrossneg_subs(s) = mean(acrossneg);
end
                
                


%%

[H,P,CI,STATS] = ttest2(allwithin_subs(tbivec>0)',allwithin_subs(tbivec==0)')

figure;
plotSpread([allwithin_subs(tbivec==0)';allwithin_subs(tbivec==1)';allwithin_subs(tbivec==2)'],'categoryIdx',[zeros(4,1);ones(9,1);ones(3,1)*2],'categoryColors',{'b','r','g'},'MarkerSize',40)
set(gcf,'Color',[1 1 1])
set(gca,'FontSize',15)
export_fig(gca,['MeanWithinNetworkConnectivity_subspecific.pdf'])




[H,P,CI,STATS] = ttest2(allacrosspos_subs(tbivec>0)',allacrosspos_subs(tbivec==0)')

figure;
plotSpread([allacrosspos_subs(tbivec==0)';allacrosspos_subs(tbivec==1)';allacrosspos_subs(tbivec==2)'],'categoryIdx',[zeros(4,1);ones(9,1);ones(3,1)*2],'categoryColors',{'b','r','g'},'MarkerSize',40)
set(gcf,'Color',[1 1 1])
set(gca,'FontSize',15)
export_fig(gca,['MeanPositiveAcrossNetworkConnectivity_subspecific.pdf'])





[H,P,CI,STATS] = ttest2(allacrossneg_subs(tbivec>0)',allacrossneg_subs(tbivec==0)')

figure;
plotSpread([allacrossneg_subs(tbivec==0)';allacrossneg_subs(tbivec==1)';allacrossneg_subs(tbivec==2)'],'categoryIdx',[zeros(4,1);ones(9,1);ones(3,1)*2],'categoryColors',{'b','r','g'},'MarkerSize',40)
set(gcf,'Color',[1 1 1])
set(gca,'FontSize',15)
export_fig(gca,['MeanNegativeAcrossNetworkConnectivity_subspecific.pdf'])