%%

alltestednames = {'C1_Baby_28sub_nosmooth_L_threshperc0.45_minparcel20watershedmerge_0.45.mat','C1_Baby_28sub_nosmooth_R_threshperc0.45_minparcel20watershedmerge_0.45.mat';...
    'C2_Baby_28sub_nosmooth_L_threshperc0.45_minparcel20watershedmerge_0.45.mat','C2_Baby_28sub_nosmooth_R_threshperc0.45_minparcel20watershedmerge_0.45.mat';...
    'C1_adult_parcels_L.mat','C1_adult_parcels_R.mat';...
    'C2_adult_parcels_L.mat','C2_adult_parcels_R.mat'};


%{'Baby_28sub_nosmooth_L_threshperc0.45_minparcel20watershedmerge_0.45.mat','Baby_28sub_nosmooth_R_threshperc0.45_minparcel20watershedmerge_0.45.mat';...
%    '/data/cn4/evan/RestingState/FC_Mapping_120/Homogeneity_testing_in_120/ReplaceNan/120_108_combined_watershedmerge_0.35_tweaked.mat','/data/cn4/evan/RestingState/FC_Mapping_120/Homogeneity_testing_in_120/ReplaceNan/120_108_combined_R_watershedmerge_0.35_tweaked.mat'};
%{'Baby_28sub_nosmooth_L_threshperc0.45_minparcel20watershedmerge_0.45.mat','Baby_28sub_nosmooth_R_threshperc0.45_minparcel20watershedmerge_0.45.mat'};
%{'C1_Baby_28sub_nosmooth_L_threshperc0.45_minparcel20watershedmerge_0.45.mat','C1_Baby_28sub_nosmooth_R_threshperc0.45_minparcel20watershedmerge_0.45.mat';...
%    'C2_Baby_28sub_nosmooth_L_threshperc0.45_minparcel20watershedmerge_0.45.mat','C2_Baby_28sub_nosmooth_R_threshperc0.45_minparcel20watershedmerge_0.45.mat'};
%{'Parcels_avg_corrmat_edgethresh_0.35_L.mat','Parcels_avg_corrmat_edgethresh_0.35_R.mat'};
%{'Baby_28sub_nosmooth_L_threshperc0.45_minparcel20watershedmerge_0.45.mat','Baby_28sub_nosmooth_R_threshperc0.45_minparcel20watershedmerge_0.45.mat';...
%     'Parcels_L.mat','Parcels_R.mat';...
%     'AAL222_MNI_L.mat','AAL222_MNI_R.mat'};

names = {'C1Babies','C2Babies','C1Adults','C2Adults'};
%{'C1','C2'};
%{'Babies','Adults','AAL'};
%{'Boundary map','Power ROIs','Craddock','Shen','Power Communities','Infomap','Yeo','Brodmann','AAL'};



resultsfolder = pwd;

plot_indfigs = 0;

hems = {'L','R'};

allhomogeneities = [];
all_sizes = [];
allparcelationIDs = [];
allrealrandIDs = [];
all_roteigvals = [];

f = 1;

%colors = {'r','g','c','m','b','y'};
colors = {[1 0 0],[0 1 0],[.3 0 .6],[.2 1 1]};%,[1 .5 0],[0 0 1],[1 1 0],[1 .6 1],[1 1 1]};

order = 1:size(alltestednames,1);
%order = [1:7];
%order = [6 5 4 3 2 1];
%order = [1];

alltestednames = alltestednames(order,:);
colors = colors(order);
names = names(order);



for parcelset = 1:size(alltestednames,1)
    parcelfilenames = alltestednames(parcelset,:);
    
    
    realeigvals_per_first = [];
    realsizes = [];
    rotated_eigvals = [];
    
    
    allranks{parcelset} = [];
    
    Zscores = [];
    Ranks = [];
    for hemnum = 1:length(hems)
        hem = hems{hemnum};
        
        
        thishem = load(alltestednames{parcelset,hemnum});
        
        for k = 1:size(thishem.rotated_eigvals,1)
            naninds = logical(isnan(thishem.rotated_eigvals(k,:)));
            thishem.rotated_eigvals(k,naninds) = nanmean(thishem.rotated_eigvals(k,:));
        end
        
        realeigvals_per_first = [realeigvals_per_first thishem.realeigvals_per_first];
        realsizes = [realsizes thishem.realsizes];
        rotated_eigvals = [rotated_eigvals; thishem.rotated_eigvals];
        
        
        for parcelnum = 1:length(thishem.realsizes)
            Zscores(end+1) = (thishem.realeigvals_per_first(parcelnum) - nanmean(thishem.rotated_eigvals(parcelnum,:))) / std(thishem.rotated_eigvals(parcelnum,:));
            Ranks(end+1) = nnz(thishem.rotated_eigvals(parcelnum,:) < thishem.realeigvals_per_first(parcelnum));
           
        end
        clear H tval
        
        allranks{parcelset} = Ranks;
        all_realeigvals_per_first{parcelset} = realeigvals_per_first;
        all_realsizes{parcelset} = realsizes;
        
        
        
    end
    all_meanrotdata(:,parcelset) = nanmean(rotated_eigvals,1)';
    all_meanrealdata(1,parcelset) = mean(realeigvals_per_first);
    
    %%
    if strcmp(names{parcelset},'Clustering-derived')
        inds = find(realeigvals_per_first<30);
        realeigvals_per_first(inds) = [];
        realsizes(inds) = [];
        rotated_eigvals(inds,:) = [];
        Zscores(inds) = [];
    end
    
    disp(names{parcelset})
    bigsizes = repmat(realsizes',1,size(rotated_eigvals,2));
    [P,H,STATS] = signrank(allranks{parcelset},size(rotated_eigvals,2)/2);
    disp(['Mean homogeneity: ' num2str(mean(realeigvals_per_first)) ' +/- ' num2str(std(realeigvals_per_first)) ' (std dev), +/- ' num2str(std(realeigvals_per_first) / sqrt(length(realeigvals_per_first))) ' (std err)'])
    disp(['Mean homogeneity of rotated parcellations: ' num2str(mean(mean(rotated_eigvals,1),2)) ' +/- ' num2str(std(mean(rotated_eigvals,1),[],2))])
    disp(['Z score of homogeneity vs random: ' num2str((mean(realeigvals_per_first) - mean(mean(rotated_eigvals,1),2)) / std(mean(rotated_eigvals,1),[],2))])
    disp(['Homogeneity is better than ' num2str(nnz(mean(rotated_eigvals,1) < mean(realeigvals_per_first))) ' of ' num2str(size(rotated_eigvals,2)) ' random parcellations'])
    

    [parcellationsizes{parcelset} sorti] = sort(realsizes');
    
    all_sizes = [all_sizes realsizes];
    all_roteigvals = [all_roteigvals nanmean(rotated_eigvals,2)'];
    
    
    
    randadd = ((rand(size(realsizes))*2) -1) * .001;
    datain = [(realsizes + randadd)', realeigvals_per_first'];
    evalc('[dataout lowerLimit upperLimit xy] = lowess(datain,f,0);');
    realfittedvals{parcelset} = xy(:,2);
    sortedsizes{parcelset} = xy(:,1);
    
    datain = [(realsizes + randadd)', nanmean(rotated_eigvals,2)];
    datain(isnan(datain)) = nanmean(datain(:,2));
    evalc('[dataout lowerLimit upperLimit xy] = lowess(datain,f,0);');
    rotfittedvals{parcelset} = xy(:,2);
    
    
    
    
    if plot_indfigs
        
        figure
        plot(reshape(bigsizes,1,numel(bigsizes)),reshape(rotated_eigvals,1,numel(rotated_eigvals)),'.','Color',[.5 .5 .5],'MarkerSize',5)
        hold on
        plot(realsizes,nanmedian(rotated_eigvals,2)','k.','MarkerSize',20)%,'Color',[0 0 1])
        plot(realsizes,realeigvals_per_first,'.','MarkerSize',20,'Color',colors{parcelset})
        
                
        plot(sortedsizes{parcelset},rotfittedvals{parcelset},'k-','LineWidth',5)
        
        
        plot(sortedsizes{parcelset},realfittedvals{parcelset},'k-','LineWidth',5)
        plot(sortedsizes{parcelset},realfittedvals{parcelset},'-','LineWidth',3,'Color',colors{parcelset})
        
        set(gcf,'Color',[1 1 1])
        ylim([30 100])
        
    end
    disp(' ')
end

figure;
plotSpread(all_meanrotdata,'xNames',names,'distributionColors',repmat({'k'},1,length(names)),'MarkerSize',10)
hold on
plotSpread(all_meanrealdata,'xNames',names,'distributionColors',repmat({'r'},1,length(names)),'MarkerSize',50)

set(gcf,'Color',[1 1 1])
set(gca,'FontSize',15)



