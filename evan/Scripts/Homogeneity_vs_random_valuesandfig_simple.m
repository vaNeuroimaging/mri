%%

alltestednames = {'Poldrome_84_subsurf_edge_L_watershedmerge_0.45_tweaked.mat','Poldrome_84_subsurf_edge_R_watershedmerge_0.45_tweaked.mat';
    'Parcels_L.mat','Parcels_R.mat';
    'AAL222_MNI_L.mat','AAL222_MNI_R.mat'};



% {'/data/cn4/evan/RestingState/FC_Mapping_120/Cohorts/C1_Lwatershedmerge_0.45.func.gii','/data/cn4/evan/RestingState/FC_Mapping_120/Cohorts/C1_Rwatershedmerge_0.45.func.gii';...
%     '/data/cn4/evan/ROIs/264_surfvert_ROIs_L.func.gii','/data/cn4/evan/ROIs/264_surfvert_ROIs_R.func.gii';...
%     '/data/cn4/evan/RestingState/FC_Mapping_120/ConsensusParcels_L.func.gii','/data/cn4/evan/RestingState/FC_Mapping_120/ConsensusParcels_R.func.gii';...
%     '/data/cn4/evan/ROIs/Craddock/Craddock_350_L.func.gii','/data/cn4/evan/ROIs/Craddock/Craddock_350_R.func.gii';...
%     '/data/cn4/evan/RestingState/FC_Mapping_120/Yeo_17_L_parcels.func.gii','/data/cn4/evan/RestingState/FC_Mapping_120/Yeo_17_R_parcels.func.gii';...
%     '/data/cn4/evan/ROIs/mode.L.aparc.a2009s.32k_fs_LR.func.gii','/data/cn4/evan/ROIs/mode.R.aparc.a2009s.32k_fs_LR.func.gii'};%...


% allparcelfilenames = {'/data/cn4/evan/RestingState/FC_Mapping_120/120_L_wateredgethresh_watershedmerge_0.45.func.gii','/data/cn4/evan/RestingState/FC_Mapping_120/120_R_wateredgethresh_watershedmerge_0.45.func.gii';...
%     '/data/cn4/evan/RestingState/FC_Mapping_120/120_L_wateredgethresh_watershedmerge_0.25.func.gii','/data/cn4/evan/RestingState/FC_Mapping_120/120_R_wateredgethresh_watershedmerge_0.25.func.gii';...
%     '/data/cn4/evan/RestingState/FC_Mapping_120/120_L_wateredgethresh_watershedmerge_0.65.func.gii','/data/cn4/evan/RestingState/FC_Mapping_120/120_R_wateredgethresh_watershedmerge_0.65.func.gii'};

names = {'Subject','Group','AAL'};
%{'Edge-derived','Functional','Infomap','NCUTS-derived','Clustering-derived','Anatomical'};
%names = {'Edge-derived .45','Edge-derived .25','Edge-derived .65'};


resultsfolder = '/data/cn4/evan/RestingState/Ind_variability/Poldrome/Homogeneity_testing/ReplaceNan';
%'/data/cn4/evan/RestingState/FC_Mapping_120/Homogeneity_testing_in_60/CorrectRot/';
%resultsfolder = '/data/cn4/evan/RestingState/FC_Mapping_120/Homogeneity_testing_in_120/CorrectRot/';

plot_indfigs = 1;


%colors = {'r','g','c','m','b','y'};
colors = {[1 0 0],[0 1 0],[.3 0 .6],[.2 1 1],[0 0 1],[1 1 0]};

order = 1:length(names);
%order = [1 3 5 6 2 4];
%order = [1 2 3];


allhomogeneities = [];
all_sizes = [];
allparcelationIDs = [];
allrealrandIDs = [];
all_roteigvals = [];

f = 1;


alltestednames = alltestednames(order,:);
colors = colors(order);
names = names(order);

for parcelset = 1:size(alltestednames,1)
    
    
    hems = {'L','R'};
    
    realeigvals_per_first = [];
    realsizes = [];
    rotated_eigvals = [];
    
    
    
    totalparcels = 0;
    allranks{parcelset} = [];
    Zscores = [];
    Ranks = [];
    for hemnum = 1:length(hems)
        hem = hems{hemnum};
        
        thishem = load(alltestednames{parcelset,hemnum});
        
        
        realeigvals_per_first = [realeigvals_per_first thishem.realeigvals_per_first];
        realsizes = [realsizes thishem.realsizes];
        rotated_eigvals = [rotated_eigvals; thishem.rotated_eigvals];
        
        
        
        
        for parcelnum = 1:length(thishem.realsizes)
            Zscores(end+1) = (thishem.realeigvals_per_first(parcelnum) - mean(thishem.rotated_eigvals(parcelnum,:))) / std(thishem.rotated_eigvals(parcelnum,:));
            Ranks(end+1) = nnz(thishem.rotated_eigvals(parcelnum,:) < thishem.realeigvals_per_first(parcelnum));
            
        end
        clear H tval
        if hemnum > 1
            allranks{parcelset} = Ranks;
            all_realeigvals_per_first{parcelset} = realeigvals_per_first;
            all_realsizes{parcelset} = realsizes;
        end
    end
    %%
    if strcmp(names{parcelset},'Clustering-derived')
        inds = find(realeigvals_per_first<30);
        realeigvals_per_first(inds) = [];
        realsizes(inds) = [];
        rotated_eigvals(inds,:) = [];
        Ranks(inds) = [];
    end
    
    
    
    disp(names{parcelset})
    bigsizes = repmat(realsizes',1,size(rotated_eigvals,2));
    [P,H,STATS] = signrank(allranks{parcelset},size(rotated_eigvals,2)/2);
    %disp(['Mean homogeneity: ' num2str(mean(realeigvals_per_first)) ' +/- ' num2str(std(realeigvals_per_first) / sqrt(length(realeigvals_per_first)))])
    disp(['Mean homogeneity: ' num2str(mean(realeigvals_per_first)) ' +/- ' num2str(std(realeigvals_per_first))])
    disp(['Mean rank vs random: ' num2str(mean(allranks{parcelset})) ' +/- ' num2str(std(allranks{parcelset}) / sqrt(length(allranks{parcelset}))) ', signedrank = ' num2str(STATS.signedrank) ', p = ' num2str(P)])
    %disp(['Mean rank vs random: ' num2str(mean(allranks{parcelset})) ' +/- ' num2str(std(allranks{parcelset})) ', signedrank = ' num2str(STATS.signedrank) ', p = ' num2str(P)])
    
    
    [parcellationsizes{parcelset} sorti] = sort(realsizes');
    
    all_sizes = [all_sizes realsizes];
    all_roteigvals = [all_roteigvals mean(rotated_eigvals,2)'];
    
    
    
    randadd = ((rand(size(realsizes))*2) -1) * .001;
    datain = [(realsizes + randadd)', realeigvals_per_first'];
    evalc('[dataout lowerLimit upperLimit xy] = lowess(datain,f,0);');
    realfittedvals{parcelset} = xy(:,2);
    sortedsizes{parcelset} = xy(:,1);
    
    %randadd = ((rand(size(realsizes))*2) -1) * .001;
    datain = [(realsizes + randadd)', mean(rotated_eigvals,2)];
    evalc('[dataout lowerLimit upperLimit xy] = lowess(datain,f,0);');
    rotfittedvals{parcelset} = xy(:,2);
    
    %disp(['Real log slope: ' num2str(realslope)])
    %disp(['Rotated log slope: ' num2str(rotslope)])
    
    
    if plot_indfigs
        
        figure
        plot(reshape(bigsizes,1,numel(bigsizes)),reshape(rotated_eigvals,1,numel(rotated_eigvals)),'.','Color',[.5 .5 .5],'MarkerSize',5)
        hold on
        %plot(realsizes,mean(rotated_eigvals,2)','.','MarkerSize',30,'Color',[0 0 1])
        plot(realsizes,median(rotated_eigvals,2)','k.','MarkerSize',20)%,'Color',[0 0 1])
        %plot(realsizes,realeigvals_per_first,[colors{parcelset} '.'],'MarkerSize',20)
        plot(realsizes,realeigvals_per_first,'.','MarkerSize',20,'Color',colors{parcelset})
        
        % [dataout ign ign2 ign3] = lowess([realsizes' mean(rotated_eigvals,2)],.7,0);
        % plot(dataout(:,1)',dataout(:,3)','k-','LineWidth',5)
        % plot(dataout(:,1)',dataout(:,3)','b-','LineWidth',3)
        
        %FO = fit(realsizes', mean(rotated_eigvals,2), 'poly1');
        %FO = fit(realsizes', mean(rotated_eigvals,2), 'poly2');
        
        %fittedvals = ((realsizes.^2)*FO.p1 + realsizes*FO.p2 + FO.p3);
        %[sortedsizes sorti] = sort(realsizes);
        
        %plot(sortedsizes,(rotfittedvals(sorti)),'k-','LineWidth',5)
        %plot(sortedsizes,(rotfittedvals(sorti)),'b-','LineWidth',3)
        
        
        
        plot(sortedsizes{parcelset},rotfittedvals{parcelset},'k-','LineWidth',5)
        %plot(sortedsizes{parcelset},rotfittedvals{parcelset},'b-','LineWidth',3)
        
        
        %FO = fit(realsizes', realeigvals_per_first', 'poly1');
        %FO = fit(realsizes', realeigvals_per_first', 'poly2');
        
        
        %fittedvals = (realsizes*FO.p1 + FO.p2);
        %fittedvals = ((realsizes.^2)*FO.p1 + realsizes*FO.p2 + FO.p3);
        
        
        % plot(sortedsizes,(realfittedvals(sorti)),'k-','LineWidth',5)
        % plot(sortedsizes,(realfittedvals(sorti)),'r-','LineWidth',3)
        
        plot(sortedsizes{parcelset},realfittedvals{parcelset},'k-','LineWidth',5)
        %plot(sortedsizes{parcelset},realfittedvals{parcelset},[colors{parcelset} '-'],'LineWidth',3)
        plot(sortedsizes{parcelset},realfittedvals{parcelset},'-','LineWidth',3,'Color',colors{parcelset})
        
        %set(gca,'Color',[.7 .7 .7])
        set(gcf,'Color',[1 1 1])
        %xlim([0 1500])
        ylim([30 100])
        %xlabel('Parcel Size')
        %ylabel('Parcel Homogeneity')
        %set(gca,'XTick',[])
        %set(gca,'YTick',[])
        %set(gcf,'Position',[1 26 1652 789])
        %pause
        %export_fig(gcf,[resultsfolder 'RealvRot_fig_' shortfilename(1:end-9) '.pdf'])
    end
    disp(' ')
end

f = .5;
randadd = ((rand(size(all_sizes))*2) -1) * .001;
datain = [(all_sizes + randadd)', all_roteigvals'];
evalc('[dataout lowerLimit upperLimit xy] = lowess(datain,f,0);');
allrotfittedvals = xy(:,2);
allsortedsizes = xy(:,1);



figure
hold on

plot(allsortedsizes,allrotfittedvals,'k-','LineWidth',7)
for parcelset = 1:size(alltestednames,1)
    plot(sortedsizes{parcelset},realfittedvals{parcelset},'k-','LineWidth',7)
    %plot(sortedsizes{parcelset},realfittedvals{parcelset},[colors{parcelset} '-'],'LineWidth',4)
    plot(sortedsizes{parcelset},realfittedvals{parcelset},'-','LineWidth',4,'Color',colors{parcelset})
end
%legend(names)
set(gcf,'Color',[1 1 1])
xlim([0 3100])









for parcelset1 = 1:size(alltestednames,1)
    for parcelset2 = 1:size(alltestednames,1)
        if parcelset1 < parcelset2
            
            %[p, H, stats] = ranksum(allranks{parcelset1},allranks{parcelset2});
            
            stats=mwwtest(all_realeigvals_per_first{parcelset1},all_realeigvals_per_first{parcelset2});
            disp(['Difference in homogeneities between ' names{parcelset1} ' and ' names{parcelset2} ':'])
            disp(['Homogeneity difference: U = ' num2str(stats.U) ', p = ' num2str(stats.p)])
            
            stats=mwwtest(allranks{parcelset1},allranks{parcelset2});
            disp(['Difference in ranks v random between ' names{parcelset1} ' and ' names{parcelset2} ':'])
            disp(['Rank difference: U = ' num2str(stats.U) ', p = ' num2str(stats.p)])
            
            
        end
    end
    disp(' ')
end


%[P,T,STATS,TERMS]=anovan(allhomogeneities,{allparcelationIDs allrealrandIDs log10(allsizes)},'continuous',[3],'model','full')
