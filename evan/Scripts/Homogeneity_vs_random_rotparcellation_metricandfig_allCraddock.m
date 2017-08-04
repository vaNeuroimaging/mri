%%

craddocknumbers = [20:10:300 350:50:1000];

for i = 4:42
    
    allparcelfilenames(i-3,:) = {['/data/cn4/evan/ROIs/Craddock/Craddock' sprintf('%04i',i) '_L.func.gii'],['/data/cn4/evan/ROIs/Craddock/Craddock' sprintf('%04i',i) '_R.func.gii']};
    %names{i-3} = ['Craddock' sprintf('%04i',i)];
    names{i-3} = num2str(craddocknumbers(i));
end
% allparcelfilenames(end+1,:) = {'/data/cn4/evan/ROIs/Shen/Shen_200_L.func.gii','/data/cn4/evan/ROIs/Shen/Shen_200_R.func.gii'};
% names{end+1} = 'Shen_200';
% allparcelfilenames(end+1,:) = {'/data/cn4/evan/ROIs/Shen/Shen_100_L.func.gii','/data/cn4/evan/ROIs/Shen/Shen_100_R.func.gii'};
% names{end+1} = 'Shen_100';

allparcelfilenames = {'/data/cn4/evan/ROIs/Shen/Shen_100_L.func.gii','/data/cn4/evan/ROIs/Shen/Shen_100_R.func.gii';...
    '/data/cn4/evan/ROIs/Shen/Shen_200_L.func.gii','/data/cn4/evan/ROIs/Shen/Shen_200_R.func.gii'};
names = {'100','200'};


resultsfolder = '/data/cn4/evan/RestingState/FC_Mapping_120/Homogeneity_testing_in_108/Craddock/ReplaceNan/';

testslopes = 0;

plot_indfigs = 0;

mediansizesplit = 0;
sizesplit = 100;

hems = {'L','R'};

allhomogeneities = [];
all_sizes = [];
allparcelationIDs = [];
allrealrandIDs = [];
all_roteigvals = [];

f = 1;


order = 1:size(allparcelfilenames,1);
%order = [1 2 4 5 3 6 7];
%order = [7 3 6 5 4 2 1];
%order = [1];

allparcelfilenames = allparcelfilenames(order,:);
names = names(order);

for parcelset = 1:size(allparcelfilenames,1)
    parcelfilenames = allparcelfilenames(parcelset,:);
    
    
    %parcelfilenames = {'/data/cn4/evan/RestingState/FC_Mapping_120/120_L_wateredgethresh_watershedmerge_0.45.func.gii','/data/cn4/evan/RestingState/FC_Mapping_120/120_R_wateredgethresh2_watershedmerge_0.45.func.gii'};
    %parcelfilenames = {'/data/cn4/evan/ROIs/mode.L.aparc.a2009s.32k_fs_LR.func.gii','/data/cn4/evan/ROIs/mode.R.aparc.a2009s.32k_fs_LR.func.gii'};
    %parcelfilenames = {'/data/cn4/evan/ROIs/Craddock/Craddock_350_L.func.gii','/data/cn4/evan/ROIs/Craddock/Craddock_350_R.func.gii'};
    %parcelfilenames = {'/data/cn4/evan/ROIs/264_surfvert_ROIs_L.func.gii','/data/cn4/evan/ROIs/264_surfvert_ROIs_R.func.gii'};
    %parcelfilenames = {'/data/cn4/evan/RestingState/FC_Mapping_120/Random_matchedto_120_L_wateredgethresh_watershedmerge_0.45.func.gii','/data/cn4/evan/RestingState/FC_Mapping_120/Random_matchedto_120_R_wateredgethresh2_watershedmerge_0.45.func.gii'};
    
    
    
    
    
    realeigvals_per_first = [];
    realsizes = [];
    rotated_eigvals = [];
    
    
    bufsize=16384;
    caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
    % Read in node neighbor file generated from caret -surface-topology-neighbors
    [neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
        neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
        textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
        'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
    neighbors = neighbors+1;
    
    
    totalparcels = 0;
    allranks{parcelset} = [];
    for hemnum = 1:length(hems)
        hem = hems{hemnum};
        
        
        maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii'];
        mask = gifti(maskname);
        mask = mask.cdata;
        
        minsize = 15;
        
        parcels = gifti(parcelfilenames{hemnum}); parcels = parcels.cdata;
        
%         gooddataname = ['/data/cn4/evan/ROIs/all_meanimage_' hem '.func.gii'];
%         gooddata = gifti(gooddataname); gooddata = gooddata.cdata;
%         gooddata = gooddata>750;
%         
%         allgoodindices = find(mask==0 .* gooddata);
%         allbadvertices = find(logical(mask + (~gooddata)));
        
        baddataname = ['/data/cn4/evan/RestingState/FC_Mapping_120/Homogeneity_testing_in_60/RotParcellation/Baddata_bigcluster_' hem '.func.gii'];
        baddata = gifti(baddataname); baddata = baddata.cdata;
        gooddata = (mask==0) .* (~baddata);
        
        allgoodindices = find((mask==0) .* (~baddata));
        allbadvertices = find(logical(mask + baddata));
        
        parcelIDs = unique(parcels); parcelIDs(parcelIDs==0) = [];
        for parcelID = parcelIDs'
            if nnz(parcels==parcelID) < (nnz(parcels(allbadvertices)==parcelID) + minsize)
                parcels(parcels==parcelID) = 0;
            end
        end
        
        parcelIDs = unique(parcels); parcelIDs(parcelIDs==0) = [];
        
        totalparcels = totalparcels + length(parcelIDs);
    end
    Zscores = [];
    Ranks = [];
    for hemnum = 1:length(hems)
        hem = hems{hemnum};
        
        
        maskname = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii'];
        mask = gifti(maskname);
        mask = mask.cdata;
        
        minsize = 15;
        
        parcels = gifti(parcelfilenames{hemnum}); parcels = parcels.cdata;
        
%         gooddataname = ['/data/cn4/evan/ROIs/all_meanimage_' hem '.func.gii'];
%         gooddata = gifti(gooddataname); gooddata = gooddata.cdata;
%         gooddata = gooddata>750;
%         
%         allgoodindices = find(mask==0 .* gooddata);
%         allbadvertices = find(logical(mask + (~gooddata)));

        baddataname = ['/data/cn4/evan/RestingState/FC_Mapping_120/Homogeneity_testing_in_60/RotParcellation/Baddata_bigcluster_' hem '.func.gii'];
        baddata = gifti(baddataname); baddata = baddata.cdata;
        gooddata = (mask==0) .* (~baddata);
        
        allgoodindices = find((mask==0) .* (~baddata));
        allbadvertices = find(logical(mask + baddata));
        
        parcelIDs = unique(parcels); parcelIDs(parcelIDs==0) = [];
        for parcelID = parcelIDs'
            if nnz(parcels==parcelID) < (nnz(parcels(allbadvertices)==parcelID) + minsize)
                parcels(parcels==parcelID) = 0;
            end
        end
        
        parcelIDs = unique(parcels); parcelIDs(parcelIDs==0) = [];
        
        slashes = strfind(parcelfilenames{hemnum},'/');
        shortfilename = parcelfilenames{hemnum}(slashes(end)+1:end);
        
        
        thishem = load([resultsfolder shortfilename(1:end-9) '.mat']);
        
        realeigvals_per_first = [realeigvals_per_first thishem.realeigvals_per_first];
        realsizes = [realsizes thishem.realsizes];
        rotated_eigvals = [rotated_eigvals; thishem.rotated_eigvals];
        
        outmetric = zeros(32492,1);
        ranksoutmetric = zeros(32492,1);
        
        for parcelnum = 1:length(parcelIDs)
            %[ign pval(parcelnum) ign2 Zscore(parcelnum)] = ztest(thishem.realeigvals_per_first(parcelnum),mean(thishem.rotated_eigvals(parcelnum,:)),std(thishem.rotated_eigvals(parcelnum,:)));
            Zscores(end+1) = (thishem.realeigvals_per_first(parcelnum) - nanmean(thishem.rotated_eigvals(parcelnum,:))) / std(thishem.rotated_eigvals(parcelnum,:));
            Ranks(end+1) = nnz(thishem.rotated_eigvals(parcelnum,:) < thishem.realeigvals_per_first(parcelnum));
            %[H(parcelnum) ign ign2 stats] = ttest(thishem.rotated_eigvals(parcelnum,:),thishem.realeigvals_per_first(parcelnum),(.05/totalparcels));
            %tval(parcelnum) = stats.tstat;
            
            %if H(parcelnum)==1
            outmetric(parcels==parcelIDs(parcelnum)) = Zscores(end);
            ranksoutmetric(parcels==parcelIDs(parcelnum)) = Ranks(end) - 500;
            %     else
            %         verts = find(parcels==parcelIDs(parcelnum));
            %         for vert = verts(:)'
            %             vertneighs = neighbors(vert,2:7); vertneighs(isnan(vertneighs))=[];
            %             if any(parcels(vertneighs) ~= parcels(vert))
            %                 outmetric(vert) = .01;
            %             end
            %         end
            %     end
        end
        clear H tval
        %save(gifti(single(outmetric)),[resultsfolder 'RealvRot_Zs_' shortfilename])
        %save(gifti(single(ranksoutmetric)),[resultsfolder 'RealvRot_Ranks_' shortfilename])
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
    end
    
    disp(names{parcelset})
    bigsizes = repmat(realsizes',1,size(rotated_eigvals,2));
    [P,H,STATS] = signrank(allranks{parcelset},size(rotated_eigvals,2)/2);
    disp(['Mean homogeneity: ' num2str(mean(realeigvals_per_first)) ' +/- ' num2str(std(realeigvals_per_first)) ' (std dev), +/- ' num2str(std(realeigvals_per_first) / sqrt(length(realeigvals_per_first))) ' (std err)'])
    %disp(['Mean rank vs random: ' num2str(mean(allranks{parcelset})) ' +/- ' num2str(std(allranks{parcelset})) ', signedrank = ' num2str(STATS.signedrank) ', p = ' num2str(P)])
    disp(['Mean homogeneity of rotated parcellations: ' num2str(mean(mean(rotated_eigvals,1),2)) ' +/- ' num2str(std(mean(rotated_eigvals,1),[],2))])
    disp(['Z score of homogeneity vs random: ' num2str((mean(realeigvals_per_first) - mean(mean(rotated_eigvals,1),2)) / std(mean(rotated_eigvals,1),[],2))])
    disp(['Homogeneity is better than ' num2str(nnz(mean(rotated_eigvals,1) < mean(realeigvals_per_first))) ' of ' num2str(size(rotated_eigvals,2)) ' random parcellations'])
    
%     if mediansizesplit
%         sizesplit = median(realsizes);
%     end
%     
%     sizeindices = logical(realsizes<=sizesplit);
%     if nnz(sizeindices) > 2
%         [P,H,STATS] = signrank(allranks{parcelset}(sizeindices),size(rotated_eigvals,2)/2);
%         % disp(['Mean homogeneity smaller than ' num2str(sizesplit) ': ' num2str(mean(realeigvals_per_first(sizeindices))) ' +/- ' num2str(std(realeigvals_per_first(sizeindices)) / sqrt(length(realeigvals_per_first(sizeindices)))) ', n = ' num2str(nnz(sizeindices))])
%         %disp(['Mean rank vs random smaller than ' num2str(sizesplit) ': ' num2str(mean(allranks{parcelset}(sizeindices))) ' +/- ' num2str(std(allranks{parcelset}(sizeindices)) / sqrt(length(allranks{parcelset}(sizeindices)))) ', signedrank = ' num2str(STATS.signedrank) ', p = ' num2str(P)])
%         disp(['Mean homogeneity smaller than ' num2str(sizesplit) ': ' num2str(mean(realeigvals_per_first(sizeindices))) ' +/- ' num2str(std(realeigvals_per_first(sizeindices))) ', n = ' num2str(nnz(sizeindices))])
%         disp(['Mean rank vs random smaller than ' num2str(sizesplit) ': ' num2str(mean(allranks{parcelset}(sizeindices))) ' +/- ' num2str(std(allranks{parcelset}(sizeindices))) ', signedrank = ' num2str(STATS.signedrank) ', p = ' num2str(P)])
%     end
%     
%     sizeindices = logical(realsizes>sizesplit);
%     if nnz(sizeindices) > 2
%         [P,H,STATS] = signrank(allranks{parcelset}(sizeindices),size(rotated_eigvals,2)/2);
%         % disp(['Mean homogeneity larger than ' num2str(sizesplit) ': ' num2str(mean(realeigvals_per_first(sizeindices))) ' +/- ' num2str(std(realeigvals_per_first(sizeindices)) / sqrt(length(realeigvals_per_first(sizeindices)))) ', n = ' num2str(nnz(sizeindices))])
%         % disp(['Mean rank vs random larger than ' num2str(sizesplit) ': ' num2str(mean(allranks{parcelset}(sizeindices))) ' +/- ' num2str(std(allranks{parcelset}(sizeindices)) / sqrt(length(allranks{parcelset}(sizeindices)))) ', signedrank = ' num2str(STATS.signedrank) ', p = ' num2str(P)])
%         disp(['Mean homogeneity larger than ' num2str(sizesplit) ': ' num2str(mean(realeigvals_per_first(sizeindices))) ' +/- ' num2str(std(realeigvals_per_first(sizeindices))) ', n = ' num2str(nnz(sizeindices))])
%         disp(['Mean rank vs random larger than ' num2str(sizesplit) ': ' num2str(mean(allranks{parcelset}(sizeindices))) ' +/- ' num2str(std(allranks{parcelset}(sizeindices))) ', signedrank = ' num2str(STATS.signedrank) ', p = ' num2str(P)])
%     end
    
    [parcellationsizes{parcelset} sorti] = sort(realsizes');
    
    all_sizes = [all_sizes realsizes];
    all_roteigvals = [all_roteigvals nanmean(rotated_eigvals,2)'];
    
    % figure
    % [rotslope, realintercept,MSE, R2, S] = logfit(reshape(bigsizes,1,numel(bigsizes))',reshape(rotated_eigvals,1,numel(rotated_eigvals))','logx');
    % close(gcf)
    %
    %
    % rotfittedvals = (log10(realsizes)*rotslope + realintercept);
    % %
    % figure
    % [realslope, rotintercept,MSE, R2, S] = logfit(realsizes', realeigvals_per_first','logx');
    % close(gcf)
    % %
    %  realfittedvals = (log10(realsizes)*realslope + rotintercept);
    %  %close all
    %  if testslopes
    %  realslopes{parcelset} = realslope;
    %  rotslopes{parcelset} = rotslope;
    %
    %  realvrandomdiffs{parcelset} = realeigvals_per_first(sorti)' - mean(rotated_eigvals(sorti),2);
    %  equationdiff{parcelset} = [realslope - rotslope , realintercept - rotintercept];
    %
    %  allhomogeneities = [realeigvals_per_first' ; mean(rotated_eigvals,2)];
    % allsizes = [realsizes' ; realsizes'];
    % %allparcelationIDs = [allparcelationIDs ; repmat(parcelset,length(realsizes),1) ; repmat(parcelset,numel(bigsizes),1)];
    % allrealrandIDs = [ones(numel(realsizes),1) ; ones(numel(realsizes),1)*2];
    %
    % regressionresults = regstats(allhomogeneities,[allrealrandIDs log10(allsizes) log10(allsizes).*allrealrandIDs]);
    % disp(['Difference in slopes: Slope difference = ' num2str(realslope - rotslope) ' +/- ' num2str(regressionresults.tstat.se(end)) ', t(' num2str(regressionresults.tstat.dfe) ') = ' num2str(regressionresults.tstat.t(end)) ', p = ' num2str(regressionresults.tstat.pval(end))])
    %
    %
    % end
    
    randadd = ((rand(size(realsizes))*2) -1) * .001;
    datain = [(realsizes + randadd)', realeigvals_per_first'];
    evalc('[dataout lowerLimit upperLimit xy] = lowess(datain,f,0);');
    realfittedvals{parcelset} = xy(:,2);
    sortedsizes{parcelset} = xy(:,1);
    
    %randadd = ((rand(size(realsizes))*2) -1) * .001;
    datain = [(realsizes + randadd)', nanmean(rotated_eigvals,2)];
    datain(isnan(datain)) = nanmean(datain(:,2));
    evalc('[dataout lowerLimit upperLimit xy] = lowess(datain,f,0);');
    rotfittedvals{parcelset} = xy(:,2);
    
    %disp(['Real log slope: ' num2str(realslope)])
    %disp(['Rotated log slope: ' num2str(rotslope)])
    
    
    
    if plot_indfigs
        
        figure
        plot(reshape(bigsizes,1,numel(bigsizes)),reshape(rotated_eigvals,1,numel(rotated_eigvals)),'.','Color',[.5 .5 .5],'MarkerSize',5)
        hold on
        %plot(realsizes,mean(rotated_eigvals,2)','.','MarkerSize',30,'Color',[0 0 1])
        plot(realsizes,nanmedian(rotated_eigvals,2)','k.','MarkerSize',20)%,'Color',[0 0 1])
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
        export_fig(gcf,[resultsfolder 'RealvRot_fig_' shortfilename(1:end-9) '.pdf'])
    end
    disp(' ')
end

figure;
plotSpread(all_meanrotdata,'xNames',names,'distributionColors',repmat({'k'},1,length(names)),'MarkerSize',10)
hold on
%plotSpread(all_meanrealdata,'xNames',names,'distributionColors',colors)
plotSpread(all_meanrealdata,'xNames',names,'distributionColors',repmat({'r'},1,length(names)),'MarkerSize',50)

set(gcf,'Color',[1 1 1])
set(gca,'FontSize',20)


f = .5;
randadd = ((rand(size(all_sizes))*2) -1) * .001;
datain = [(all_sizes + randadd)', all_roteigvals'];
datain(isnan(datain)) = nanmean(datain(:,2));
evalc('[dataout lowerLimit upperLimit xy] = lowess(datain,f,0);');
allrotfittedvals = xy(:,2);
allsortedsizes = xy(:,1);


% if testslopes
%
% figure
% hold on
% colors = {'r','b','g','c','m'};
% for parcelset = 1:size(allparcelfilenames,1)
%
%     fittedvals = equationdiff{parcelset}(1)*log10(parcellationsizes{parcelset});
%
%     %plot(parcellationsizes{parcelset},fittedvals,'k-','LineWidth',5)
%     %plot(parcellationsizes{parcelset},fittedvals,[colors{parcelset} '-'],'LineWidth',5)
% end
% %set(gcf,'Color',[1 1 1])
% %legend('Edge-derived','Anatomical','NCUTS-derived','Functional','Random')
% end

figure
hold on

% for parcelset = 1:size(allparcelfilenames,1)
%     %plot(sortedsizes{parcelset},realfittedvals{parcelset},'k-','LineWidth',5)
%     plot(all_realsizes{parcelset},all_realeigvals_per_first{parcelset},[colors{parcelset} '.'],'MarkerSize',30)
% end

% plot(allsortedsizes,allrotfittedvals,'k-','LineWidth',7)
% for parcelset = 1:size(allparcelfilenames,1)
%     plot(sortedsizes{parcelset},realfittedvals{parcelset},'k-','LineWidth',7)
%     %plot(sortedsizes{parcelset},realfittedvals{parcelset},[colors{parcelset} '-'],'LineWidth',4)
%     plot(sortedsizes{parcelset},realfittedvals{parcelset},'-','LineWidth',4,'Color',colors{parcelset})
% end
% %legend(names)
% set(gcf,'Color',[1 1 1])
% xlim([0 3100])



% figure
% hold on
%
% % for parcelset = 1:size(allparcelfilenames,1)
% %     %plot(sortedsizes{parcelset},realfittedvals{parcelset},'k-','LineWidth',5)
% %     plot(all_realsizes{parcelset},all_realeigvals_per_first{parcelset},[colors{parcelset} '.'],'MarkerSize',30)
% % end
%
%
% for parcelset = 1:size(allparcelfilenames,1)
%     %plot(sortedsizes{parcelset},rotfittedvals{parcelset},[colors{parcelset} '-'],'LineWidth',6)
%     plot(sortedsizes{parcelset},rotfittedvals{parcelset},'-','LineWidth',6,'Color',colors{parcelset})
%
% end
% namesplus = names; namesplus{end+1} = 'Random';
% plot(allsortedsizes,allrotfittedvals,'k-','LineWidth',7)
% legend(namesplus)
% set(gcf,'Color',[1 1 1])
% xlim([0 4000])



% 
% 
% 
% for parcelset1 = 1:size(allparcelfilenames,1)
%     for parcelset2 = 1:size(allparcelfilenames,1)
%         if parcelset1 < parcelset2
%             
%             %[p, H, stats] = ranksum(allranks{parcelset1},allranks{parcelset2});
%             
%             stats=mwwtest(all_realeigvals_per_first{parcelset1},all_realeigvals_per_first{parcelset2});
%             disp(['Difference in homogeneities between ' names{parcelset1} ' and ' names{parcelset2} ':'])
%             disp(['Homogeneity difference: U = ' num2str(stats.U) ', p = ' num2str(stats.p)])
%             
%             %             stats=mwwtest(allranks{parcelset1},allranks{parcelset2});
%             %             disp(['Difference in ranks v random between set ' num2str(parcelset1) ' and ' num2str(parcelset2) ':'])
%             %             disp(['Rank difference: U = ' num2str(stats.U) ', p = ' num2str(stats.p)])
%             
%             if testslopes
%                 Y = [realvrandomdiffs{parcelset1} ; realvrandomdiffs{parcelset2}];
%                 parcellationID = [ones(length(realvrandomdiffs{parcelset1}),1) ; ones(length(realvrandomdiffs{parcelset2}),1)*2];
%                 thesesizes = [parcellationsizes{parcelset1} ; parcellationsizes{parcelset2}];
%                 X = [parcellationID log10(thesesizes) parcellationID.*log10(thesesizes)];
%                 regressionresults = regstats(Y,X);
%                 disp(['Difference in slopes v random between set ' num2str(parcelset1) ' and ' num2str(parcelset2) ':'])
%                 disp(['Slope difference = ' num2str(equationdiff{parcelset1}(1) - equationdiff{parcelset2}(1)) ' +/- ' num2str(regressionresults.tstat.se(end)) ', t(' num2str(regressionresults.tstat.dfe) ') = ' num2str(regressionresults.tstat.t(end)) ', p = ' num2str(regressionresults.tstat.pval(end))])
%             end
%         end
%     end
%     disp(' ')
% end


%[P,T,STATS,TERMS]=anovan(allhomogeneities,{allparcelationIDs allrealrandIDs log10(allsizes)},'continuous',[3],'model','full')
