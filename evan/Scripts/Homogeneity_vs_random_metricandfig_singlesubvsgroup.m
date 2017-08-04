%%

cohortfile = '/data/cn4/evan/RestingState/Ind_variability/Subjects/LFRS_DATALIST.txt';
[subjects surfdatafile] = textread(cohortfile,'%s %s');

f = 1;

for s = 1%:length(subjects)
    subname = subjects{s};
    
    disp(['Subject ' subname])

allparcelfilenames = {['/data/hcp-bluearc/home/laumannt/LFRS_parcellation/' subname '/' subname '_L_watershedmerge_0.45.func.gii'],['/data/hcp-bluearc/home/laumannt/LFRS_parcellation/' subname '/' subname '_R_watershedmerge_0.45.func.gii'];...
    '/data/cn4/evan/RestingState/FC_Mapping_120/120_L_wateredgethresh_watershedmerge_0.45.func.gii','/data/cn4/evan/RestingState/FC_Mapping_120/120_R_wateredgethresh2_watershedmerge_0.45.func.gii'};

rotationresultsfilenames = {['/data/hcp-bluearc/home/laumannt/LFRS_parcellation/' subname '/' subname '_L_.mat'],['/data/hcp-bluearc/home/laumannt/LFRS_parcellation/' subname '/' subname '_R_.mat'];...
    ['/data/hcp-bluearc/home/laumannt/LFRS_parcellation/' subname '/Group_L_.mat'],['/data/hcp-bluearc/home/laumannt/LFRS_parcellation/' subname '/Group_R_.mat']};

resultsfolder = ['/data/hcp-bluearc/home/laumannt/LFRS_parcellation/' subname '/'];


%Poldrome

% allparcelfilenames = {'/data/cn4/evan/RestingState/Ind_variability/Poldrome/Poldrome_L_.01watershedmerge_0.45.func.gii','/data/cn4/evan/RestingState/Ind_variability/Poldrome/Poldrome_R_.01watershedmerge_0.45.func.gii';...
%     '/data/cn4/evan/RestingState/FC_Mapping_120/120_L_wateredgethresh_watershedmerge_0.45.func.gii','/data/cn4/evan/RestingState/FC_Mapping_120/120_R_wateredgethresh2_watershedmerge_0.45.func.gii'};
% 
% rotationresultsfilenames = {'/data/cn4/evan/RestingState/Ind_variability/Poldrome/Poldrome_L_.mat','/data/cn4/evan/RestingState/Ind_variability/Poldrome/Poldrome_R_.mat';...
%     '/data/cn4/evan/RestingState/Ind_variability/Poldrome/Groupparcels_L_.mat','/data/cn4/evan/RestingState/Ind_variability/Poldrome/Groupparcels_R_.mat'};
% 
% resultsfolder = '/data/cn4/evan/RestingState/Ind_variability/Poldrome/';

testslopes = 0;

allhomogeneities = [];
allsizes = [];
allparcelationIDs = [];
allrealrandIDs = [];

for parcelset = 1%:size(allparcelfilenames,1)
    parcelfilenames = allparcelfilenames(parcelset,:);


%parcelfilenames = {'/data/cn4/evan/RestingState/FC_Mapping_120/120_L_wateredgethresh_watershedmerge_0.45.func.gii','/data/cn4/evan/RestingState/FC_Mapping_120/120_R_wateredgethresh2_watershedmerge_0.45.func.gii'};
%parcelfilenames = {'/data/cn4/evan/ROIs/mode.L.aparc.a2009s.32k_fs_LR.func.gii','/data/cn4/evan/ROIs/mode.R.aparc.a2009s.32k_fs_LR.func.gii'};
%parcelfilenames = {'/data/cn4/evan/ROIs/Craddock/Craddock_350_L.func.gii','/data/cn4/evan/ROIs/Craddock/Craddock_350_R.func.gii'};
%parcelfilenames = {'/data/cn4/evan/ROIs/264_surfvert_ROIs_L.func.gii','/data/cn4/evan/ROIs/264_surfvert_ROIs_R.func.gii'};
%parcelfilenames = {'/data/cn4/evan/RestingState/FC_Mapping_120/Random_matchedto_120_L_wateredgethresh_watershedmerge_0.45.func.gii','/data/cn4/evan/RestingState/FC_Mapping_120/Random_matchedto_120_R_wateredgethresh2_watershedmerge_0.45.func.gii'};

hems = {'L','R'};



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

gooddataname = ['/data/cn4/evan/ROIs/all_meanimage_' hem '.func.gii'];
gooddata = gifti(gooddataname); gooddata = gooddata.cdata;
gooddata = gooddata>750;

allgoodindices = find(mask==0 .* gooddata);
allbadvertices = find(logical(mask + (~gooddata)));

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

gooddataname = ['/data/cn4/evan/ROIs/all_meanimage_' hem '.func.gii'];
gooddata = gifti(gooddataname); gooddata = gooddata.cdata;
gooddata = gooddata>750;

allgoodindices = find(mask==0 .* gooddata);
allbadvertices = find(logical(mask + (~gooddata)));

parcelIDs = unique(parcels); parcelIDs(parcelIDs==0) = [];
for parcelID = parcelIDs'
    if nnz(parcels==parcelID) < (nnz(parcels(allbadvertices)==parcelID) + minsize)
        parcels(parcels==parcelID) = 0;
    end
end

parcelIDs = unique(parcels); parcelIDs(parcelIDs==0) = [];

slashes = strfind(parcelfilenames{hemnum},'/');
shortfilename = parcelfilenames{hemnum}(slashes(end)+1:end);
pathname = parcelfilenames{hemnum}(1:slashes(end));
%copyfile([pathname 'PCA_eigval_per_first_' shortfilename],['/data/cn4/evan/RestingState/Ind_variability/Subjects/' subname '/PCA_eigval_per_first_' shortfilename]);


thishem = load(rotationresultsfilenames{parcelset,hemnum});

realeigvals_per_first = [realeigvals_per_first thishem.realeigvals_per_first];
realsizes = [realsizes thishem.realsizes];
rotated_eigvals = [rotated_eigvals; thishem.rotated_eigvals];

outmetric = zeros(32492,1);
ranksoutmetric = zeros(32492,1);

for parcelnum = 1:length(parcelIDs)
    %[ign pval(parcelnum) ign2 Zscore(parcelnum)] = ztest(thishem.realeigvals_per_first(parcelnum),mean(thishem.rotated_eigvals(parcelnum,:)),std(thishem.rotated_eigvals(parcelnum,:)));
    Zscores(end+1) = (thishem.realeigvals_per_first(parcelnum) - mean(thishem.rotated_eigvals(parcelnum,:))) / std(thishem.rotated_eigvals(parcelnum,:));
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
save(gifti(single(ranksoutmetric)),[resultsfolder 'RealvRot_Ranks_' shortfilename])
allranks{parcelset} = [allranks{parcelset} Ranks];
end
%%
bigsizes = repmat(realsizes',1,size(rotated_eigvals,2));
[P,H,STATS] = signrank(allranks{parcelset},size(rotated_eigvals,2)/2);
disp(['Mean rank vs random: ' num2str(mean(allranks{parcelset})) ' +/- ' num2str(std(allranks{parcelset}) / sqrt(length(allranks{parcelset}))) ', signedrank = ' num2str(STATS.signedrank) ', p = ' num2str(P)])

 [parcellationsizes{parcelset} sorti] = sort(realsizes');


%[rotslope, intercept,MSE, R2, S] = logfit(realsizes', mean(rotated_eigvals,2),'logx');
% figure
% [rotslope, realintercept,MSE, R2, S] = logfit(reshape(bigsizes,1,numel(bigsizes))',reshape(rotated_eigvals,1,numel(rotated_eigvals))','logx');
% close(gcf)
% 
% 
%  rotfittedvals = (log10(realsizes)*rotslope + realintercept);
% % 
% figure
% [realslope, rotintercept,MSE, R2, S] = logfit(realsizes', realeigvals_per_first','logx');
% close(gcf)
% % 
%  realfittedvals = (log10(realsizes)*realslope + rotintercept);
 %close all
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
 
 %disp(['Real log slope: ' num2str(realslope)])
 %disp(['Rotated log slope: ' num2str(rotslope)])


 randadd = ((rand(size(realsizes))*2) -1) * .001;
datain = [(realsizes + randadd)', realeigvals_per_first'];
evalc('[dataout lowerLimit upperLimit xy] = lowess(datain,f,0);');
realfittedvals{s} = xy(:,2);
sortedsizes{s} = xy(:,1);

%randadd = ((rand(size(realsizes))*2) -1) * .001;
datain = [(realsizes + randadd)', mean(rotated_eigvals,2)];
evalc('[dataout lowerLimit upperLimit xy] = lowess(datain,f,0);');
rotfittedvals{s} = xy(:,2);



figure
plot(reshape(bigsizes,1,numel(bigsizes)),reshape(rotated_eigvals,1,numel(rotated_eigvals)),'.','Color',[.5 .5 .5],'MarkerSize',5)
hold on
%plot(realsizes,mean(rotated_eigvals,2)','.','MarkerSize',30,'Color',[0 0 1])
plot(realsizes,median(rotated_eigvals,2)','.','MarkerSize',20,'Color',[0 0 0])
plot(realsizes,realeigvals_per_first,'r.','MarkerSize',20)

% [dataout ign ign2 ign3] = lowess([realsizes' mean(rotated_eigvals,2)],.7,0);
% plot(dataout(:,1)',dataout(:,3)','k-','LineWidth',5)
% plot(dataout(:,1)',dataout(:,3)','b-','LineWidth',3)

%FO = fit(realsizes', mean(rotated_eigvals,2), 'poly1');
%FO = fit(realsizes', mean(rotated_eigvals,2), 'poly2');

%fittedvals = ((realsizes.^2)*FO.p1 + realsizes*FO.p2 + FO.p3);
%[sortedsizes sorti] = sort(realsizes);

plot(sortedsizes{s},(rotfittedvals{s}),'k-','LineWidth',5)
%plot(sortedsizes,(rotfittedvals(sorti)),'b-','LineWidth',3)


%FO = fit(realsizes', realeigvals_per_first', 'poly1');
%FO = fit(realsizes', realeigvals_per_first', 'poly2');


%fittedvals = (realsizes*FO.p1 + FO.p2);
%fittedvals = ((realsizes.^2)*FO.p1 + realsizes*FO.p2 + FO.p3);


plot(sortedsizes{s},(realfittedvals{s}),'k-','LineWidth',5)
plot(sortedsizes{s},(realfittedvals{s}),'r-','LineWidth',3)


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
disp(' ')
end

if testslopes

figure
hold on
colors = {'r','b','g','c','m'};
for parcelset = 1:size(allparcelfilenames,1)
    
    fittedvals = equationdiff{parcelset}(1)*log10(parcellationsizes{parcelset});
    
    %plot(parcellationsizes{parcelset},fittedvals,'k-','LineWidth',5)
    plot(parcellationsizes{parcelset},fittedvals,[colors{parcelset} '-'],'LineWidth',5)
end
set(gcf,'Color',[1 1 1])
legend('Subject parcels','Group parcels')
end

% for parcelset1 = 1:size(allparcelfilenames,1)
%     for parcelset2 = 1:size(allparcelfilenames,1)
%         if parcelset1 < parcelset2
%             
%             stats=mwwtest(allranks{parcelset1},allranks{parcelset2});
%             disp(['Difference in ranks v random between set ' num2str(parcelset1) ' and ' num2str(parcelset2) ':'])
%             disp(['Rank difference: U = ' num2str(stats.U) ', p = ' num2str(stats.p)])
%             
%             if testslopes
%             Y = [realvrandomdiffs{parcelset1} ; realvrandomdiffs{parcelset2}];
%             parcellationID = [ones(length(realvrandomdiffs{parcelset1}),1) ; ones(length(realvrandomdiffs{parcelset2}),1)*2];
%             thesesizes = [parcellationsizes{parcelset1} ; parcellationsizes{parcelset2}];
%             X = [parcellationID log10(thesesizes) parcellationID.*log10(thesesizes)];
%             regressionresults = regstats(Y,X);
%             disp(['Difference in slopes v random between set ' num2str(parcelset1) ' and ' num2str(parcelset2) ':'])
%             disp(['Slope difference = ' num2str(equationdiff{parcelset1}(1) - equationdiff{parcelset2}(1)) ' +/- ' num2str(regressionresults.tstat.se(end)) ', t(' num2str(regressionresults.tstat.dfe) ') = ' num2str(regressionresults.tstat.t(end)) ', p = ' num2str(regressionresults.tstat.pval(end))])
%             end
%         end
%     end
%     disp(' ')
% end

end
            

%[P,T,STATS,TERMS]=anovan(allhomogeneities,{allparcelationIDs allrealrandIDs log10(allsizes)},'continuous',[3],'model','full')
