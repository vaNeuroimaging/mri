%%

allparcelfilenames = {'/data/cn4/evan/RestingState/FC_Mapping_120/Cohorts/C1_Lwatershedmerge_0.45.func.gii','/data/cn4/evan/RestingState/FC_Mapping_120/Cohorts/C1_Rwatershedmerge_0.45.func.gii';...
'/data/cn4/evan/ROIs/264_surfvert_ROIs_L.func.gii','/data/cn4/evan/ROIs/264_surfvert_ROIs_R.func.gii';...
'/data/cn4/evan/RestingState/FC_Mapping_120/ConsensusParcels_L.func.gii','/data/cn4/evan/RestingState/FC_Mapping_120/ConsensusParcels_R.func.gii';...
'/data/cn4/evan/ROIs/Craddock/Craddock_350_L.func.gii','/data/cn4/evan/ROIs/Craddock/Craddock_350_R.func.gii';...
'/data/cn4/evan/RestingState/FC_Mapping_120/Yeo_17_L_parcels.func.gii','/data/cn4/evan/RestingState/FC_Mapping_120/Yeo_17_R_parcels.func.gii';...
'/data/cn4/evan/ROIs/mode.L.aparc.a2009s.32k_fs_LR.func.gii','/data/cn4/evan/ROIs/mode.R.aparc.a2009s.32k_fs_LR.func.gii'};

names = {'Edge-derived','Functional','Infomap','NCUTS-derived','Clustering-derived','Anatomical'};

resultsfolder = '/data/cn4/evan/RestingState/FC_Mapping_120/Homogeneity_testing_in_60/';

testslopes = 0;

plot_indfigs = 0;

mediansizesplit = 1;
sizesplit = 200;

allhomogeneities = [];
all_sizes = [];
allparcelationIDs = [];
allrealrandIDs = [];
all_roteigvals = [];

f = 1;

%colors = {'r','g','c','m','b','y'};
colors = {[1 0 0],[0 1 0],[.3 0 .6],[.2 1 1],[0 0 1],[1 1 0]};

order = 1:length(colors);
%order = [5 1 2 3 4 6];

allparcelfilenames = allparcelfilenames(order,:);
colors = colors(order);
names = names(order);

for parcelset = 1:size(allparcelfilenames,1)
    parcelfilenames = allparcelfilenames(parcelset,:);

hems = {'L','R'};

realeigvals_per_first = [];
realsizes = [];


bufsize=16384;
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;


totalparcels = 0;

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


%thishem = load([resultsfolder shortfilename(1:end-9) '.mat']);

PCAmap = gifti([resultsfolder 'PCA_eigval_per_first_' shortfilename]);

clear thishem
for parcelnum = 1:length(parcelIDs)
    
    thishem.realsizes(parcelnum) = nnz(parcels==parcelIDs(parcelnum));
    thishem.realeigvals_per_first(parcelnum) = mean(PCAmap.cdata(parcels==parcelIDs(parcelnum)));
end

realeigvals_per_first = [realeigvals_per_first thishem.realeigvals_per_first];
realsizes = [realsizes thishem.realsizes];


end
%%
% if parcelset==1
%     inds = find(realeigvals_per_first<30);
%     realeigvals_per_first(inds) = [];
%     realsizes(inds) = [];
% end
all_realeigvals_per_first{parcelset} = realeigvals_per_first;
all_realsizes{parcelset} = realsizes;

disp(names{parcelset})
disp(['Mean homogeneity: ' num2str(mean(realeigvals_per_first)) ' +/- ' num2str(std(realeigvals_per_first) / sqrt(length(realeigvals_per_first)))])


sizeindices = logical(realsizes<=sizesplit);
if nnz(sizeindices) > 2
disp(['Mean homogeneity smaller than ' num2str(sizesplit) ': ' num2str(mean(realeigvals_per_first(sizeindices))) ' +/- ' num2str(std(realeigvals_per_first(sizeindices)) / sqrt(length(realeigvals_per_first(sizeindices)))) ', n = ' num2str(nnz(sizeindices))])
end

sizeindices = logical(realsizes>sizesplit);
if nnz(sizeindices) > 2
disp(['Mean homogeneity larger than ' num2str(sizesplit) ': ' num2str(mean(realeigvals_per_first(sizeindices))) ' +/- ' num2str(std(realeigvals_per_first(sizeindices)) / sqrt(length(realeigvals_per_first(sizeindices)))) ', n = ' num2str(nnz(sizeindices))])
end

[parcellationsizes{parcelset} sorti] = sort(realsizes');
 
all_sizes = [all_sizes realsizes];


randadd = ((rand(size(realsizes))*2) -1) * .001;
datain = [(realsizes + randadd)', realeigvals_per_first'];
evalc('[dataout lowerLimit upperLimit xy] = lowess(datain,f,0);');
realfittedvals{parcelset} = xy(:,2);
sortedsizes{parcelset} = xy(:,1);




if plot_indfigs

figure
%plot(reshape(bigsizes,1,numel(bigsizes)),reshape(rotated_eigvals,1,numel(rotated_eigvals)),'.','Color',[.5 .5 .5],'MarkerSize',5)
hold on
plot(realsizes,realeigvals_per_first,'.','MarkerSize',20,'Color',colors{parcelset})

plot(sortedsizes{parcelset},realfittedvals{parcelset},'k-','LineWidth',5)
plot(sortedsizes{parcelset},realfittedvals{parcelset},'-','LineWidth',3,'Color',colors{parcelset})

set(gcf,'Color',[1 1 1])
ylim([30 100])
%set(gca,'XTick',[])
set(gca,'YTick',[])
%export_fig(gcf,[resultsfolder 'RealvRot_fig_' shortfilename(1:end-9) '.pdf'])
end

clear realsizes real_eigvals_per_first

end
disp(' ')


figure
hold on


%plot(allsortedsizes,allrotfittedvals,'k-','LineWidth',7)
for parcelset = 1:size(allparcelfilenames,1)
    plot(sortedsizes{parcelset},realfittedvals{parcelset},'k-','LineWidth',7)
    %plot(sortedsizes{parcelset},realfittedvals{parcelset},[colors{parcelset} '-'],'LineWidth',4)
    plot(sortedsizes{parcelset},realfittedvals{parcelset},'-','LineWidth',4,'Color',colors{parcelset})
end
%legend(names)
set(gcf,'Color',[1 1 1])
xlim([0 3100])






for parcelset1 = 1:size(allparcelfilenames,1)
    for parcelset2 = 1:size(allparcelfilenames,1)
        if parcelset1 < parcelset2
            
            %[p, H, stats] = ranksum(allranks{parcelset1},allranks{parcelset2});
            
            stats=mwwtest(all_realeigvals_per_first{parcelset1},all_realeigvals_per_first{parcelset2});
            disp(['Difference in homogeneities between ' names{parcelset1} ' and ' names{parcelset2} ':'])
            disp(['Homogeneity difference: U = ' num2str(stats.U) ', p = ' num2str(stats.p)])
            
%             stats=mwwtest(allranks{parcelset1},allranks{parcelset2});
%             disp(['Difference in ranks v random between set ' num2str(parcelset1) ' and ' num2str(parcelset2) ':'])
%             disp(['Rank difference: U = ' num2str(stats.U) ', p = ' num2str(stats.p)])
            
            
        end
    end
    disp(' ')
end
            

%[P,T,STATS,TERMS]=anovan(allhomogeneities,{allparcelationIDs allrealrandIDs log10(allsizes)},'continuous',[3],'model','full')
