function Homogeneity_values_and_figures(parcellation_names,homogeneity_results,plot_indfigs)
%Homogeneity_values_and_figures(parcellation_names,rotation_results,plot_indfigs)
%
% Determines Z-scores and significance of previously calculated parcellation
% homogeneities relative to matched null model parcellation homogeneities.
% Plots the results using scatter and "beeswarm" graphs.
%
% Inputs:
% parcellation_names - a cell array containing one string per parcellation
%  tested; the strings are the names of the parcellations that will be put
%  onto the plots. E.g., 
%    {'Adult cohort 1','Adult cohort 2','Baby cohort 1','Baby cohort 2'};
% rotation_results - a [n parcellations] x 2 cell array of strings
%  describing the locations of the left and right hemisphere homogeneity
%  testing results on each row. E.g.,
%    {'Adult_1_L.mat','Adult_1_R.mat'; 
%     'Adult_2_L.mat', 'Adult_2_R.mat';
%     'Baby_1_L.mat','Baby_1_R.mat'; 
%     'Baby_2_L.mat','Baby_2_R.mat'}
%  Rows must be in the same order as in parcellation_names.
% plot-indfigs - a binary variable specifying whether to plot a separate
%  scatterplot figure for each parcellation containing the homogeneities of
%  every single real and null model parcel. 1 or omit = yes; 0 = no
%
%EMG 08/26/15

if ischar(parcellation_names)
    parcellation_names = {parcellation_names};
end

if ~exist('plot_indfigs','var')
    plot_indfigs = 1;
end


hems = {'L','R'};

all_sizes = [];
all_roteigvals = [];

f = 1;

for parcelset = 1:size(homogeneity_results,1)
        
    realeigvals_per_first = [];
    realsizes = [];
    rotated_eigvals = [];
    
    
    allranks{parcelset} = [];
    
    Zscores = [];
    Ranks = [];
    for hemnum = 1:length(hems)
        hem = hems{hemnum};
        
        
        thishem = load(homogeneity_results{parcelset,hemnum});
        
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
    
    
    disp(parcellation_names{parcelset})
    bigsizes = repmat(realsizes',1,size(rotated_eigvals,2));
    [P,H,STATS] = signrank(allranks{parcelset},size(rotated_eigvals,2)/2);
    disp(['Mean homogeneity: ' num2str(mean(realeigvals_per_first)) ' +/- ' num2str(std(realeigvals_per_first)) ' (std dev of parcels)'])
    disp(['Mean homogeneity of rotated parcellations: ' num2str(mean(mean(rotated_eigvals,1),2)) ' +/- ' num2str(std(mean(rotated_eigvals,1),[],2)) ' (std dev of parcellation mean homogeneity values)'])
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
        plot(realsizes,realeigvals_per_first,'.','MarkerSize',20,'Color',[1 0 0])
        
                
        plot(sortedsizes{parcelset},rotfittedvals{parcelset},'k-','LineWidth',5)
        
        
        plot(sortedsizes{parcelset},realfittedvals{parcelset},'k-','LineWidth',5)
        plot(sortedsizes{parcelset},realfittedvals{parcelset},'-','LineWidth',3,'Color',[1 0 0])
        
        set(gcf,'Color',[1 1 1])
        ylim([30 100])
        
    end
    disp(' ')
end

figure;
plotSpread(all_meanrotdata,'xNames',parcellation_names,'distributionColors',repmat({'k'},1,length(parcellation_names)),'MarkerSize',10)
hold on
plotSpread(all_meanrealdata,'xNames',parcellation_names,'distributionColors',repmat({'r'},1,length(parcellation_names)),'MarkerSize',50)

set(gcf,'Color',[1 1 1])
set(gca,'FontSize',15)



