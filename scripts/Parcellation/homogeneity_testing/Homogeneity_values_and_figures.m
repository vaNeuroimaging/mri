function Homogeneity_values_and_figures(homogeneity_results,plot_indfigs)
%Homogeneity_values_and_figures(rotation_results,plot_indfigs)
%
% Determines Z-scores and significance of previously calculated parcellation
% homogeneities relative to matched null model parcellation homogeneities.
% Plots the results using scatter and "beeswarm" graphs.
%
% Inputs:
% rotation_results - a cell array of strings describing the locations of
%  the homogeneity testing results.
% plot-indfigs - a binary variable specifying whether to plot a separate
%  scatterplot figure for each parcellation containing the homogeneities of
%  every single real and null model parcel. 1 = yes; 0 or omit = no
%
%EMG 08/26/15


if ~exist('plot_indfigs','var')
    plot_indfigs = 0;
end

if ischar(homogeneity_results)
    homogeneity_results = {homogeneity_results};
end


all_sizes = [];
all_roteigvals = [];

f = 1;

for parcelset = 1:length(homogeneity_results)
    
    [~,name,~] = fileparts(homogeneity_results{parcelset});
    name(strfind(name,'_')) = ' ';
    names{parcelset} = name;
        
    allranks{parcelset} = [];
    
    load(homogeneity_results{parcelset});
    
    Zscores = zeros(size(rotated_eigvals,1),1);
    Ranks = zeros(size(rotated_eigvals,1),1);
    
    for parcelnum = 1:size(rotated_eigvals,1)
        naninds = logical(isnan(rotated_eigvals(parcelnum,:)));
        rotated_eigvals(parcelnum,naninds) = nanmean(rotated_eigvals(parcelnum,:));
        
        Zscores(parcelnum) = (realeigvals_per_first(parcelnum) - nanmean(rotated_eigvals(parcelnum,:))) / std(rotated_eigvals(parcelnum,:));
        Ranks(parcelnum) = nnz(rotated_eigvals(parcelnum,:) < realeigvals_per_first(parcelnum));
    end
    
    allranks{parcelset} = Ranks;
    all_realeigvals_per_first{parcelset} = realeigvals_per_first;
    all_realsizes{parcelset} = realsizes;
    
    
    all_meanrotdata(:,parcelset) = nanmean(rotated_eigvals,1)';
    all_meanrealdata(1,parcelset) = mean(realeigvals_per_first);
    
    %%
    
    
    disp(names{parcelset})
    bigsizes = repmat(realsizes',1,size(rotated_eigvals,2));
    [P,H,STATS] = signrank(allranks{parcelset},size(rotated_eigvals,2)/2);
    disp(['Mean homogeneity: ' num2str(mean(realeigvals_per_first)) ' +/- ' num2str(std(realeigvals_per_first)) ' (std dev of parcels)'])
    disp(['Mean homogeneity of rotated parcellations: ' num2str(mean(mean(rotated_eigvals,1),2)) ' +/- ' num2str(std(mean(rotated_eigvals,1),[],2)) ' (std dev of parcellation mean homogeneity values)'])
    disp(['Z score of homogeneity vs random: ' num2str((mean(realeigvals_per_first) - mean(mean(rotated_eigvals,1),2)) / std(mean(rotated_eigvals,1),[],2))])
    disp(['Homogeneity is better than ' num2str(nnz(mean(rotated_eigvals,1) < mean(realeigvals_per_first))) ' of ' num2str(size(rotated_eigvals,2)) ' random parcellations'])
    

    [parcellationsizes{parcelset} sorti] = sort(realsizes');
    
    all_sizes = [all_sizes realsizes];
    all_roteigvals = [all_roteigvals nanmean(rotated_eigvals,2)'];
    
    
    
    if plot_indfigs
        
        randadd = ((rand(size(realsizes))*2) -1) * .001;
        datain = [(realsizes + randadd)', realeigvals_per_first'];
        evalc('[dataout lowerLimit upperLimit xy] = lowess(datain,f,0);');
        realfittedvals{parcelset} = xy(:,2);
        sortedsizes{parcelset} = xy(:,1);
        
        datain = [(realsizes + randadd)', nanmean(rotated_eigvals,2)];
        datain(isnan(datain)) = nanmean(datain(:,2));
        evalc('[dataout lowerLimit upperLimit xy] = lowess(datain,f,0);');
        rotfittedvals{parcelset} = xy(:,2);
        
        
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
plotSpread(all_meanrotdata,'xNames',names,'distributionColors',repmat({'k'},1,length(names)),'MarkerSize',10)
hold on
plotSpread(all_meanrealdata,'xNames',names,'distributionColors',repmat({'r'},1,length(names)),'MarkerSize',50)

set(gcf,'Color',[1 1 1])
set(gca,'FontSize',15)



