threshs = [1.7 : .1 : 3.2];


for threshnum = 1:length(threshs)
    thresh = threshs(threshnum);
    load(['Poldrome_wateredge_L_' num2str(thresh) '_watershedmerge_randommatched2.mat'])
    
    realskew(threshnum) = skewness(realsizes);
    realkurt(threshnum) = kurtosis(realsizes);
    randomskew(threshnum) = skewness(reshape(randomsizes,numel(randomsizes),1));
    randomkurt(threshnum) = kurtosis(reshape(randomsizes,numel(randomsizes),1));
    
    realeigval(threshnum) = realmeaneigval;
    randomeigval(threshnum) = mean(randommeaneigvals);
    realsim(threshnum) = realmeansimval;
    randomsim(threshnum) = mean(randommeansimvals);
    
    realminrandomeigval(threshnum) = realeigval(threshnum) - randomeigval(threshnum);
    Zrealminrandomeigval(threshnum) = (realeigval(threshnum) - randomeigval(threshnum)) / std(randommeaneigvals);
    realminrandomsim(threshnum) = realsim(threshnum) - randomsim(threshnum);
    Zrealminrandomsim(threshnum) = (realsim(threshnum) - randomsim(threshnum)) / std(randommeansimvals);
    
    P = polyfit(realsizes',realeigvals_per_first',1);
    slope_realeigval(threshnum) = P(1); int_realeigval(threshnum) = (P(1) * 10) + P(2); 
    P = polyfit(realsizes',realsimilarities',1);
    slope_realsim(threshnum) = P(1); int_realsim(threshnum) = (P(1) * 10) + P(2); 
    
    P = polyfit(reshape(randomsizes,numel(randomsizes),1),reshape(randomeigvals_per_first,numel(randomeigvals_per_first),1),1);
    slope_randomeigval(threshnum) = P(1); int_randomeigval(threshnum) = (P(1) * 10) + P(2); 
    P = polyfit(reshape(randomsizes,numel(randomsizes),1),reshape(randomsimilarities,numel(randomeigvals_per_first),1),1);
    slope_randomsim(threshnum) = P(1); int_randomsim(threshnum) = (P(1) * 10) + P(2); 
    
    
%     [B_realeigval(threshnum),Bint_realeigval(threshnum),resid_realeigvals] = regress(realeigvals_per_first',realsizes',ones(numel(realsizes)));
%     resid_realeigval(threshnum) = mean(resid_allrealeigvals);
%     [B_realsim(threshnum),Bint_realsim(threshnum),resid_realsims] = regress(realsimilarities',realsizes',ones(numel(realsizes)));
%     resid_realsim(threshnum) = mean(resid_allrealsims);
    
%     [B_randomeigval(threshnum),Bint_randomeigval(threshnum),resid_randomeigval] = regress(reshape(randomeigvals_per_first,numel(randomeigvals_per_first),1),reshape(randomsizes,numel(randomsizes),1),ones(numel(randomsizes),1));
%     [B_randomsim(threshnum),Bint_randomsim(threshnum),resid_randomsim] = regress(reshape(randomeigvals_per_first,numel(randomeigvals_per_first),1),reshape(randomsizes,numel(randomsizes),1),ones(numel(randomsizes),1));
    
    
    slope_realminrandomeigval(threshnum) = slope_realeigval(threshnum) - slope_randomeigval(threshnum);
    slope_realminrandomsim(threshnum) = slope_realsim(threshnum) - slope_randomsim(threshnum);
    
    int_realminrandomeigval(threshnum) = int_realeigval(threshnum) - int_randomeigval(threshnum);
    int_realminrandomsim(threshnum) = int_realsim(threshnum) - int_randomsim(threshnum);
    
    %resid_realminrandomeigval(threshnum) = mean(resid_realeigval) - mean(resid_randomeigval);
    %resid_realminrandomsim(threshnum) = mean(resid_realsim) - mean(resid_randomsim);
    

    
end


    

% plot(threshs,realeigval);
% title('Real and Random Homogeneities');
% set(gca,'XTick',threshs);
% hold
% plot(threshs,randomeigval,'r');
% set(gcf,'Nextplot','new');
% 
% plot(threshs,realsim);
% title('Real and Random Similarities');
% set(gca,'XTick',threshs);
% hold
% plot(threshs,randomsim,'r');
% set(gcf,'Nextplot','new');
% 
% plot(threshs,Zrealminrandomeigval);
% title('Real - Random Homogeneities');
% set(gca,'XTick',threshs);
% set(gcf,'Nextplot','new');
% 
% plot(threshs,Zrealminrandomsim);
% title('Real - Random Similarities');
% set(gca,'XTick',threshs);
% set(gcf,'Nextplot','new');



plot(threshs,realskew);
title('Real and Random Size Skewness');
set(gca,'XTick',threshs);
hold
plot(threshs,randomskew,'r');
set(gcf,'Nextplot','new');

plot(threshs,realkurt);
title('Real and Random Size Kurtosis');
set(gca,'XTick',threshs);
hold
plot(threshs,randomkurt,'r');
set(gcf,'Nextplot','new');



% plot(threshs,slope_realeigval);
% title('Real and Random Homogeneity size dependence');
% set(gca,'XTick',threshs);
% hold
% plot(threshs,slope_randomeigval,'r');
% set(gcf,'Nextplot','new');
% 
% plot(threshs,slope_realsim);
% title('Real and Random Residual Similarity size dependence');
% set(gca,'XTick',threshs);
% hold
% plot(threshs,slope_randomsim,'r');
% set(gcf,'Nextplot','new');
% 
% plot(threshs,slope_realminrandomeigval);
% title('Real - Random Homogeneity size dependence');
% set(gca,'XTick',threshs);
% set(gcf,'Nextplot','new');
% 
% plot(threshs,slope_realminrandomsim);
% title('Real - Random Similarity size dependence');
% set(gca,'XTick',threshs);
% set(gcf,'Nextplot','new');





% plot(threshs,int_realeigval);
% title('Real and Random Homogeneity v Size intercept');
% set(gca,'XTick',threshs);
% hold
% plot(threshs,int_randomeigval,'r');
% set(gcf,'Nextplot','new');
% 
% plot(threshs,int_realsim);
% title('Real and Random Similarity v Size intercept');
% set(gca,'XTick',threshs);
% hold
% plot(threshs,int_randomsim,'r');
% set(gcf,'Nextplot','new');
% 
% plot(threshs,int_realminrandomeigval);
% title('Real - Random Homogeneity v Size intercept');
% set(gca,'XTick',threshs);
% set(gcf,'Nextplot','new');
% 
% plot(threshs,int_realminrandomsim);
% title('Real - Random Similarity v Size intercept');
% set(gca,'XTick',threshs);
% set(gcf,'Nextplot','new');




% plot(threshs,resid_realeigval);
% title('Real and Random Residual Homogeneities');
% set(gca,'XTick',threshs);
% hold
% plot(threshs,resid_randomeigval,'r');
% set(gcf,'Nextplot','new');
% 
% plot(threshs,resid_realsim);
% title('Real and Random Residual Similarities');
% set(gca,'XTick',threshs);
% hold
% plot(threshs,resid_randomsim,'r');
% set(gcf,'Nextplot','new');
% 
% plot(threshs,resid_realminrandomeigval);
% title('Real - Random Residual Homogeneities');
% set(gca,'XTick',threshs);
% set(gcf,'Nextplot','new');
% 
% plot(threshs,resid_realminrandomsim);
% title('Real - Random Residual Similarities');
% set(gca,'XTick',threshs);
% set(gcf,'Nextplot','new');