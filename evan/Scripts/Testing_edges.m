outputfolder = '/data/cn4/evan/RestingState/FC_Mapping_120/Comparing_edges/';

hem = 'L';

numparcels_totest = [200 210 220 230 240 250];

edges_totest = {'/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/gradients_cifti_erode_wateredge/C1/avgcorrofcorr_allgrad_L_smooth2.55_wateredge_avg.func.gii',...
    '/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/gradients_cifti_erode_wateredge/C1/avgcorrofcorr_allgrad_L_smooth2.55_edge_avg.func.gii',...
    '/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/gradients_cifti_erode_concat_wateredge/C1/avgcorrofcorr_allgrad_L_smooth2.55_wateredge_avg.func.gii',...
    '/data/cn4/laumannt/fcMapping_redux/oldproc_wateredge_test/C1/avgcorrofcorr_smooth2.55_allgrad_L_smooth2.55_wateredge_avg.func.gii',...
    '/data/cn4/laumannt/fcMapping_redux/oldproc_wateredge_test/C1/avgcorrofcorr_smooth2.55_allgrad_L_smooth2.55_edge_avg.func.gii',...
    '/data/cn4/evan/RestingState/FC_Mapping_120/Comparing_edges/Newproc_smallwall_avg_waterthresh_fixed.func.gii',...
    '/data/hcp-bluearc/home/laumannt/120_parcellation/gradients_cifti_smallwall_wateredge/C1/avgcorrofcorr_allgrad_L_smooth2.55_wateredge_avg.func.gii'};

shortnames = {'Newproc_avg_water','Newproc_avg','Newproc_concat_water','Oldproc_avg_water','Oldproc_avg','Newproc_smallwall_avg_waterthresh','Newproc_smallwall_avg_water'};



for numparcels = numparcels_totest
    
    for edgemapnum = 1:length(edges_totest)
        edgemapname = edges_totest{edgemapnum};
        
        if ~exist([outputfolder '/' shortnames{edgemapnum} '_' num2str(numparcels) '_watershedmerge_noratio.func.gii'])
        disp([shortnames{edgemapnum} ' ' num2str(numparcels)])
        watershed_algorithm_merge_andthresh_targetnumber_noratio(edgemapname,outputfolder,[shortnames{edgemapnum} '_' num2str(numparcels) '_'],hem,numparcels)
        end
        
    end
end
%%
% avgcorr = gifti('/data/cn4/evan/RestingState/FC_Mapping_120/Subavg_corr_L.func.gii'); avgcorr = avgcorr.cdata;
% avgcorr(isnan(avgcorr)) = 0;
% load(['/data/cn4/evan/fsaverage_LR32k/Surface_distances_' hem '.mat'])

%%
for numparcels = numparcels_totest
    
    for edgemapnum = 1:length(edges_totest)
        %edgemapname = edges_totest{edgemapnum};
        %if ~exist([outputfolder '/' shortnames{edgemapnum} '_' num2str(numparcels) '.mat'])
            disp([shortnames{edgemapnum} ' ' num2str(numparcels)])
            generate_rotated_parcels_andPCA6([outputfolder '/' shortnames{edgemapnum} '_' num2str(numparcels) '_watershedmerge_noratio.func.gii'],100,cov_corr,2,hem,geo_distances,[outputfolder '/' shortnames{edgemapnum} '_' num2str(numparcels) '_corrected'])
        %end
    end
end

%%
colors = {'r','b','g','y','c','k','m'};
sizethresh = 0;
tstring = 'plot(';
meanstring = 'plot(';
tofZstring = 'plot(';
meanZstring = 'plot(';
tvals = zeros(length(edges_totest),length(numparcels_totest));
meandiffs = zeros(length(edges_totest),length(numparcels_totest));
tofZs = zeros(length(edges_totest),length(numparcels_totest));
meanZs = zeros(length(edges_totest),length(numparcels_totest));
for edgemapnum = 1:length(edges_totest)
    
    for  iparcelcount = 1:length(numparcels_totest)
    numparcels = numparcels_totest(iparcelcount);
    
        %edgemapname = edges_totest{edgemapnum};
        if exist([outputfolder '/' shortnames{edgemapnum} '_' num2str(numparcels) '_corrected.mat'])
            load([outputfolder '/' shortnames{edgemapnum} '_' num2str(numparcels) '_corrected.mat'])
            
            bigparcels = realsizes>sizethresh;
            
            [H,P,CI,STATS] = ttest(realeigvals_per_first,mean(rotated_eigvals,2)');
            tvals(edgemapnum,iparcelcount) = STATS.tstat;
            meandiffs(edgemapnum,iparcelcount) = mean(realeigvals_per_first - mean(rotated_eigvals,2)');
            
            Zs = zeros(length(realeigvals_per_first),1);
            for i = 1:length(realeigvals_per_first)
                Zs(i) = (realeigvals_per_first(i) - mean(rotated_eigvals(i,:)))/std(rotated_eigvals(i,:));
            end
            [H,P,CI,STATS] = ttest(Zs);
            tofZs(edgemapnum,iparcelcount) = STATS.tstat; 
            meanZs(edgemapnum,iparcelcount) = mean(Zs);
            
        end
    end
    
    tstring = [tstring 'numparcels_totest,tvals(' num2str(edgemapnum) ',:),''' colors{edgemapnum} '-'','];
    meanstring = [meanstring 'numparcels_totest,meandiffs(' num2str(edgemapnum) ',:),''' colors{edgemapnum} '-'','];
    tofZstring = [tofZstring 'numparcels_totest,tofZs(' num2str(edgemapnum) ',:),''' colors{edgemapnum} '-'','];
    meanZstring = [meanZstring 'numparcels_totest,meanZs(' num2str(edgemapnum) ',:),''' colors{edgemapnum} '-'','];
    nound_shortnames{edgemapnum} = shortnames{edgemapnum};
    nound_shortnames{edgemapnum}(nound_shortnames{edgemapnum}=='_') = [];
end
figure
eval([tstring(1:end-1) ')'])
legend(nound_shortnames)
title('T scores vs random')
figure
eval([meanstring(1:end-1) ')'])
legend(nound_shortnames)
title('Mean differences from random')
figure
eval([meanZstring(1:end-1) ')'])
legend(nound_shortnames)
title('Mean of Z scores vs random')
figure
eval([tofZstring(1:end-1) ')'])
legend(nound_shortnames)
title('1-sample T scores of Z scores vs random')

%%

load Newproc_smallwall_avg_waterthresh_230_corrected.mat
rep_realsizes = repmat(realsizes,100,1)';
plot(rep_realsizes(:),rotated_eigvals(:),'b.',realsizes,mean(rotated_eigvals,2)','g.',realsizes,realeigvals_per_first,'r.')
legend('Rotated parcels','Mean of rotated parcels','Real parcels')