function how_much_data_summedweightedPC_MSC

run_new = 1;

PC_thresh = .75;
thresholdarray = [.01 : .01 : .2];
MSCnums = 1:10;

outfolder = ['/data/nil-bluearc/GMT/Evan/MSC/Hub_detection/Session_Homogenous_Infomap_summedPC_convergence'];
mkdir(outfolder);

xdistance = 30;
networksizeminimum = 5;
iterations = 50;
sessionnums_totest = [1/3 2/3 1:5];
sesslength = 30;
dataquant = sessionnums_totest * sesslength;

%if run_new; pool = parpool(10); end

for MSCnum = MSCnums
    MSCname = ['MSC' sprintf('%02i',MSCnum)];
    
    if run_new
    
    ciftiparcelsdir = ['/data/nil-bluearc/GMT/Evan/MSC/subjects/' MSCname '/parcellation/'];
    parcels_LR = [ciftiparcelsdir '/RSFC_parcels_edgethresh_0.5_homogenous.dtseries.nii'];
    parcels_struct = ft_read_cifti_mod(parcels_LR);
    parcels = parcels_struct.data;
    parcelIDs = unique(parcels); parcelIDs(parcelIDs<1) = [];
    
    
    %Calculate parcel distances
    distancesfile = ['/data/nil-bluearc/GMT/Laumann/MSC/' MSCname '/normalwall_distmat_333_native_freesurf/distmat_surf_geodesic_vol_euc_uint8.mat'];
    distances = smartload(distancesfile);
    parcel_centers = zeros(length(parcelIDs),1);
    for IDnum = 1:length(parcelIDs)
        parcelinds = find(parcels==parcelIDs(IDnum));
        distancesum = sum(distances(parcelinds,parcelinds),2);
        [~,mini] = min(distancesum);
        parcel_centers(IDnum) = parcelinds(mini);
    end
    parcel_distances = distances(parcel_centers,parcel_centers);
    save([ciftiparcelsdir '/parcel_distances.mat'],'parcel_distances')
    clear distances
    
    
    
    parcel_distances = smartload([ciftiparcelsdir '/parcel_distances.mat']);
    parcel_distances_tooclose = (parcel_distances < xdistance);
    
    
    
    
    
    
%     tmaskfile = ['/data/nil-bluearc/GMT/Laumann/MSC/' MSCname '/' MSCname '_TMASKLIST.txt'];
%     [subjectlist, tmask_list] = textread(tmaskfile,'%s %s');
%     %subjectlist = subjectlist(2:end); tmask_list = tmask_list(2:end);
%     cifti_all = [];
%     tmask_all = [];
%     sessnum = [];
%     withinsess_ind = [];
%     for s = 1 : length(subjectlist)
%         ciftifiles{s} = ['/data/nil-bluearc/GMT/Laumann/MSC/' MSCname '/Functionals/FCPROCESS_SCRUBBED_UWRPMEAN/cifti_timeseries_normalwall_native_freesurf/' subjectlist{s} '_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii'];
%         tmask = load(tmask_list{s});
%         tmask_all = [tmask_all ; tmask];
%         sessnum = [sessnum ; (ones(size(tmask,1),1) * s)];
%         withinsess_ind = [withinsess_ind ; [1:length(tmask)]'];
%         temp_data = ft_read_cifti_mod(ciftifiles{s});
%         cifti_all = [cifti_all temp_data.data];
%     end
%     tmask_all = logical(tmask_all);
%     cifti_all = cifti_all(:,tmask_all);
%     sessnum = sessnum(tmask_all);
%     withinsess_ind = withinsess_ind(tmask_all);
    
    
  
    
%    numsessions = length(subjectlist);
    
    end
    
    corrmat_similarity = zeros(iterations,length(sessionnums_totest));
    PC_similarity = zeros(iterations,length(sessionnums_totest));
    Hub_similarity = zeros(iterations,length(sessionnums_totest));
    
    prevstring = [];
    
    for amountnum = 1:length(sessionnums_totest)
        for iter = 1:iterations
        
            half_outfolder = [outfolder '/' MSCname '_' num2str(dataquant(amountnum)) '_iter' num2str(iter) '_half/'];
            this_outfolder = [outfolder '/' MSCname '_' num2str(dataquant(amountnum)) '_iter' num2str(iter) '/'];
            
            if run_new
                disp([MSCname ': iteration ' num2str(iter) ', ' num2str(dataquant(amountnum)) ' minutes of data'])
                mkdir(this_outfolder);
                cd(this_outfolder)
                
%                 randorder = randperm(numsessions);
%                 
%                 thishalf = randorder(1:floor(numsessions/2));
%                 otherhalf = setdiff(randorder,thishalf);
%                 
%                 
%                 amount = sessionnums_totest(amountnum);
%                 if amount < 1;
%                     chosenrun = thishalf(1);
%                     inds_touse = find(sessnum==chosenrun);
%                     dataamount = floor(sessionnums_totest(amountnum) * length(inds_touse));
%                     inds_touse = inds_touse(1:dataamount);
%                 else
%                     chosenruns = thishalf(1:sessionnums_totest(amountnum));
%                     inds_touse = [];
%                     for r = 1:length(chosenruns)
%                         inds_touse = [inds_touse; find(sessnum==chosenruns(r))];
%                     end
%                 end
%                 
%                 thisiter_cifti = cifti_all(:,inds_touse);
%                 temp_data.data = thisiter_cifti;
%                 parcel_infomap_singlesub_knowndistances(parcels_struct,temp_data,[ciftiparcelsdir '/parcel_distances.mat'],this_outfolder,thresholdarray,xdistance,networksizeminimum,1)
                
                
                
                all_parcel_corrmat = smartload([this_outfolder '/corrmat.mat']);
                communities = load([this_outfolder '/rawassn_minsize' num2str(networksizeminimum) '.txt']);
                mean_PCs = PC_calc(all_parcel_corrmat,parcel_distances_tooclose,communities,thresholdarray);
                save([this_outfolder '/mean_wPCs.mat'],'mean_PCs')
                
                sorted = sort(mean_PCs,'ascend');
                PC_threshold_number = sorted(round(length(sorted) * PC_thresh));
                Hubs = (mean_PCs > PC_threshold_number);
                save([this_outfolder '/Hubs_wPCthr' num2str(PC_thresh) '.mat'],'Hubs')
                
                
%                 mkdir(half_outfolder);
%                 cd(half_outfolder)
%                 disp([MSCname ': iteration ' num2str(iter) ', ' num2str(dataquant(amountnum)) ' minutes of data, random half'])
%                 
%                 
%                 
%                 inds_touse = [];
%                 for r = 1:length(otherhalf)
%                     inds_touse = [inds_touse; find(sessnum==otherhalf(r))];
%                 end
%                 
%                 otherhalf_cifti = cifti_all(:,inds_touse);
%                 temp_data.data = otherhalf_cifti;
%                 parcel_infomap_singlesub_knowndistances(parcels_struct,temp_data,[ciftiparcelsdir '/parcel_distances.mat'],half_outfolder,thresholdarray,xdistance,networksizeminimum,1)
                
                
                all_parcel_corrmat = smartload([half_outfolder '/corrmat.mat']);
                communities = load([half_outfolder '/rawassn_minsize' num2str(networksizeminimum) '.txt']);
                mean_PCs = PC_calc(all_parcel_corrmat,parcel_distances_tooclose,communities,thresholdarray);
                mean_PCs(isnan(mean_PCs)) = 0;
                save([half_outfolder '/mean_wPCs.mat'],'mean_PCs')
                
                sorted = sort(mean_PCs,'ascend');
                PC_threshold_number = sorted(round(length(sorted) * PC_thresh));
                Hubs = (mean_PCs > PC_threshold_number);
                save([half_outfolder '/Hubs_wPCthr' num2str(PC_thresh) '.mat'],'Hubs')
                
                
                
            else
                
                string = [MSCname ': iteration ' num2str(iter) ', ' num2str(dataquant(amountnum)) ' minutes of data'];
                fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
                prevstring = string;
            end
            
            
            
            %%
            half_corrmat = smartload([half_outfolder '/corrmat.mat']);
            half_PC = smartload([half_outfolder '/mean_wPCs.mat']);
            half_PC(isnan(half_PC)) = 0;
           if exist([half_outfolder '/Hubs_wPCthr' num2str(PC_thresh) '.mat'])
               half_Hubs = smartload([half_outfolder '/Hubs_wPCthr' num2str(PC_thresh) '.mat']);
           else
                sorted = sort(half_PC,'ascend');
                PC_threshold_number = sorted(round(length(sorted) * PC_thresh));
                Hubs = (half_PC > PC_threshold_number);
                save([half_outfolder '/Hubs_wPCthr' num2str(PC_thresh) '.mat'],'Hubs')
                half_Hubs = Hubs;
           end
            
            iter_corrmat = smartload([this_outfolder '/corrmat.mat']);
            iter_PC = smartload([this_outfolder '/mean_wPCs.mat']);
            iter_PC(isnan(iter_PC)) = 0;
          if exist([this_outfolder '/Hubs_wPCthr' num2str(PC_thresh) '.mat'])
              iter_Hubs = smartload([this_outfolder '/Hubs_wPCthr' num2str(PC_thresh) '.mat']);
          else
                sorted = sort(iter_PC,'ascend');
                PC_threshold_number = sorted(round(length(sorted) * PC_thresh));
                Hubs = (iter_PC > PC_threshold_number);
                save([this_outfolder '/Hubs_wPCthr' num2str(PC_thresh) '.mat'],'Hubs')
                iter_Hubs = Hubs;
          end
                
            
            corrmat_similarity(iter,amountnum) = paircorr_mod(half_corrmat(triu(true(size(half_corrmat)),1)),iter_corrmat(triu(true(size(half_corrmat)),1)));
            PC_similarity(iter,amountnum) = corr(half_PC,iter_PC,'type','Spearman');
            Hub_similarity(iter,amountnum) = (2*nnz(half_Hubs & iter_Hubs)) / (nnz(half_Hubs) + nnz(iter_Hubs));
                        
            
        end
    end
    disp(' ')
    
    
    h = figure('Color','white','position',[1982 478 1352 804]);
    hold on
    
    meancorrmat = mean(corrmat_similarity,1);
    stdcorrmat = std(corrmat_similarity,0,1);
    plot(dataquant,meancorrmat,'Color','k','LineWidth',3)
    %plot(dataquant,meancorrmat - stdcorrmat,'--k','LineWidth',2)
    %plot(dataquant,meancorrmat + stdcorrmat,'--k','LineWidth',2)
    
    
    meanPC = mean(PC_similarity,1);
    stdPC = std(PC_similarity,0,1);
    plot(dataquant,meanPC,'Color','b','LineWidth',3)
    %plot(dataquant,meanPC - stdPC,'--b','LineWidth',2)
    %plot(dataquant,meanPC + stdPC,'--b','LineWidth',2)
    
    meanHub = mean(Hub_similarity,1);
    stdHub = std(Hub_similarity,0,1);
    plot(dataquant,meanHub,'Color','r','LineWidth',3)
    %plot(dataquant,meanHub - stdHub,'--r','LineWidth',2)
    %plot(dataquant,meanHub + stdHub,'--r','LineWidth',2)
    
    %ylim([.65 .9])
    
    set(gca,'FontSize',14,'FontWeight','bold','LineWidth',2)
    title(gca,['Convergence Plots: ' MSCname])
    xlabel('Time (minutes)','Fontweight','bold','FontSize',14)
    ylabel('Similarity to other half','Fontweight','bold','FontSize',14)
    set(gca,'XTick',[0:sesslength:dataquant(end)])
    ylim([.3 1])
    
    export_fig(gca,[outfolder '/Convergence Plots: ' MSCname ', wPC' num2str(PC_thresh) '.pdf'])
    
end

end

function mean_PCs = PC_calc(all_parcel_corrmat,parcel_distances_tooclose,communities,thresholdarray)




all_parcel_corrmat(logical(parcel_distances_tooclose)) = 0;
parcel_correlation_vec = all_parcel_corrmat(triu(true(size(all_parcel_corrmat)),1));
sorted_vals = sort(parcel_correlation_vec,'descend');

all_PCs = zeros(size(all_parcel_corrmat,1),length(thresholdarray));

%PCvals = zeros(size(communities));

for threshnum = 1:length(thresholdarray)
    kdenthresh = thresholdarray(threshnum);
    rthresh = sorted_vals(ceil(length(sorted_vals) .* kdenthresh));
    corrmat_weighted = all_parcel_corrmat .* (all_parcel_corrmat > rthresh);
    
    thresh_communities = unique(communities(:,threshnum));
    thresh_communities(thresh_communities<1) = [];
    
    parcel_degree = sum(corrmat_weighted,2);
    moduleconnectionssum = zeros(size(corrmat_weighted,1),1);
    
    for communitynum = 1:length(thresh_communities)
        communityID = thresh_communities(communitynum);
        communityindices = communities(:,threshnum)==communityID;
        parcel_community_degree = sum(corrmat_weighted(:,communityindices),2);
        parcel_community_ratio = (parcel_community_degree ./ parcel_degree);
        moduleconnectionssum = moduleconnectionssum + (parcel_community_ratio.^2);
    end
    
    all_PCs(:,threshnum) = 1-moduleconnectionssum;
end

mean_PCs = nanmean(all_PCs,2);
mean_PCs(isnan(mean_PCs)) = 0;

end
