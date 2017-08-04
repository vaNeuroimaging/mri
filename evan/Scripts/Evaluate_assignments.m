%% Assignment characteristics

group_parcels = cifti_read('/data/cn4/evan/Published_parcels/Parcels_LR.dtseries.nii');
group_parcelIDs = unique(group_parcels); group_parcelIDs(group_parcelIDs==0) = [];

assigned_parcels = cifti_read('Assignments_xd15_1.0thr_merged.dtseries.nii');
assigned_parcelIDs = unique(assigned_parcels); assigned_parcelIDs(assigned_parcelIDs<1) = [];

% group_networks = cifti_read('/data/cn4/evan/RestingState/Consensus/120_LR_minsize400_recolored_manualconsensus_LR_cleaned.dtseries.nii');
% networkIDs = unique(group_networks); networkIDs(networkIDs<1) = [];
% groupnetworks_noassign = find(group_networks < 1);
% 
% sub_networks = cifti_read('Template_match_systems_Poldrome_MSC.dtseries.nii');
% sub_networks(groupnetworks_noassign,:) = 0;



orig_parcels = cifti_read('Parcels_Poldrome_MSC.dtseries.nii');

for s = 1:size(assigned_parcels,2)
    parcels_unassignedmask = single(logical(orig_parcels(:,s)>0)) - single(logical(assigned_parcels(:,s)>0));
    parcels_unassignedmask(parcels_unassignedmask<1) = 0;
    numparcels_unassigned = nnz(unique(orig_parcels(:,s) .* parcels_unassignedmask)) - 1;
    disp(['Subject ' num2str(s) ': ' num2str(numparcels_unassigned) ' parcels unassigned out of ' num2str((nnz(unique(orig_parcels(:,s)))-1)) ' total parcels, covering ' num2str(nnz(parcels_unassignedmask) / numel(parcels_unassignedmask) * 100) '% of vertices'])
end
    

MSCnames = {'MSC01','MSC02','MSC03','MSC04'};
%% Resting State Similarity Calculations
group_parcelmat = zeros(length(group_parcelIDs), length(group_parcelIDs), size(assigned_parcels,2));
assigned_parcelmat = zeros(length(assigned_parcelIDs), length(assigned_parcelIDs), size(assigned_parcels,2));
assigned_parcels_presentvecs = zeros(length(assigned_parcelIDs), size(assigned_parcels,2));

%group_networkmat = zeros(length(networkIDs),length(networkIDs), size(assigned_parcels,2));
%sub_networkmat = zeros(length(networkIDs),length(networkIDs), size(assigned_parcels,2));

for s = 1:size(assigned_parcels,2)
    disp(['Subject ' num2str(s)])
    if s==1
        subtimeseries = cifti_read('/data/cn4/evan/RestingState/Ind_variability/Poldrome/allsubs_LR_timeseries.dtseries.nii');
        tmask = load('/data/cn4/evan/RestingState/Ind_variability/Poldrome/allsubs_total_tmask.txt');
        
    else
        MSCname = MSCnames{s-1};
        tmaskfile = ['/data/nil-bluearc/GMT/Laumann/MSC/' MSCname '/' MSCname '_TMASKLIST.txt'];
        [subjects tmasks] = textread(tmaskfile,'%s %s');
        for i = 1:length(subjects)
            if i == 1
                subtimeseries = cifti_read(['/data/nil-bluearc/GMT/Laumann/MSC/' MSCname '/Functionals/FCPROCESS_SCRUBBED_UWRPMEAN/cifti_timeseries_normalwall/' subjects{i} '_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii']);
                tmask = load(tmasks{i});
            else
                subtimeseries = [subtimeseries cifti_read(['/data/nil-bluearc/GMT/Laumann/MSC/' MSCname '/Functionals/FCPROCESS_SCRUBBED_UWRPMEAN/cifti_timeseries_normalwall/' subjects{i} '_LR_surf_subcort_333_32k_fsLR_smooth2.55.dtseries.nii'])];
                tmask = [tmask; load(tmasks{i})];
            end
        end
        
    end
    subtimeseries = subtimeseries(:,logical(tmask));
    
    
    parceltimeseries = zeros(size(subtimeseries,2),length(assigned_parcelIDs));
    for parcelnum = 1:length(assigned_parcelIDs)
        if any(assigned_parcels(:,s)==assigned_parcelIDs(parcelnum));
            
            assigned_parcels_presentvecs(parcelnum,s) = 1;
            
            parceltimeseries(:,parcelnum) = mean(subtimeseries(logical(assigned_parcels(:,s)==assigned_parcelIDs(parcelnum)),:),1);
        end
    end
    assigned_parcelmat(:,:,s) = FisherTransform(paircorr_mod(parceltimeseries));
    
    
    parceltimeseries = zeros(size(subtimeseries,2),length(group_parcelIDs));
    for parcelnum = 1:length(group_parcelIDs)
        parceltimeseries(:,parcelnum) = mean(subtimeseries(logical(group_parcels==group_parcelIDs(parcelnum)),:),1);
    end
    group_parcelmat(:,:,s) = FisherTransform(paircorr_mod(parceltimeseries));
    
%     group_networktimeseries = zeros(size(subtimeseries,2),length(networkIDs));
%     sub_networktimeseries = zeros(size(subtimeseries,2),length(networkIDs));
%     for networknum = 1:length(networkIDs)
%         group_networktimeseries(:,networknum) = mean(subtimeseries(logical(group_networks==networkIDs(networknum)),:),1);
%         sub_networktimeseries(:,networknum) = mean(subtimeseries(logical(sub_networks(:,s)==networkIDs(networknum)),:),1);
%     end
%     group_networkmat(:,:,s) = paircorr_mod(group_networktimeseries);
%     sub_networkmat(:,:,s) = paircorr_mod(sub_networktimeseries);
    
end

%% Resting State Similarity Results

grouptogroup = zeros(size(assigned_parcels,2));
assignedtoassigned = zeros(size(assigned_parcels,2));
%groupnettogroupnet = zeros(size(assigned_parcels,2));
%subnettosubnet = zeros(size(assigned_parcels,2));

for s1 = 1:size(assigned_parcels,2)
    for s2 = 1:size(assigned_parcels,2)
        
        if s1 < s2
            
            gs1mat = group_parcelmat(:,:,s1);
            gs1vals = gs1mat(logical(triu(ones(size(group_parcelmat(:,:,1))),1)));
            gs2mat = group_parcelmat(:,:,s2);
            gs2vals = gs2mat(logical(triu(ones(size(group_parcelmat(:,:,1))),1)));
            
            grouptogroup(s1,s2) = FisherTransform(paircorr_mod(gs1vals(:),gs2vals(:)));
            
            
            
            commonparcelinds = logical(assigned_parcels_presentvecs(:,s1) .* assigned_parcels_presentvecs(:,s2));
            %disp(['Subject ' num2str(s1) ' vs Subject ' num2str(s2) ': ' num2str(nnz(commonparcelinds)) ' parcels matched'])
            commons1mat = assigned_parcelmat(commonparcelinds,commonparcelinds,s1);
            commons2mat = assigned_parcelmat(commonparcelinds,commonparcelinds,s2);
            as1vals = commons1mat(logical(triu(ones(size(commons1mat)),1)));
            as2vals = commons2mat(logical(triu(ones(size(commons2mat)),1)));
            
            assignedtoassigned(s1,s2) = FisherTransform(paircorr_mod(as1vals(:),as2vals(:)));
            
            
%             gn1mat = group_networkmat(:,:,s1);
%             gn1vals = gn1mat(logical(triu(ones(size(group_networkmat(:,:,1))),1)));
%             gn2mat = group_networkmat(:,:,s2);
%             gn2vals = gn2mat(logical(triu(ones(size(group_networkmat(:,:,1))),1)));
%             
%             groupnettogroupnet(s1,s2) = paircorr_mod(gn1vals(:),gn2vals(:));
%             
%             
%             sn1mat = sub_networkmat(:,:,s1);
%             sn1vals = sn1mat(logical(triu(ones(size(sub_networkmat(:,:,1))),1)));
%             sn2mat = sub_networkmat(:,:,s2);
%             sn2vals = sn2mat(logical(triu(ones(size(sub_networkmat(:,:,1))),1)));
%             
%             subnettosubnet(s1,s2) = paircorr_mod(sn1vals(:),sn2vals(:));
            
            
            
            
        end
    end
end

triumat = logical(triu(ones(size(grouptogroup)),1));

disp(['Mean group parcel matrix REST similarities: ' num2str(mean(grouptogroup(triumat))) ' +/- ' num2str(std(grouptogroup(triumat)))])
disp(['Mean assigned parcel matrix REST similarities: ' num2str(mean(assignedtoassigned(triumat))) ' +/- ' num2str(std(assignedtoassigned(triumat)))])
[H,P,CI,STATS] = ttest(assignedtoassigned(triumat),grouptogroup(triumat));
disp(['Paired t-test REST results: T=' num2str(STATS.tstat) ', p=' num2str(P)])

% disp(['Mean group network matrix REST similarities: ' num2str(mean(groupnettogroupnet(triumat))) ' +/- ' num2str(std(groupnettogroupnet(triumat)))])
% disp(['Mean subject network matrix REST similarities: ' num2str(mean(subnettosubnet(triumat))) ' +/- ' num2str(std(subnettosubnet(triumat)))])
% [H,P,CI,STATS] = ttest(subnettosubnet(triumat),groupnettogroupnet(triumat));
% disp(['Paired t-test REST results: T=' num2str(STATS.tstat) ', p=' num2str(P)])


figure;imagesc(grouptogroup,[.4 .75]); colormap hot; title('Subject REST similarities: group parcels'); colorbar
figure;imagesc(assignedtoassigned,[.4 .75]); colormap hot; title('Subject REST similarities: assigned parcels'); colorbar
%figure;imagesc(groupnettogroupnet,[.4 .75]); colormap hot; title('Subject REST similarities: group networks'); colorbar
%figure;imagesc(subnettosubnet,[.4 .75]); colormap hot; title('Subject REST similarities: subject networks'); colorbar

figure;imagesc(assignedtoassigned - grouptogroup,[-.2 .2]); 
colormap hot; hotmap = colormap; coolmap = [hotmap(:,3),hotmap(:,2),hotmap(:,1)];
hotmap = hotmap(1:end-3, :); coolmap = coolmap(1:end-3, :);
combined = [flipdim(coolmap,1); zeros(10,3); hotmap];
colormap(combined);
title('Subject REST similarities: assigned parcels > group parcels'); colorbar

%% Task Similarity Calculations




task_maps = {[1 2 5 7 9 11 13 15 17],[1],[1:2:23]};

assigned_activations = ones(length(assigned_parcelIDs),length(cell2mat(task_maps)),size(assigned_parcels,2)) * NaN;
group_activations = ones(length(group_parcelIDs),length(cell2mat(task_maps)),size(assigned_parcels,2)) * NaN;
assigned_parcels_presentvecs = zeros(length(assigned_parcelIDs), size(assigned_parcels,2));

%groupnetwork_activations = ones(length(networkIDs),length(cell2mat(task_maps)),size(assigned_parcels,2)) * NaN;
%subnetwork_activations = ones(length(networkIDs),length(cell2mat(task_maps)),size(assigned_parcels,2)) * NaN;

for s = 2:size(assigned_parcels,2)
    MSCname = MSCnames{s-1};
    task_files = {['/data/nil-bluearc/GMT/MSC/' MSCname '/master_spec/' MSCname '.motor_taskmaps.dscalar.nii'],['/data/nil-bluearc/GMT/MSC/' MSCname '/master_spec/' MSCname '.memory_taskmaps.dscalar.nii'],['/data/nil-bluearc/GMT/MSC/' MSCname '/master_spec/' MSCname '.mixed_taskmaps.dscalar.nii']};
    contrastcounter = 0;
    for filenum = 1:length(task_files)
        data = cifti_read(task_files{filenum});
        for contrastnum = 1:length(task_maps{filenum})
            contrast = task_maps{filenum}(contrastnum);
            contrastcounter = contrastcounter +1;
            
            for parcelIDnum = 1:length(assigned_parcelIDs)
                parcelID = assigned_parcelIDs(parcelIDnum);
                
                if any(assigned_parcels(:,s)==assigned_parcelIDs(parcelIDnum));
                    assigned_parcels_presentvecs(parcelIDnum,s) = 1;
                    assigned_activations(parcelIDnum,contrastcounter,s) = mean(data(assigned_parcels(:,s)==parcelID,contrast));
                end
            end
            
            for parcelIDnum = 1:length(group_parcelIDs)
                parcelID = group_parcelIDs(parcelIDnum);
                group_activations(parcelIDnum,contrastcounter,s) = mean(data(group_parcels==parcelID,contrast));
            end
            
%             for networknum = 1:length(networkIDs)
%                 networkID = networkIDs(networknum);
%                 groupnetwork_activations(networknum,contrastcounter,s) = mean(data(group_networks==networkID,contrast));
%                 subnetwork_activations(networknum,contrastcounter,s) = mean(data(sub_networks(:,s)==networkID,contrast));
%             end
            
        end
    end
end


%% Task Similarity Results

grouptogroup = zeros(size(assigned_parcels,2));
assignedtoassigned = zeros(size(assigned_parcels,2));
%groupnettogroupnet = zeros(size(assigned_parcels,2));
%subnettosubnet = zeros(size(assigned_parcels,2));

for s1 = 2:size(assigned_parcels,2)
    for s2 = 2:size(assigned_parcels,2)
        
        if s1 < s2
            
            gs1mat = group_activations(:,:,s1);
            gs2mat = group_activations(:,:,s2);
            
            grouptogroup(s1,s2) = FisherTransform(paircorr_mod(gs1mat(:),gs2mat(:)));
            
            commonparcelinds = logical(assigned_parcels_presentvecs(:,s1) .* assigned_parcels_presentvecs(:,s2));
            as1mat = assigned_activations(commonparcelinds,:,s1);
            as2mat = assigned_activations(commonparcelinds,:,s2);
            
            assignedtoassigned(s1,s2) = FisherTransform(paircorr_mod(as1mat(:),as2mat(:)));
            
%             gn1mat = groupnetwork_activations(:,:,s1);
%             gn2mat = groupnetwork_activations(:,:,s2);
%             
%             groupnettogroupnet(s1,s2) = paircorr_mod(gn1mat(:),gn2mat(:));
%             
%             sn1mat = subnetwork_activations(:,:,s1);
%             sn2mat = subnetwork_activations(:,:,s2);
%             
%             subnettosubnet(s1,s2) = paircorr_mod(sn1mat(:),sn2mat(:));
            
        end
    end
end

triumat = logical(triu(ones(size(grouptogroup)),1));
triumat(1,:) = 0;

disp(['Mean group parcel matrix TASK similarities: ' num2str(mean(grouptogroup(triumat))) ' +/- ' num2str(std(grouptogroup(triumat)))])
disp(['Mean assigned parcel matrix TASK similarities: ' num2str(mean(assignedtoassigned(triumat))) ' +/- ' num2str(std(assignedtoassigned(triumat)))])
[H,P,CI,STATS] = ttest(assignedtoassigned(triumat),grouptogroup(triumat));
disp(['Paired t-test results for TASK: T=' num2str(STATS.tstat) ', p=' num2str(P)])

% disp(['Mean group network matrix TASK similarities: ' num2str(mean(groupnettogroupnet(triumat))) ' +/- ' num2str(std(groupnettogroupnet(triumat)))])
% disp(['Mean subject network matrix TASK similarities: ' num2str(mean(subnettosubnet(triumat))) ' +/- ' num2str(std(subnettosubnet(triumat)))])
% [H,P,CI,STATS] = ttest(subnettosubnet(triumat),groupnettogroupnet(triumat));
% disp(['Paired t-test TASK results: T=' num2str(STATS.tstat) ', p=' num2str(P)])


figure;imagesc(grouptogroup,[.4 .75]); colormap hot; title('Subject REST similarities: group parcels'); colorbar
figure;imagesc(assignedtoassigned,[.4 .75]); colormap hot; title('Subject REST similarities: assigned parcels'); colorbar
%figure;imagesc(groupnettogroupnet,[.4 .75]); colormap hot; title('Subject REST similarities: group networks'); colorbar
%figure;imagesc(subnettosubnet,[.4 .75]); colormap hot; title('Subject REST similarities: subject networks'); colorbar

figure;imagesc(assignedtoassigned - grouptogroup,[-.2 .2]); 
colormap hot; hotmap = colormap; coolmap = [hotmap(:,3),hotmap(:,2),hotmap(:,1)];
hotmap = hotmap(1:end-3, :); coolmap = coolmap(1:end-3, :);
combined = [flipdim(coolmap,1); zeros(10,3); hotmap];
colormap(combined);
title('Subject TASK similarities: assigned parcels > group parcels'); colorbar