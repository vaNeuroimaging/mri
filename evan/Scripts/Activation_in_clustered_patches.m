%function Activation_in_clustered_patches(probmapfile)

probmapfile = 'Cluster_probability_maps_sorted_10mm_27sub_selected.dscalar.nii';
tmaskfile = '/data/cn4/evan/RestingState/Ind_variability/HCP/HCP_80_TMASKLIST.txt';
tasks = {'EMOTION','GAMBLING','LANGUAGE','MOTOR','RELATIONAL','SOCIAL','WM'};

ncortverts = 59412;

[subs ign] = textread(tmaskfile,'%s%s');
%%

probmaps_template = ft_read_cifti_mod(probmapfile);
probmap = probmaps_template.data;
probmaps_template.data = [];
prevstring = [];

% inds = false(size(probmap,2),1);
% for col = 1:size(probmap,2)
%     tokens = tokenize(probmaps_template.mapname{col},':');
%     networkname = tokens{1};
%     tokens = tokenize(probmaps_template.mapname{col},'=');
%     tokens2 = tokenize(tokens{2},';');
%     patchsize = str2num(tokens2{1}(1:end-3));
%     if patchsize > 250
%         inds(col) = true;
%     end
% end
% probmap = probmap(:,inds);

load System_clustering.mat
simple_assigns = load('rawassn_minsize2.txt');

communities = unique(simple_assigns); communities(communities<1) = [];

networklabels = {'Default','Visual','FrontoPar','FrontoPar 2','DorsalAttn','DorsalAttn2','VentAttn','Salience','CingOperc','MortorHand','MotorMouth','Auditory','MTL1','MTL2','MedPar','ParOccip'};

%rs = zeros(length(communities),1);

matchedpatches_bysub = false(size(probmap,1),length(subs),size(probmap,2));

for column = 1:size(probmap,2)
    
    string = ['Getting subject patches from patch cluster ' num2str(column) ' of ' num2str(size(probmap,2))];
    fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
    prevstring = string;
    
    this_probmap = probmap(:,column);
    
    
    for communitynum = 1:length(communities)
        community = communities(communitynum);
        
        communitymap = zeros(size(subnetworks_bycluster,1),1);
        
        patchnums_incommunity = find(simple_assigns==community);
        
        for p = 1:length(patchnums_incommunity)
            patchnum = patchnums_incommunity(p);
            sub = clusters_subs_IDs_dices_SAs(patchnum,1);
            communitymap = communitymap + (subnetworks_bycluster(:,sub)==patchnum);
        end
        
        rval = paircorr_mod(communitymap,this_probmap);
        
        if rval>.999
            submaps = false(size(subnetworks_bycluster));
            for p = 1:length(patchnums_incommunity)
                patchnum = patchnums_incommunity(p);
                sub = clusters_subs_IDs_dices_SAs(patchnum,1);
                submaps(:,sub) = submaps(:,sub) | (subnetworks_bycluster(:,sub)==patchnum);
            end
            
            matchedpatches_bysub(:,:,column) = submaps;
            
            communities(communitynum) = [];
            
            break
        else
            clear submaps
        end
    end
    if ~exist('submaps')
        error(['No matching probability map for column ' num2str(column)])
    end
end
disp(' ')
%%
prevstring = [];

groupsystems = ft_read_cifti_mod('/data/cn4/evan/RestingState/Consensus/120_LR_minsize400_recolored_manualconsensus_LR_cleaned.dtseries.nii');
groupsystems = groupsystems.data;
IDs = unique(groupsystems); IDs(IDs<1) = [];

clusteredgroupsystems = false(size(groupsystems,1),0);
for ID = IDs(:)'
    outputclusters = cifti_cluster(groupsystems,ID-.5,ID+.5,0);
    clusteredgroupsystems(:,end+1:end+size(outputclusters,2)) = logical(outputclusters);
end


conditioncounts = cell(length(tasks),1);
conditioncounter = 0;
for tasknum = 1:length(tasks)
    taskdata = ft_read_cifti_mod(['/data/cn4/evan/HCP_tasks/' subs{1} '_tfMRI_' tasks{tasknum} '_level2_hp200_s4.dscalar.nii']);
    
    conditioncounts{tasknum} = [(conditioncounter+1) : (conditioncounter + (size(taskdata.data,2) / 2))];
    conditioncounter = conditioncounter + (size(taskdata.data,2) / 2);
    
end

activationmatrix = zeros(size(probmap,2),conditioncounter,length(subs));
groupactivationmatrix= zeros(size(clusteredgroupsystems,2),conditioncounter,length(subs));

for subnum = 1:length(subs)
    
    for tasknum = 1:length(tasks)
        
        string = ['Calculating activation for subject ' num2str(subnum) ' of ' num2str(length(subs)) ': ' tasks{tasknum}];
        fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
        prevstring = string;
        
        if exist(['/data/cn4/evan/HCP_tasks/' subs{subnum} '_tfMRI_' tasks{tasknum} '_level2_hp200_s4.dscalar.nii'])
            taskdata = ft_read_cifti_mod(['/data/cn4/evan/HCP_tasks/' subs{subnum} '_tfMRI_' tasks{tasknum} '_level2_hp200_s4.dscalar.nii']);
            taskdata = taskdata.data(1:ncortverts,1:(size(taskdata.data,2) / 2));
            for conditionnum = 1:size(taskdata,2)
                conditionindex = conditioncounts{tasknum}(conditionnum);
                for column = 1:size(probmap,2)
                    if any(matchedpatches_bysub(1:ncortverts,subnum,column))
                        activationmatrix(column,conditionindex,subnum) = mean(taskdata(matchedpatches_bysub(1:ncortverts,subnum,column),conditionnum),1);
                    else
                        activationmatrix(column,conditionindex,subnum) = NaN;
                    end
                end
                
                for column = 1:size(clusteredgroupsystems,2)
                    groupactivationmatrix(column,conditionindex,subnum) = mean(taskdata(clusteredgroupsystems(1:ncortverts,column),conditionnum),1);
                end
                
            end
        else
            activationmatrix(:,conditioncounts{tasknum},subnum) = NaN;
            groupactivationmatrix(:,conditioncounts{tasknum},subnum) = NaN;
        end
        
    end
end
%figure;imagesc(nanmean(activationmatrix,3))

thresh = 1;

activation_similarity = zeros(length(subs));
activation_overlap = zeros(length(subs));
groupactivation_similarity = zeros(length(subs));
groupactivation_overlap = zeros(length(subs));
for i = 1:length(subs)
    for j = i:length(subs)
        activationi = activationmatrix(:,:,i); activationi = activationi(:);
        activationj = activationmatrix(:,:,j); activationj = activationj(:);
        both_inds = logical((~isnan(activationi)) .* (~isnan(activationj)));
        activation_similarity(i,j) = paircorr_mod(activationi(both_inds),activationj(both_inds));
        activation_overlap(i,j) = nnz((activationi(both_inds)>thresh) & (activationj(both_inds)>thresh)) / nnz((activationi(both_inds)>thresh) | (activationj(both_inds)>thresh));
        
        groupactivationi = groupactivationmatrix(:,:,i); groupactivationi = groupactivationi(:);
        groupactivationj = groupactivationmatrix(:,:,j); groupactivationj = groupactivationj(:);
        both_inds = logical((~isnan(groupactivationi)) .* (~isnan(groupactivationj)));
        groupactivation_similarity(i,j) = paircorr_mod(groupactivationi(both_inds),groupactivationj(both_inds));
        groupactivation_overlap(i,j) = nnz((groupactivationi(both_inds)>thresh) & (groupactivationj(both_inds)>thresh)) / nnz((groupactivationi(both_inds)>thresh) | (groupactivationj(both_inds)>thresh));
    end
end

figure; imagesc(activation_similarity,[.2 .7])
figure; imagesc(groupactivation_similarity,[.2 .7])
figure; imagesc(activation_similarity-groupactivation_similarity,[-.2 .2])


figure; imagesc(activation_overlap,[.2 .7])
figure; imagesc(groupactivation_overlap,[.2 .7])
figure; imagesc(activation_overlap-groupactivation_overlap,[-.2 .2])
        

disp(' ')







