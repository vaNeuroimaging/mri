%function Activation_in_clustered_patches(probmapfile)

probmapfile = 'Cluster_probability_maps_sorted_10mm_27sub_selected.dscalar.nii';
tmaskfile = '/data/cn4/evan/RestingState/Ind_variability/HCP/HCP_80_TMASKLIST.txt';
tasks = {'EMOTION','GAMBLING','LANGUAGE','MOTOR','RELATIONAL','SOCIAL','WM'};

ncortverts = 59412;

[subs ign] = textread(tmaskfile,'%s%s');


baddata = ft_read_cifti_mod('/data/cn4/evan/ROIs/Baddata_bigcluster_LR.dtseries.nii');
baddata = logical(baddata.data(1:ncortverts,:));

%%


load System_clustering.mat
load match_matrix.mat

subnetworks_bycluster = subnetworks_bycluster(1:ncortverts,:);

%%

probmaps_template = ft_read_cifti_mod(probmapfile);
probmap = probmaps_template.data;
probmaps_template.data = [];
prevstring = [];

simple_assigns = load('rawassn_minsize2.txt');
assigns_inprobmaps = zeros(size(simple_assigns));

communities = unique(simple_assigns); communities(communities<1) = [];

matchedpatches_bysub = false(size(probmap,1),length(subs),size(probmap,2));

for column = 1:size(probmap,2)
    
    string = ['Getting subject patches from patch cluster ' num2str(column) ' of ' num2str(size(probmap,2))];
    fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
    prevstring = string;
    
    this_probmap = probmap(1:ncortverts,column);
    
    
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
            assigns_inprobmaps(simple_assigns==community) = community;
            
            break
        
        end
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

%activationmatrix = zeros(size(probmap,2),conditioncounter,length(subs));
%groupactivationmatrix= zeros(size(clusteredgroupsystems,2),conditioncounter,length(subs));

%figure;imagesc(nanmean(activationmatrix,3))

thresh = 1;

activation_similarity = zeros(length(subs));
groupactivation_similarity = zeros(length(subs));

prevstring = [];
for subnum1 = 23%1:length(subs)
    
    for tasknum = 1:length(tasks)
        if exist(['/data/cn4/evan/HCP_tasks/' subs{subnum1} '_tfMRI_' tasks{tasknum} '_level2_hp200_s4.dscalar.nii'])
            taskdata = ft_read_cifti_mod(['/data/cn4/evan/HCP_tasks/' subs{subnum1} '_tfMRI_' tasks{tasknum} '_level2_hp200_s4.dscalar.nii']);
            taskdata_s1{tasknum} = taskdata.data(1:ncortverts,1:(size(taskdata.data,2) / 2 -1));
        else
            taskdata_s1{tasknum} = [];
        end
    end
    
    
    for subnum2 = 64%(subnum1+1):length(subs)
        
            
            string = ['Calculating activation in patches matched between subjects ' num2str(subnum1) ' and ' num2str(subnum2)];
            fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
            prevstring = string;
            
            patchcounter = 0;
            
            s1patchnums = find(clusters_subs_IDs_dices_SAs(:,1)==subnum1);
            used_s1patchnums = [];
            s2patchnums = find(clusters_subs_IDs_dices_SAs(:,1)==subnum2);
            
            s1ROIs = false(size(subnetworks_bycluster,1),0);
            s2ROIs = false(size(subnetworks_bycluster,1),0);
            
            for i = 1:length(s1patchnums)
                these_s1patches = s1patchnums(i);
                these_s2patches = intersect(s2patchnums,find(match_matrix(these_s1patches,:)));
                if (~any(used_s1patchnums==these_s1patches)) && (~isempty(these_s2patches))
                    
                    these_s1patches = unique([these_s1patches; intersect(s1patchnums,find(sum(match_matrix(these_s2patches,:),1)))]);
                    used_s1patchnums = [used_s1patchnums these_s1patches(:)'];
                    
                    if all(assigns_inprobmaps([these_s1patches(:) ; these_s2patches(:)]))
                    
                    patchcounter = patchcounter +1;
                    
                    s1ROIs(:,patchcounter) = false;
                    s2ROIs(:,patchcounter) = false;
                                       
                    for j = these_s1patches(:)'
                        s1ROIs(:,patchcounter) = s1ROIs(:,patchcounter) | (subnetworks_bycluster(:,subnum1)==j);
                    end
                    
                    for j = these_s2patches(:)'
                        s2ROIs(:,patchcounter) = s2ROIs(:,patchcounter) | (subnetworks_bycluster(:,subnum2)==j);
                    end
                    
                    if (nnz(baddata & (s1ROIs(:,patchcounter))) / nnz(s1ROIs(:,patchcounter)) > .1) || (nnz(baddata & (s2ROIs(:,patchcounter))) / nnz(s2ROIs(:,patchcounter)) > .1)
                        s1ROIs(:,patchcounter) = [];
                        s2ROIs(:,patchcounter) = [];
                        patchcounter = patchcounter - 1;
                    end
                    end
                    
                end
            end
            
            activations_thesetwosubs = zeros(patchcounter,conditioncounter,2);
            groupactivations_thesetwosubs= zeros(size(clusteredgroupsystems,2),conditioncounter,2);
            
            for tasknum = 1:length(tasks)
                
                string = ['Calculating activation in patches matched between subjects ' num2str(subnum1) ' and ' num2str(subnum2) ': ' tasks{tasknum}];
                fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
                prevstring = string;
                
                if exist(['/data/cn4/evan/HCP_tasks/' subs{subnum1} '_tfMRI_' tasks{tasknum} '_level2_hp200_s4.dscalar.nii']) && exist(['/data/cn4/evan/HCP_tasks/' subs{subnum2} '_tfMRI_' tasks{tasknum} '_level2_hp200_s4.dscalar.nii'])
                    %taskdata = ft_read_cifti_mod(['/data/cn4/evan/HCP_tasks/' subs{subnum1} '_tfMRI_' tasks{tasknum} '_level2_hp200_s4.dscalar.nii']);
                    %taskdata = taskdata.data(1:ncortverts,1:(size(taskdata.data,2) / 2 -1));
                    taskdata = taskdata_s1{tasknum};
                    for conditionnum = 1:size(taskdata,2)
                        conditionindex = conditioncounts{tasknum}(conditionnum);
                        for patchnum = 1:patchcounter
                            activations_thesetwosubs(patchnum,conditionindex,1) = mean(taskdata(s1ROIs(:,patchnum),conditionnum),1);
                        end
                        for patchnum = 1:size(clusteredgroupsystems,2)
                            groupactivations_thesetwosubs(patchnum,conditionindex,1) = mean(taskdata(clusteredgroupsystems(1:ncortverts,patchnum),conditionnum),1);
                        end
                    end
                    
                    taskdata = ft_read_cifti_mod(['/data/cn4/evan/HCP_tasks/' subs{subnum2} '_tfMRI_' tasks{tasknum} '_level2_hp200_s4.dscalar.nii']);
                    taskdata = taskdata.data(1:ncortverts,1:(size(taskdata.data,2) / 2 -1));
                    for conditionnum = 1:size(taskdata,2)
                        conditionindex = conditioncounts{tasknum}(conditionnum);
                        for patchnum = 1:patchcounter
                            activations_thesetwosubs(patchnum,conditionindex,2) = mean(taskdata(s2ROIs(:,patchnum),conditionnum),1);
                        end
                        for patchnum = 1:size(clusteredgroupsystems,2)
                            groupactivations_thesetwosubs(patchnum,conditionindex,2) = mean(taskdata(clusteredgroupsystems(1:ncortverts,patchnum),conditionnum),1);
                        end
                    end
                else
                    activations_thesetwosubs(:,conditioncounts{tasknum},:) = NaN;
                    groupactivations_thesetwosubs(:,conditioncounts{tasknum},:) = NaN;
                end
            end
            
            nonnaninds = (~isnan(activations_thesetwosubs(:,:,1)));
            
            s1activation = activations_thesetwosubs(:,:,1);
            s1activation = s1activation(nonnaninds);
            s2activation = activations_thesetwosubs(:,:,2);
            s2activation = s2activation(nonnaninds);
            
            if ~isempty(s1activation);
            
                activation_similarity(subnum1,subnum2) = corr(s1activation,s2activation);
                
            end
            
            
            s1activation = groupactivations_thesetwosubs(:,:,1);
            s1activation = s1activation(:);
            s2activation = groupactivations_thesetwosubs(:,:,2);
            s2activation = s2activation(:);
            
            groupactivation_similarity(subnum1,subnum2) = corr(s1activation,s2activation);
            
       
    end
end

figure; imagesc(activation_similarity,[.2 .7])
figure; imagesc(groupactivation_similarity,[.2 .7])
figure; imagesc(activation_similarity-groupactivation_similarity,[-.2 .2])


%figure; imagesc(activation_overlap,[.2 .7])
%figure; imagesc(groupactivation_overlap,[.2 .7])
%figure; imagesc(activation_overlap-groupactivation_overlap,[-.2 .2])
        

disp(' ')







