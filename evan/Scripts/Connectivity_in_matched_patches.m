
probmapfile = 'Cluster_probability_maps_sorted_10mm_103sub_selected.dscalar.nii';
tmaskfile = '/data/cn5/selfRegulation/V4Process_nosmooth/Finaltmasklist_120_108.txt';
datalistfile = '/data/cn5/selfRegulation/V4Process_nosmooth/Surfdatalist_120_108.txt';
hcptmaskfile = '/data/cn4/evan/RestingState/Ind_variability/HCP/HCP_80_TMASKLIST.txt';
hcpdatalistfile = '/data/cn4/evan/RestingState/Ind_variability/HCP/HCP_80_DATALIST.txt';

ncortverts = 59412;

[subs tmasks] = textread(tmaskfile,'%s%s');
[subs ciftifiles] = textread(datalistfile,'%s%s');

[hcpsubs hcptmasks] = textread(hcptmaskfile,'%s%s');
[hcpsubs hcpciftifiles] = textread(hcpdatalistfile,'%s%s');

subs = [subs ; hcpsubs];
tmasks = [tmasks ; hcptmasks];
ciftifiles = [ciftifiles ; hcpciftifiles];


%substouse = 1:120;
%subs = subs(substouse); tmasks = tmasks(substouse); ciftifiles = ciftifiles(substouse);


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



thresh = 1;

connectivity_similarity = zeros(length(subs));
groupconnectivity_similarity = zeros(length(subs));

prevstring = [];
for subnum1 = 1:length(subs)
    
    s1data = ft_read_cifti_mod(ciftifiles{subnum1});
    tmask = load(tmasks{subnum1});
    s1data = s1data.data(1:ncortverts,logical(tmask));
    
    
    for subnum2 = (subnum1+1):length(subs)
        
            
            string = ['Calculating pairwise connectivity in patches matched between subjects ' num2str(subnum1) ' and ' num2str(subnum2)];
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
            
            s1_timecourses = zeros(size(s1data,2),size(s1ROIs,2));
            for roinum = 1:size(s1ROIs,2)
                s1_timecourses(:,roinum) = mean(s1data(logical(s1ROIs(:,roinum)),:),1);
            end
            s1_connectivity = paircorr_mod(s1_timecourses);
            s1_timecourses = zeros(size(s1data,2),size(clusteredgroupsystems,2));
            for roinum = 1:size(clusteredgroupsystems,2)
                s1_timecourses(:,roinum) = mean(s1data(logical(clusteredgroupsystems(:,roinum)),:),1);
            end
            s1_groupconnectivity = paircorr_mod(s1_timecourses);
                        
            s2data = ft_read_cifti_mod(ciftifiles{subnum2});
            tmask = load(tmasks{subnum2});
            s2data = s2data.data(1:ncortverts,logical(tmask));
            s2_timecourses = zeros(size(s2data,2),size(s2ROIs,2));
            for roinum = 1:size(s2ROIs,2)
                s2_timecourses(:,roinum) = mean(s2data(logical(s2ROIs(:,roinum)),:),1);
            end
            s2_connectivity = paircorr_mod(s2_timecourses);
            s2_timecourses = zeros(size(s2data,2),size(clusteredgroupsystems,2));
            for roinum = 1:size(clusteredgroupsystems,2)
                s2_timecourses(:,roinum) = mean(s2data(logical(clusteredgroupsystems(:,roinum)),:),1);
            end
            s2_groupconnectivity = paircorr_mod(s2_timecourses);
            
            connectivitymask = triu(true(size(s1_connectivity)),1);
            groupconnectivitymask = triu(true(size(s1_groupconnectivity)),1);
            
            connectivity_similarity(subnum1,subnum2) = paircorr_mod(s1_connectivity(connectivitymask),s2_connectivity(connectivitymask));
            groupconnectivity_similarity(subnum1,subnum2) = paircorr_mod(s1_groupconnectivity(groupconnectivitymask),s2_groupconnectivity(groupconnectivitymask));
            
            clear s1_timecourses s1_connectivity s1_groupconnectivity s2data s2_timecourses s2_connectivity s2_groupconnectivity 
            
    end
end

figure; imagesc(connectivity_similarity,[.2 .7])
figure; imagesc(groupconnectivity_similarity,[.2 .7])
figure; imagesc(connectivity_similarity-groupconnectivity_similarity,[-.2 .2])


%figure; imagesc(activation_overlap,[.2 .7])
%figure; imagesc(groupactivation_overlap,[.2 .7])
%figure; imagesc(activation_overlap-groupactivation_overlap,[-.2 .2])
        

disp(' ')







