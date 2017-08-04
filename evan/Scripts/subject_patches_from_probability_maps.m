function subject_patches_from_probability_maps(probmapfile,mapnums,patch_probability_map_params)
%subject_patches_from_probability_maps(probmapfile,mapnums,patch_probability_map_params)

%Find params file
[paramspath,paramsname,paramsextension] = fileparts(patch_probability_map_params);
origpath = pwd;
if ~isempty(paramspath)
    cd(paramspath)
end

%Load parameters
params = feval(paramsname);
varnames = fieldnames(params);
for i = 1:length(varnames)
    evalc([varnames{i} ' = params.' varnames{i}]);
end
clear varnames params



probmaps = ft_read_cifti_mod(probmapfile);
cifti_template = probmaps; cifti_template.data = [];


load(system_clustering_names)


simple_assigns = load('rawassn_minsize2.txt');

communities = unique(simple_assigns); communities(communities<1) = [];


for mapnum = mapnums

    rs = zeros(length(communities),1);

    probmap = probmaps.data(:,mapnum);
    
    
    
    for communitynum = 1:length(communities)
        community = communities(communitynum);
        
        communitymap = zeros(size(subnetworks_bycluster,1),1);
        
        patchnums_incommunity = find(simple_assigns==community);
        
        for p = 1:length(patchnums_incommunity)
            patchnum = patchnums_incommunity(p);
            sub = clusters_subs_IDs_dices_SAs(patchnum,1);
            communitymap = communitymap + (subnetworks_bycluster(:,sub)==patchnum);
        end
        
        rs(communitynum) = paircorr_mod(communitymap,probmap);
        
        if rs(communitynum)>.999
            submaps = false(size(subnetworks_bycluster));
            for p = 1:length(patchnums_incommunity)
                patchnum = patchnums_incommunity(p);
                sub = clusters_subs_IDs_dices_SAs(patchnum,1);
                submaps(:,sub) = submaps(:,sub) | (subnetworks_bycluster(:,sub)==patchnum);
            end
            tokens = tokenize(probmaps.mapname{mapnum},':');
            columnnetworkname = tokens{1};
            color = find(strcmp(columnnetworkname,networklabels));
            submaps = submaps .* color;
            break
        else
            clear submaps
        end
    end
    
    if ~exist('submaps')
        error('No matching probability map')
    end
    
    cifti_template.mapname = [];
    for s = 1:size(submaps,1);
        cifti_template.mapname{s} = ['Subject ' num2str(s)];
    end
    cifti_template.data = submaps;
    
    ft_write_cifti_mod(['Submaps_probmap' num2str(mapnum) '_' probmapfile],cifti_template);
    
end



        
        
        
    