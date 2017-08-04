%function label = watershed_algorithm_merge2(edgemetricname,minimametricname,outputdir,filestem,subject,hem,thresh)

edgemetricname = '/data/cn4/laumannt/left_hem_edge/vc32347_avg_edge_avg_smooth_L_noalone.func.gii';
minimametricname = '/data/cn4/evan/RestingState/Ind_variability/vc32347_minima1.func.gii';
% stepnum = 200;
% fracmaxh = 1;
outputdir = '/data/cn4/evan/RestingState/Ind_variability/';
filestem = 'vc32347_minima1';
subject = 'vc32347';
hem = 'L';
thresh = .8;
stepnum = 100;
fracmaxh = 1;

bufsize=16384;
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;

medialmaskdata = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hem '.32k_fs_LR.func.gii']);
medialmaskdata = medialmaskdata.cdata;
corticalindices = find(medialmaskdata==0);
medialindices = find(medialmaskdata);

edgemetric = gifti(edgemetricname); edgemetric = edgemetric.cdata;
minimametric = gifti(minimametricname); minimametric = minimametric.cdata;

sortedge = unique(sort(edgemetric,'descend'));
%[sortedge sortedgepos] = sort(edgemetric,'descend');
%Label initial markers with unique value
labelnum = nnz(minimametric);
labelpos = find(minimametric==1);
label = zeros(size(minimametric));
randomlabels = randperm(labelnum);
for j = 1:labelnum;
    label(labelpos(j)) = randomlabels(j);
end

minh = sortedge(1);
maxh = sortedge(end);

stoph = sortedge(round(length(sortedge)*fracmaxh));
step = (maxh-minh)/stepnum;
hiter = minh:step:stoph;

%for i = 1:length(sortedge)
for i = 1:length(hiter);
    
    string{i} = ['Running iteration ' num2str(i) ' out of ' num2str(length(hiter))];
    if i==1; fprintf('%s',string{i}); else fprintf([repmat('\b',1,length(string{i-1})) '%s'],string{i}); end
    
    %maskpos = find(edgemetric<sortedge(i)); % Take values in metric less than current iteration
    maskpos = find(edgemetric<hiter(i)); % Take values in metric less than current iteration
    
    for m = 1:length(maskpos) %For all nodes at this threshold
        nodeneigh = neighbors(maskpos(m),2:7);
        %nodeneigh = neighbors(sortedgepos(i),2:7);
        nodeneigh(isnan(nodeneigh)) = [];
        nodeneighlab = label(nodeneigh);
        
        if nnz(nodeneighlab)>0 % If there are neighbors who are labeled
            nodeneighlab(nodeneighlab==0)=[];
            if length(unique(nodeneighlab))>1 % If neighbors have more than one label, then watershed node
                label(maskpos(m)) = 0;
                %label(sortedgepos(i)) = 0;
            else % If neighbors only have one label than join them
                label(maskpos(m)) = unique(nodeneighlab);
                %label(sortedgepos(i)) = unique(nodeneighlab);
            end
        end
    end
end

disp(' ')
label(medialindices) = 0;
origlabel = label;

disp('Calculating subject connectivities')

cohortfile = '/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/NEW_nokids_TMASKLIST.txt';
[subjects tmasks] = textread(cohortfile,'%s %s');
subjectnum = find(strcmp(subject,subjects));

subjectdata = ['/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/cifti_timeseries/' subject '_BOLD_L_surf_subcort_smooth2.55_32k_fsLR.dtseries.nii'];
copyfile(subjectdata,'/data/cn4/evan/Temp/Temp.dtseries.nii');
subjectdata = cifti_read('/data/cn4/evan/Temp/Temp.dtseries.nii');
tmask = load(tmasks{subjectnum});
subjectdata = subjectdata(:,logical(tmask));

subcorrelmat = FisherTransform(paircorr_mod(subjectdata'));
subcorrelmat(isnan(subcorrelmat)) = 0;

save(gifti(label),[outputdir '/' filestem 'watershedmerge_orig.func.gii']);

nomerges = 0;
iteration = 0;

while nomerges==0
    nomerges=1;
    iteration = iteration+1;
    disp(['Merging watersheds: iteration ' num2str(iteration)]);
    
    
    watersheds = unique(label);
    watersheds(watersheds==0) = [];
    watershedcorticaldata = label(corticalindices);
    
    for watershednum = 1:length(watersheds)
        
        watershed = watersheds(watershednum);
        watershedindices{watershednum} = find(watershedcorticaldata==watershed);
        
        watershedconnectivity(:,watershednum) = mean(subcorrelmat(:,watershedindices{watershednum}),2);
    end
    
    watershedsimilarity = paircorr_mod(watershedconnectivity);
    watershedsimilarity = (watershedsimilarity>thresh);
    
    borderindices = intersect(find(label==0),corticalindices)';
    merges = 0;
    for bordervertex = borderindices
        
        borderneighs = neighbors(bordervertex,2:7);
        borderneighs(isnan(borderneighs)) = [];
        borderneighvals = label(borderneighs);
        borderneighvals(borderneighvals==0) = [];
        
        watershedneighbors = unique(borderneighvals);
        
        if length(watershedneighbors)>1
            
            for neighbornumi = 1:length(watershedneighbors);
                for neighbornumj = 1:length(watershedneighbors);
                    watershedneighborindex(neighbornumj) = find(watersheds==watershedneighbors(neighbornumj));
                    if neighbornumi<neighbornumj && watershedsimilarity(watershedneighborindex(neighbornumi),watershedneighborindex(neighbornumj));
                        label(label==watershedneighbors(neighbornumi)) = watershedneighbors(neighbornumj);
                        watershedsimilarity(watershedneighborindex(neighbornumi),watershedneighborindex(neighbornumj)) = 0;
                        merges = merges+1;
                        nomerges = 0;
                    end
                end
            end
        end
        clear watershedneighborindex
    end
    
    disp([num2str(merges) ' watersheds merged'])
    
    for bordervertex = borderindices
        
        borderneighs = neighbors(bordervertex,2:7);
        borderneighs(isnan(borderneighs)) = [];
        borderneighvals = label(borderneighs);
        borderneighvals(borderneighvals==0) = [];
        if length(unique(borderneighvals)) == 1;
            label(bordervertex) = unique(borderneighvals);
        end
    end
    
    clear watershedconnectivity watershedindices
    
    save(gifti(label),[outputdir '/' filestem 'watershedmerge_iter' num2str(iteration) '.func.gii']);
    
end


save(gifti(label),[outputdir '/' filestem 'watershedmerge.func.gii']);
