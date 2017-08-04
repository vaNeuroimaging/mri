%function label = watershed_algorithm_merge4(edgemetricname,outputdir,filestem,hem,thresh)

 edgemetricname = 'PoldromeC3_L_smoothed.func.gii';
% %'/data/cn4/laumannt/left_hem_edge/vc32347_avg_edge_avg_smooth_L_noalone.func.gii';
outputdir = './';
filestem = 'PoldromeC3_L_smoothed';
hem = 'L';

thresh = 1.78;
stepnum = 200;
fracmaxh = 1;
neighdist = 1;
minimathresh = .12;
%percentvertices_insideparcel = 1;
edgeoutlierthresh = 2.5;
parceloutlierthresh = 2.5;
edgeval_hysteresis = .15;%.165;
minparcelsize = 25;

% thresh = 1.45;
% stepnum = 200;
% fracmaxh = 1;
% neighdist = 1;
% minimathresh = .12;
% percentvertices_insideparcel = 1;
% edgeval_hysteresis = .165;
% minparcelsize = 20;

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

%%

disp('Finding minima')

metric = gifti(edgemetricname);
metric = metric.cdata;

metric(1:12) = 100000;

metric(medialindices) = 100000;

%[temp metric] = upper_completion(metric);
clear temp;

%save(gifti(single(metric)),[outputdir '/avgcorrofcorr_smooth2.55_allgrad_L_smooth2.55_edge_avg_uc.func.gii']);


minimametric = zeros(size(metric));
for i = 1:length(metric)
    
    nodeneigh = i;
    newneigh = nodeneigh;
    curneigh = newneigh;
    
    for n = 1:neighdist
        
        for t = 1:length(curneigh)
            newneigh = [newneigh neighbors(curneigh(t),2:7)];
            newneigh(isnan(newneigh)) = [];
        end
        curneigh = setdiff(newneigh,nodeneigh);
        nodeneigh = union(nodeneigh,newneigh);
        
    end
    
    ind = find(nodeneigh==i); %Which nodeneigh is original point
    [minval mini] = min(metric(nodeneigh));
    minindices = find(metric(nodeneigh)==minval);
    
    
    if minval == metric(i)
        minimametric(i) = 1;
        if length(minindices) > 1
            
            minindices(logical(minindices==ind)) = [];
            metric(nodeneigh(minindices)) = minval+.00001;

        end
    end
    
end

clear metric


%%
edgemetric = gifti(edgemetricname); edgemetric = edgemetric.cdata;
origedgemetric = edgemetric;
%edgemetric(logical((edgemetric > edgeval_hysteresis) .* (edgemetric < .4))) = .4;
edgemetric(logical(edgemetric > edgeval_hysteresis)) = edgemetric(logical(edgemetric > edgeval_hysteresis)) + .1;
%minimametric = gifti(minimametricname); minimametric = minimametric.cdata;
minimametric(1:12) = 0;
minimametric(edgemetric>minimathresh) = 0;

sortedge = unique(sort(edgemetric,'descend'));
%[sortedge sortedgepos] = sort(edgemetric,'descend');
%Label initial markers with unique value
labelnum = nnz(minimametric);
labelpos = find(minimametric==1);
label = zeros(size(minimametric));
%randomlabels = randperm(labelnum);
[ign sortorder] = sort(edgemetric(labelpos));
for j = 1:labelnum;
    %label(labelpos(j)) = randomlabels(j);
    label(labelpos(j)) = sortorder(j);
end

minh = sortedge(1);
maxh = sortedge(end);

stoph = sortedge(round(length(sortedge)*fracmaxh));
step = (maxh-minh)/stepnum;
hiter = minh:step:stoph;

%for i = 1:length(sortedge)
for i = 1:length(hiter);
    
    string{i} = ['Running watershed iteration ' num2str(i) ' out of ' num2str(length(hiter))];
    if i==1; fprintf('%s',string{i}); else fprintf([repmat('\b',1,length(string{i-1})) '%s'],string{i}); end
    
    %maskpos = find(edgemetric<sortedge(i)); % Take values in metric less than current iteration
    %maskpos = find(edgemetric<hiter(i)); % Take values in metric less than current iteration
    maskpos = find((edgemetric<hiter(i)) .* (label==0)); % Take values in metric less than current iteration
    
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

% disp('Calculating subject connectivities')
% 
% cohortfile = '/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/NEW_nokids_TMASKLIST.txt';
% [subjects tmasks] = textread(cohortfile,'%s %s');
% subjectnum = find(strcmp(subject,subjects));
% 
% subjectdata = ['/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/cifti_timeseries/' subject '_BOLD_L_surf_subcort_smooth2.55_32k_fsLR.dtseries.nii'];
% copyfile(subjectdata,'/data/cn4/evan/Temp/Temp.dtseries.nii');
% subjectdata = cifti_read('/data/cn4/evan/Temp/Temp.dtseries.nii');
% tmask = load(tmasks{subjectnum});
% subjectdata = subjectdata(:,logical(tmask));
% 
% subcorrelmat = FisherTransform(paircorr_mod(subjectdata'));
% subcorrelmat(isnan(subcorrelmat)) = 0;
%save(gifti(label),[outputdir '/' filestem 'watershedmerge_orig.func.gii']);

%%

label = origlabel;

nomerges = 0;
iteration = 0;

while nomerges==0
    nomerges=1;
    iteration = iteration+1;

    fprintf('%s',['Merging watersheds: iteration ' num2str(iteration) '.....'])
%     string{iteration} = ['Merging watersheds: iteration ' num2str(iteration) '.....'];
%     if iteration==1; fprintf('%s',string{iteration}); else fprintf([repmat('\b',1,length(string{iteration-1})) '%s'],string{iteration}); end
    

watersheds = unique(label);
watersheds(watersheds==0) = [];
for waternum = 1:length(watersheds)
    
    adjacentwatersheds{waternum} = [];
%    edgeval_withinwater(waternum) = mean(edgemetric(label==watersheds(waternum)),1);
     valueswithinthiswatershed = sort(origedgemetric(label==watersheds(waternum)));
     valueswithinthiswatershed(abs(zscore(valueswithinthiswatershed))>parceloutlierthresh) = [];
     %edgeval_withinwater(waternum) = mean(valueswithinthiswatershed(1:(round(length(valueswithinthiswatershed)*percentvertices_insideparcel))),1);
     edgeval_withinwater(waternum) = mean(valueswithinthiswatershed);
    
    for waternum2 = 1:length(watersheds)
        watershedborders{waternum,waternum2} = [];
    end
    
end

borderindices = intersect(find(label==0),corticalindices)';
for bordervertex = borderindices
        
        borderneighs = neighbors(bordervertex,2:7);
        borderneighs(isnan(borderneighs)) = [];
        borderneighvals = label(borderneighs);
        borderneighvals(borderneighvals==0) = [];
        
        watershedneighbors = unique(borderneighvals);
        
        for waterneighbornum = 1:length(watershedneighbors)
            
            thiswatershedindex(waterneighbornum) = find(watersheds==watershedneighbors(waterneighbornum));
            
            otherwatershedneighbors = watershedneighbors;
            otherwatershedneighbors(otherwatershedneighbors==watershedneighbors(waterneighbornum)) = [];
            
            %watershedborders{thiswatershedindex(waterneighbornum)} = unique([watershedborders{thiswatershedindex(waterneighbornum)} bordervertex]);
            adjacentwatersheds{thiswatershedindex(waterneighbornum)} = unique([adjacentwatersheds{thiswatershedindex(waterneighbornum)} otherwatershedneighbors']);
        end
        
        if length(watershedneighbors) == 2 && length(borderneighvals)>2
        
            watershedborders{thiswatershedindex(1),thiswatershedindex(2)} = unique([watershedborders{thiswatershedindex(1),thiswatershedindex(2)} bordervertex]);
            watershedborders{thiswatershedindex(2),thiswatershedindex(1)} = unique([watershedborders{thiswatershedindex(2),thiswatershedindex(1)} bordervertex]);
        
        end
        
end
merges = 0;
for wateri = 1:length(watersheds)
    for waterj = 1:length(watersheds)
        if wateri < waterj
        %sharedborder = intersect(watershedborders{wateri},watershedborders{waterj});
        
        if length(watershedborders{wateri,waterj}) > 1
            
            edgevals_thisborder = edgemetric(watershedborders{wateri,waterj});
            edgevals_thisborder((abs(zscore(edgevals_thisborder)) > edgeoutlierthresh)) = [];
            
            edgeval_thisborder = mean(edgevals_thisborder,1);
            if (edgeval_withinwater(wateri)*thresh > edgeval_thisborder) && (edgeval_withinwater(waterj)*thresh > edgeval_thisborder)
                label(label==watersheds(waterj)) = mean(label(origlabel==watersheds(wateri)),1);
                nomerges = 0;
                merges = merges+1;
            end
        end
        end
    end
end
            
    fprintf('%s',[num2str(merges) ' watersheds merged'])
    %string{iteration} = [string{iteration} num2str(merges) ' watersheds merged'];
    disp(' ')
    
    for bordervertex = borderindices
        
        borderneighs = neighbors(bordervertex,2:7);
        borderneighs(isnan(borderneighs)) = [];
        borderneighvals = label(borderneighs);
        borderneighvals(borderneighvals==0) = [];
        if length(unique(borderneighvals)) == 1;
            label(bordervertex) = unique(borderneighvals);
        end
    end
    
    clear watershedborders adjacentwatersheds edgeval_withinwater
    
%     if merges > 0
%         save(gifti(label),[outputdir '/' filestem 'watershedmerge_iter' num2str(iteration) '.func.gii']);
%         
%     end
    
end

%disp(' ')
disp('Removing small parcels')

for repeatthis = 1:2

watersheds = unique(label);
watersheds(watersheds==0) = [];
for waternum = 1:length(watersheds)
    
    adjacentwatersheds{waternum} = [];
    %edgeval_withinwater(waternum) = mean(edgemetric(label==watersheds(waternum)),1);
    valueswithinthiswatershed = sort(origedgemetric(label==watersheds(waternum)));
    valueswithinthiswatershed(abs(zscore(valueswithinthiswatershed))>parceloutlierthresh) = [];
    %edgeval_withinwater(waternum) = mean(valueswithinthiswatershed(1:(round(length(valueswithinthiswatershed)*percentvertices_insideparcel))),1);
    edgeval_withinwater(waternum) = mean(valueswithinthiswatershed);
    
    for waternum2 = 1:length(watersheds)
        watershedborders{waternum,waternum2} = [];
    end
    
end

borderindices = intersect(find(label==0),corticalindices)';
for bordervertex = borderindices
        
        borderneighs = neighbors(bordervertex,2:7);
        borderneighs(isnan(borderneighs)) = [];
        borderneighvals = label(borderneighs);
        borderneighvals(borderneighvals==0) = [];
        
        watershedneighbors = unique(borderneighvals);
        
        for waterneighbornum = 1:length(watershedneighbors)
            
            thiswatershedindex(waterneighbornum) = find(watersheds==watershedneighbors(waterneighbornum));
            
            otherwatershedneighbors = watershedneighbors;
            otherwatershedneighbors(otherwatershedneighbors==watershedneighbors(waterneighbornum)) = [];
            
            %watershedborders{thiswatershedindex(waterneighbornum)} = unique([watershedborders{thiswatershedindex(waterneighbornum)} bordervertex]);
            adjacentwatersheds{thiswatershedindex(waterneighbornum)} = unique([adjacentwatersheds{thiswatershedindex(waterneighbornum)} otherwatershedneighbors']);
        end
        
        if length(watershedneighbors) == 2
        
            watershedborders{thiswatershedindex(1),thiswatershedindex(2)} = unique([watershedborders{thiswatershedindex(1),thiswatershedindex(2)} bordervertex]);
            watershedborders{thiswatershedindex(2),thiswatershedindex(1)} = unique([watershedborders{thiswatershedindex(2),thiswatershedindex(1)} bordervertex]);
        
        end
        
end


for waternum = 1:length(watersheds)

    if (length(find(label==watersheds(waternum))) > 0) && (length(find(label==watersheds(waternum))) < minparcelsize);
        minedge = 1;
        minwaternum = 0;
        for i = 1:length(adjacentwatersheds{waternum})
            thisadjacentwatershedindex = find(watersheds == adjacentwatersheds{waternum}(i));
            if mean(edgemetric(watershedborders{waternum,thisadjacentwatershedindex}),1) < minedge
                minwaternum = thisadjacentwatershedindex;
                minedge = mean(edgemetric(watershedborders{waternum,thisadjacentwatershedindex}),1);
            end
        end
        if minwaternum
            label(label==watersheds(waternum)) = watersheds(minwaternum);
        else
            label(label==watersheds(waternum)) = 0;
        end
    end
end

for bordervertex = borderindices
        
        borderneighs = neighbors(bordervertex,2:7);
        borderneighs(isnan(borderneighs)) = [];
        borderneighvals = label(borderneighs);
        borderneighvals(borderneighvals==0) = [];
        if length(unique(borderneighvals)) == 1;
            label(bordervertex) = unique(borderneighvals);
        end
    end
end

save(gifti(label),[outputdir '/' filestem 'watershedmerge.func.gii']);
