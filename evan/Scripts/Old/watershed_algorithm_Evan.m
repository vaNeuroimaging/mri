function label = watershed_algorithm_Evan(edgemetricname,minimametricname,stepnum,fracmaxh,outputdir,filestem)

%label = watershed_algorithm_Evan(edgemetricname,minimametricname,stepnum,fracmaxh,outputdir,filestem)

% edgemetricname = 'avgcorrofcorr_smooth2.55_allgrad_L_smooth2.55_edge_avg_smooth2.55_lOT_forwatershed.func.gii';
% minimametricname = 'avgcorrofcorr_smooth2.55_allgrad_L_smooth2.55_edge_avg_smooth2.55_minima4_lOT.func.gii';
% stepnum = 200;
% fracmaxh = 1;
% outputdir = '/data/cn4/evan/RestingState/FC_Mapping_120/lOT/Watershed/';
% filestem = 'lOT_';

bufsize=16384;
%caredir = '/data/cn4/laumannt/assignment_problem_v2/caret2/PALS_B12.BOTH-HEMS.CLEAN.73730/LEFT/18K';
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;

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
 regionnums = unique(label(logical(label>0)));
 for i = 1:length(regionnums);
     regionnum = regionnums(i);
     rois(:,regionnum) = single(label==regionnum);
 end

bordermetric = zeros(length(edgemetric),1);

for i = 1:length(edgemetric)
    
    theseneighbors = neighbors(i,2:7);
    theseneighbors(isnan(theseneighbors)) = [];
    
    if label(i) > 0 && any(label(theseneighbors) == 0)
        bordermetric(i) = 1;
    end
        
end

%save(gifti(bordermetric),[outputdir '/' filestem 'watershed_outlines.func.gii']);

save(gifti(label),[outputdir '/' filestem 'watershed.func.gii']);
%save(gifti(rois),[outputdir '/' filestem 'watershed_rois.func.gii'])
  %system(['mv ' outputdir '/' filestem 'watershed.gii ' outputdir '/' filestem 'watershed.func.gii'])