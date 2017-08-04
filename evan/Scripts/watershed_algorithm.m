function label = watershed_algorithm(edgemetric,minimametric,stepnum,fracmaxh,outputdir,filestem)

bufsize=16384;
%caredir = '/data/cn4/laumannt/assignment_problem_v2/caret2/PALS_B12.BOTH-HEMS.CLEAN.73730/LEFT/18K';
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
cd(caretdir)
% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread('node_neighbors.txt','%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;

%dir = '/data/cn4/laumannt/fcMapping_redux/fcdata_tmasked';
%dir = '/data/hcp-bluearc/home/laumannt/monkey_fcmri/OHSU';
%cd(dir)
%temp = gifti('AllC_avg_edge_avg_smooth_L.metric');
%temp = gifti('DHA_avg_edge_avg_smooth_L.metric');
%edgemetric = temp.cdata;
%temp = gifti('AllC_avg_edge_avg_smooth_extrema3_L.metric');
%temp = gifti('DHA_avg_edge_avg_smooth_extrema4_L.metric');
%minimametric = single(temp.cdata<0);

%filestem = 'DHA_extrema3_';

sortedge = unique(sort(edgemetric,'descend'));
%[sortedge sortedgepos] = sort(edgemetric,'descend');
%Label initial markers with unique value
labelnum = nnz(minimametric);
labelpos = find(minimametric==1);
label = zeros(size(minimametric));
for j = 1:labelnum;
label(labelpos(j)) = j;
end

minh = sortedge(1);
maxh = sortedge(end);

stoph = sortedge(round(length(sortedge)*fracmaxh));
step = (maxh-minh)/stepnum;
hiter = minh:step:stoph;

%for i = 1:length(sortedge)
 for i = 1:length(hiter);
    disp(['Number of iterations will be ' num2str(length(sortedge)) ', Iteration = ' num2str(i)])
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
%   
  save(gifti(label),[outputdir '/' filestem 'watershed.gii'])
  system(['mv ' outputdir '/' filestem 'watershed.gii ' outputdir '/' filestem 'watershed.func.gii'])