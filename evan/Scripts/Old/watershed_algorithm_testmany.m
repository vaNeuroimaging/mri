%function label = watershed_algorithm_vector(edgemetricname,minimametricname,stepnum,fracmaxh,outputdir,filestem)

edgemetricname = 'HCP20_avg_corrofcorr_L_smooth2.55grad_edges_avg.func.gii';
minimametricname = 'HCP20_avg_corrofcorr_L_smooth2.55grad_edges_avg_minima2.func.gii';
stepnum = 200;
fracmaxh = 1;
outputdir = './';
filestem = 'HCP20_avg_corrofcorr_L_smooth2.55grad_edges_avg_smooth2.55';
medialmask = '/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii';
medialmaskdata = gifti(medialmask);
medialind = find(medialmaskdata.cdata~=0);


bufsize=16384;
%caredir = '/data/cn4/laumannt/assignment_problem_v2/caret2/PALS_B12.BOTH-HEMS.CLEAN.73730/LEFT/18K';
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';
% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;
%neighbors(:,1) = [];



edgemetric = gifti(edgemetricname); edgemetric = edgemetric.cdata;

minimametric = gifti(minimametricname); minimametric = minimametric.cdata;

minimametric(medialind) = 0;

minimametrics = minimametric;

%load minimametrics.mat

minimametrics = zeros([length(minimametric) ceil(length(find(minimametric))*4/5)]);
minimametrics(:,1) = minimametric;

for minimaset = 2:size(minimametrics,2)
    
    string{minimaset} = ['Creating minima subset ' num2str(minimaset) ' out of ' num2str(size(minimametrics,2))];
     if minimaset==2; fprintf('%s',string{minimaset}); else fprintf([repmat('\b',1,length(string{minimaset-1})) '%s'],string{minimaset}); end

minimaindices = find(minimametrics(:,minimaset-1));

distancetonextminima = zeros(length(minimaindices),1);

for i = 1:length(minimaindices)
    
    nodeneigh = minimaindices(i);
    newneigh = nodeneigh;
    curneigh = newneigh;
    
    neighborminimafound = 0;
    neighbordist = 1;
    
    while neighborminimafound == 0
    %for n = 1:neighdist
        
        for t = 1:length(curneigh)
            newneigh = [newneigh neighbors(curneigh(t),2:7)];
            newneigh(isnan(newneigh)) = [];
        end
        curneigh = setdiff(newneigh,nodeneigh);
        nodeneigh = union(nodeneigh,newneigh);
        
        if any(minimametric(newneigh))
            distancetonextminima(i) = neighbordist;
            neighborminimafound = 1;
        end
        
    end
    
end

minimavalues = edgemetric(minimaindices);

minimaevaluation = minimavalues ./ (distancetonextminima*10);

[ign worstminimaidx] = max(minimavalues);

newminimametric = minimametrics(:,minimaset-1);
newminimametric(minimaindices(worstminimaidx)) = 0;
minimametrics(:,minimaset) = newminimametric;

end


   disp(' ')


%save('minimametrics.mat','minimametrics')


%%

sortedge = unique(sort(edgemetric,'descend'));
%[sortedge sortedgepos] = sort(edgemetric,'descend');
%Label initial markers with unique value
labelnum = nnz(minimametrics);
labelpos = find(minimametrics==1);
label = zeros(size(minimametrics));
randomlabels = randperm(labelnum);
for j = 1:labelnum;
label(labelpos(j)) = j;%randomlabels(j);
end

minh = sortedge(1);
maxh = sortedge(end);

stopnum = length(sortedge)*fracmaxh;

step = round(stopnum/stepnum);

stepindices = 1:step:round(stopnum);
hiter = sortedge(stepindices);

%stoph = sortedge(stopnum);
% step = (maxh-minh)/stepnum;
% hiter = minh:step:stoph;

for i=1:6
    neighbormat(:,:,i) = repmat(neighbors(:,i+1),[1 size(minimametrics,2)]);
end
for i = 2:size(minimametrics,2)
    neighbormat(:,i,:) = neighbormat(:,i,:) + size(neighbormat,1)*(i-1);
end

nomedialnodesmat = repmat((~logical(medialmaskdata.cdata)),[1 size(minimametrics,2)]);

 nodestomakewatershedzone = zeros(size(label));
 for i = 1:length(hiter);
     
     string{i} = ['Running iteration ' num2str(i) ' out of ' num2str(length(hiter))];
     if i==1; fprintf('%s',string{i}); else fprintf([repmat('\b',1,length(string{i-1})) '%s'],string{i}); end
     
    mask = (repmat(edgemetric,[1 size(minimametrics,2)]))<=(hiter(i)) .* nomedialnodesmat;
    
    if nnz(mask) > 0;
    
        if mask(1833)
            a=1;
        end
        
        
    nodeneighlab = zeros([size(label) 6]);
    
    for neigh=1:6
        tempneighbormat = neighbormat(:,:,neigh);
        theseneighbors = tempneighbormat(1:numel(tempneighbormat));
        
        findazero = find(label==0);
        theseneighbors(isnan(theseneighbors)) = findazero(1);
        
        tempnodeneighlab = zeros(size(label));
        tempnodeneighlab(1:numel(tempnodeneighlab)) = label(theseneighbors);
        
        nodeneighlab(:,:,neigh) = tempnodeneighlab;
        nodeneighlab(isnan(nodeneighlab)) = 0;
        
    end
    
    highestlabelneighbor = max(nodeneighlab,[],3);
    
    nonzeroavglabels = sum(nodeneighlab,3) ./ sum((nodeneighlab>0),3);
    nodestoadd = (nonzeroavglabels == highestlabelneighbor) .* mask;
    
    label(logical(nodestoadd)) = highestlabelneighbor(logical(nodestoadd));
    
    
    
    end
    
 end
 disp(' ')
 
 zerosbecomebig = nodeneighlab;
 zerosbecomebig(logical(zerosbecomebig==0)) = numel(nodeneighlab)+1;
 
 watershedzonenodes = (max(nodeneighlab,[],3) ~= min(zerosbecomebig,[],3) .* (max(nodeneighlab,[],3)>0));
 label(logical(watershedzonenodes)) = 0;

    
% nodestocheck = find(mask);
% 
%      nodenum = 1;
%  
%  for node = nodestocheck';
%     
%      if node==497
%          a=1;
%      end
%      
%      
%     zerosbecomebig = nodeneighlab(node,:,:);
%     zerosbecomebig(find(zerosbecomebig==0)) = numel(nodeneighlab)+1;
%     
%     watershedzonenodes = (max(nodeneighlab(node,:,:),[],3) ~= min(zerosbecomebig,[],3)) .* (max(nodeneighlab(node,:,:),[],3)>0);
%     
%     if nnz(watershedzonenodes)>0
%         
%         string{nodenum} = ['Finalizing borders at node ' num2str(node)];
%         if nodenum==1; fprintf('%s',string{nodenum}); else fprintf([repmat('\b',1,length(string{nodenum-1})) '%s'],string{nodenum}); end
%         
%         watershedzonematrix = zeros(size(label));
%         watershedzonematrix(node,:) = watershedzonenodes;
%         
%         label(logical(watershedzonematrix)) = 0;
%         
%         watershedzoneindices = find(watershedzonematrix);
%         neighborschanged = ismember(neighbormat,watershedzoneindices);
%         
%         nodeneighlab(logical(neighborschanged)) = 0;
%         
%         nodenum = nodenum+1;
%         
%     end
%     
%  end
%     
%     disp(' ')

 
 
 meanedgevalueinparcels = sum((repmat(edgemetric,[1 size(minimametrics,2)]) .* (label>0)),1) ./ sum((label>0),1);
 meanedgevalueinborders = sum((repmat(edgemetric,[1 size(minimametrics,2)]) .* (label==0) .* nomedialnodesmat),1) ./ sum(((label==0) .* nomedialnodesmat),1);
 [ign, bestminima] = max(meanedgevalueinborders - meanedgevalueinparcels);
 
save(gifti(label),[outputdir '/' filestem '_allwatersheds.func.gii'])
 
save(gifti(label(:,bestminima)),[outputdir '/' filestem '_optimizedwatershed.func.gii'])
save(gifti(minimametrics(:,bestminima)),[outputdir '/' filestem '_optimizedminima.func.gii'])


