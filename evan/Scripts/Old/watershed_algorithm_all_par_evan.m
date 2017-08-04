function labels = watershed_algorithm_all_par(edgemetrics,minimametrics,stepnum,fracmaxh)

bufsize=16384;
%caredir = '/data/cn4/laumannt/assignment_problem_v2/caret2/PALS_B12.BOTH-HEMS.CLEAN.73730/LEFT/18K';
caretdir = '/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k';

% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread([caretdir '/node_neighbors.txt'],'%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;

%Label initial markers with unique value
labels = zeros(size(minimametrics));
    
minh = min(edgemetrics(:));
maxh = max(edgemetrics(:));

stoph = maxh*fracmaxh;
step = (maxh-minh)/stepnum;
hiter = minh:step:stoph;
system('mkdir countingDir')
system('touch countingDir/countfile')
matlabpool open 8
numlabels = size(labels,2);
for j = 1:2
    
    if j == 1
        labelstodo = 1:floor(numlabels/2);
    elseif j == 2
        labelstodo = (floor(numlabels/2)+1):numlabels;
    end
    parfor l = labelstodo%size(labels,2)
        system(['echo ' num2str(l) ' >> countingDir/countfile']);
        [ign count] = system('cat countingDir/countfile | wc -l');
        disp(['Calculating edges on metric #:' num2str(l) ', Count = ' num2str(count)])
        label = labels(:,l);
        edgemetric = edgemetrics(:,l);
        
        labelpos = find(minimametrics(:,l)==1);
        randval = randn(length(labelpos),1);
        [randign randind] = sort(randval);
        temp = 1:length(labelpos);
        labelnums = temp(randind);
        label(labelpos) = labelnums;
        watershed_zones = zeros(size(label));
        for i = 1:length(hiter);
            %disp(['Number of iterations will be ' num2str(length(hiter)) ', Iteration = ' num2str(i)])
            %maskpos = find(edgemetric<sortedge(i)); % Take values in metric less than current iteration
            maskmetrics = edgemetric<hiter(i); % Take values in metric less than current iteration
            maskmetrics = maskmetrics & ~label>0 & ~watershed_zones;
            
            maskpos = find(sum(maskmetrics,2)>0);
            
            randnum = randperm(length(maskpos));
            maskpos = maskpos(randnum);
            
            for m = 1:length(maskpos) %For all nodes at this threshold
                
                nodeneigh = neighbors(maskpos(m),2:7);
                maskinthismetric = maskmetrics(maskpos(m),:);
                %nodeneigh = neighbors(sortedgepos(i),2:7);
                
                nodeneigh(isnan(nodeneigh)) = [];
                nodeneighlab = label(nodeneigh,:);
                
                %Find minimum value other than 0 among neighbors
                minfindnodeneighlab = nodeneighlab;
                minfindnodeneighlab(nodeneighlab==0) = 100000;
                minnodeneighlab = min(minfindnodeneighlab,[],1);
                
                %Find maximum value other than 0 among neighbors
                maxfindnodeneighlab = nodeneighlab;
                maxfindnodeneighlab(nodeneighlab==0) = -100000;
                maxnodeneighlab = max(maxfindnodeneighlab,[],1);
                
                %If min and max differ (i.e. two or more neighbor water), watershed
                %zone
                watershed_nodes = (minnodeneighlab~=maxnodeneighlab) & (minnodeneighlab~=100000) & maskinthismetric;
                watershed_zones(maskpos(m),watershed_nodes) = 1;
                label(maskpos(m),watershed_nodes) = 0;
                
                %If min and max the same but different from 0, add to neighbor
                %water
                next_to_water = (minnodeneighlab==maxnodeneighlab) & (minnodeneighlab~=100000) & maskinthismetric;
                label(maskpos(m),next_to_water) = minnodeneighlab(next_to_water);
                
                
            end
            
        end
        labels(:,l) = label;
    end
end
system('rm -r countingDir')
matlabpool close

%   
%  save(gifti(label),[outputdir '/' filestem 'watershed.gii'])
%  system(['mv ' outputdir '/' filestem 'watershed.gii ' outputdir '/' filestem 'watershed.func.gii'])