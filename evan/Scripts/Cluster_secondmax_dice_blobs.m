subjects = '/home/data/subjects/DART.txt';
analysisfolder = '/home/data/Analysis/Hub_clustering/Secondmax_Dice_TM/';

secondmaxthresh = .15;
minsize = 20;

if iscell(subjects)

elseif exist(subjects,'file')
    
    subjects = textread(subjects,'%s');
    
elseif ischar(subjects)
    
    subjects = {subjects};
    
end


for s = 1:length(subjects)
    subject = subjects{s};
    disp(subject)
    folder = ['/home/data/subjects/' subject '/template_matching/'];
    
    if s==1
        %neighbors = cifti_neighbors([folder '/Templatematch_dice_bysubject_kden0.05_dicemaxes.dtseries.nii']);
        secondmax_out_allsubs = zeros(size(neighbors,1),1);
    end
    
%     maxes = ft_read_cifti_mod([folder '/Templatematch_dice_bysubject_kden0.05_dicemaxes.dtseries.nii']);
%     diffs = ft_read_cifti_mod([folder '/Templatematch_dice_bysubject_kden0.05_dicediffs.dtseries.nii']);
%     secondmax = maxes.data - diffs.data;
%     secondmax_thresholded = secondmax > secondmaxthresh;
%     
%     secondmax_thresholded_cutlines = zeros(size(secondmax_thresholded));
%     inds = find(secondmax_thresholded);
%     for i = inds(:)'
%         vertneighs = neighbors(i,2:end); vertneighs(isnan(vertneighs)) = [];
%         vertneighs(secondmax_thresholded(vertneighs)==0) = [];
%         for j = vertneighs(:)'
%             neighofneighs = neighbors(j,2:end); neighofneighs(isnan(neighofneighs)) = [];
%             neighofneighs = intersect(neighofneighs,vertneighs);
%             if any(secondmax_thresholded(neighofneighs))
%                 secondmax_thresholded_cutlines(i) = 1;
%                 break
%             end
%         end
%     end
%     
%     
%     
%     
%     data_inthresh = find(secondmax_thresholded_cutlines);
% 
%     %initialize the metric keeping track of unique cluster identifiers
%     clustereddata = zeros(size(secondmax_thresholded_cutlines));
%     
%     for vertex = data_inthresh'
%         
%         %find the neighbors of this vertex
%         vertexneighbors = neighbors(vertex,:);
%         
%         %find which of those neighbors also pass the thresholds
%         vertexneighbors_inthresh = intersect(data_inthresh,vertexneighbors);
%         
%         %find if those neighbors have already been assigned different cluster values
%         uniqueneighborvals = unique(clustereddata(vertexneighbors_inthresh));
%         uniqueneighborvals(uniqueneighborvals==0) = [];
%         
%         %if no neighbors have cluster identifiers, assign them the number of this vertex as a unique cluster identifier
%         if isempty(uniqueneighborvals)
%             clustereddata(vertexneighbors_inthresh) = vertex;
%             %if there is only one previous cluster identifier present, make all the neighbors that value
%         elseif length(uniqueneighborvals)==1
%             clustereddata(vertexneighbors_inthresh) = uniqueneighborvals;
%             %if there are multiple cluster identifier values in the neighborhood, merge them into one
%         else
%             for valuenum = 2:length(uniqueneighborvals)
%                 clustereddata(clustereddata==uniqueneighborvals(valuenum)) = uniqueneighborvals(1);
%             end
%         end
%         
%     end
%     
%     %find out what the unique cluster identifier values are
%     uniqueclustervals = unique(clustereddata);
%     uniqueclustervals(uniqueclustervals==0) = [];
%     
%     outputclusters = zeros(size(clustereddata,1),1);
%     
%     %save each unique cluster that passes the cluster size minimum into a column of the output metric
%     clustercount = 1;
%     
%     for clusternum = 1:length(uniqueclustervals)
%         
%         if nnz(clustereddata==uniqueclustervals(clusternum)) > minsize
%             outputclusters(clustereddata == uniqueclustervals(clusternum)) = clustercount;
%             clustercount = clustercount + 1;
%         end
%         
%     end
%     
%     
%     
%     maxes.data = logical(outputclusters);
%     ft_write_cifti_mod([folder '/SecondhighestDice_thresh' num2str(secondmaxthresh) '_' num2str(minsize) 'verts'],maxes)
    %secondmax_out_allsubs = secondmax_out_allsubs + (single(logical(outputclusters)) / length(subjects));
    
    outputclusters = ft_read_cifti_mod([folder '/SecondhighestDice_thresh' num2str(secondmaxthresh) '_' num2str(minsize) 'verts.dtseries.nii']);
    secondmax_out_allsubs = secondmax_out_allsubs + (double(outputclusters.data) / length(subjects));
    
    
    
end

maxes.data = secondmax_out_allsubs;
ft_write_cifti_mod([analysisfolder '/SecondhighestDice_thresh' num2str(secondmaxthresh) '_' num2str(minsize) 'verts_allsubsprobmap'],maxes)