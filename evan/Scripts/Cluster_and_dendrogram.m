


Y = pdist(networkconnections','hamming');           % 'pdist' converts the square adjacency matrix to a
%  1 x n matrix so that the function linkage can construct the tree

clustering = linkage(squareform(Y), 'average','hamming');      % 'linkage' computes the data to construct the tree
% 'average' refers to the UPGMA algorithm

clusters = cluster(clustering, 'MaxClust', [1:80]);

for s = 1:size(networkconnections,2)
    for s2 = 1:size(networkconnections,2)
        similaritymat(s,s2) = nnz(networkconnections(:,s)==networkconnections(:,s2));
    end
end


for numclust = 1:size(clusters,2)
    Qvals(numclust) = M_calc_modularity(clusters(:,numclust),similaritymat);
end
[maxQval maxQind] = max(Qvals);

cut = clustering(end-maxQind,3);


colormap = [1 0 0;0 0 1;1 1 0; 0 .8 0; 1 .6 1;1 .5 0;1 .7 .4;0 .6 .6;.6 .2 1];
    %0 0 0;.2 1 1;;0 .2 .4;.2 1 .2;0 0 1;1 1 .8;0 .4 0;.25 .25 .25];

[H,T,perm] = dendrogram_evan(clustering, colormap,sizethresh,0, 'orientation','left', 'colorthreshold', cut+.00001 ,'ORIENTATION','left');
        
%[H,T,perm] = dendrogram(clustering, 0, 'orientation','left', 'colorthreshold', cut,'ORIENTATION','left');

cophenetic_r = cophenet(clustering, Y);
disp(cophenetic_r)