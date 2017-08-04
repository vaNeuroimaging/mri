tmaskfile = '/data/cn5/selfRegulation/V4Process_nosmooth/Finaltmasklist_120_108.txt';
[subjects tmasks] = textread(tmaskfile,'%s %s');

subjectmatches = cifti_read('Templatematch_dice_bysubject.dtseries.nii');

networkIDs = unique(subjectmatches); networkIDs(networkIDs<1) = [];

networkareas = zeros(length(subjects),length(networkIDs));

hems = {'L','R'};

maskL = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii']); maskL = maskL.cdata;
maskR = gifti(['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.R.32k_fs_LR.func.gii']); maskR = maskR.cdata;

for s = 1:length(subjects)
    subname = subjects{s};
    disp(['Subject ' num2str(s)])
    surffolder = ['/data/cn4/segmentation/freesurfer5_supercomputer/FREESURFER_fs_LR/' subname '/7112b_fs_LR/fsaverage_LR32k/'];
    
    if ~exist([surffolder subname '.L.midthickness.32k_fs_LR_surfaceareas.func.gii'])
        
        system(['wb_command -surface-vertex-areas ' surffolder subname '.L.midthickness.32k_fs_LR.surf.gii ' surffolder subname '.L.midthickness.32k_fs_LR_surfaceareas.func.gii']);
        system(['wb_command -surface-vertex-areas ' surffolder subname '.R.midthickness.32k_fs_LR.surf.gii ' surffolder subname '.R.midthickness.32k_fs_LR_surfaceareas.func.gii']);
        
    end
    
    SA_L = gifti([surffolder subname '.L.midthickness.32k_fs_LR_surfaceareas.func.gii']);
    SA_R = gifti([surffolder subname '.R.midthickness.32k_fs_LR_surfaceareas.func.gii']);
    
    SA_cifti = [SA_L.cdata(logical(maskL==0)) ; SA_R.cdata(logical(maskR==0))];
    
    for IDnum = 1:length(networkIDs)
        
        ID = networkIDs(IDnum);
        
        networkareas(s,IDnum) = sum(SA_cifti(find(subjectmatches(:,s)==ID)));
        
    end
end

%%

Y = pdist(networkareas,'correlation');           % 'pdist' converts the square adjacency matrix to a
%  1 x n matrix so that the function linkage can construct the tree

clustering = linkage(Y, 'average');      % 'linkage' computes the data to construct the tree
% 'average' refers to the UPGMA algorithm


%clusters = cluster(clustering, 'MaxClust', [1:200]);
% for numclust = 1:size(clusters,2)
%     Qvals(numclust) = M_calc_modularity(clusters(:,numclust),subject_similaritymat);
% end
% [maxQval maxQind] = max(Qvals);
% 
% cut = clustering(end-maxQind,3);
cut = .18;

colormap = [1 0 0;0 0 1;1 1 0; 0 .8 0; 1 .6 1;1 .5 0;1 .7 .4;0 .6 .6;.6 .2 1];

[H,T,perm] = dendrogram(clustering, 0, 'orientation','left', 'colorthreshold', cut+.00001,'ORIENTATION','left');

%[H,T,perm] = dendrogram_evan(clustering, colormap,sizethresh,0, 'orientation','left', 'colorthreshold', cut+.00001 ,'ORIENTATION','left');



cophenetic_r = cophenet(clustering, Y);
disp(cophenetic_r)
    