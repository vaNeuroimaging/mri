function [outmaxnodes outmaxnodes_vals] = surface_edges_all_test_faster(fullmetric,specfile,outputdir,outputname,noalone)
%SURFACE_EDGES Compute surface gradient edges
% Performs a form of non-maxima suppression, smoothing the gradient metric
% beforehand is helpful
% Input parameters available in previous versions: highthresh,lowthresh
% Output parameters in previous versions: edgenodes

%metric=gifti(metric);
%metric=metric.cdata;

[caretdir] = fileparts(specfile);
[specheader specfiles] = read_caret_spec(strcat(specfile));

startdir = cd(caretdir);
caret('-surface-topology-neighbors', specfiles.CLOSEDtopo_file{1},...
    'node_neighbors.txt');

% TODO: do we really need a larger buffer than the defalt (4096)?
% we're only reading integers here.
bufsize=16384;

% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread('node_neighbors.txt','%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;

%Identify NaN indices
nanneigh = single(isnan(neighbors(:,7)));
ind_nonan = 1:size(fullmetric,1);
ind_nonan(logical(nanneigh))=[];
ind_nan = 1:size(fullmetric,1);
ind_nan(~logical(nanneigh))=[];

%Identify neighbors
neighbors_nonan = neighbors;
neighbors_nonan(logical(nanneigh),:) = [];
just_neighbors_nonan = neighbors_nonan(:,2:7);

neighbors_nan = neighbors;
neighbors_nan(~logical(nanneigh),:) = [];
just_neighbors_nan = neighbors_nan(:,2:6);

%Find combinations of neighbors and identify non-adacent neighbors
ind = nchoosek(1:6,2);
combo_nonan = reshape(just_neighbors_nonan(:,ind),[size(just_neighbors_nonan,1) size(ind,1) 2]);
comboneigh1_nonan = reshape(neighbors(combo_nonan(:,:,1),2:7),[size(combo_nonan,1) size(combo_nonan,2) 6]);
comboneigh2_nonan = repmat(combo_nonan(:,:,2),[1 1 6]);
adjacent_nonan = sum(comboneigh1_nonan==comboneigh2_nonan,3);
nonadj_nonan = adjacent_nonan<1;

ind = nchoosek(1:5,2);
combo_nan = reshape(just_neighbors_nan(:,ind),[size(just_neighbors_nan,1) size(ind,1) 2]);
comboneigh1_nan = reshape(neighbors(combo_nan(:,:,1),2:7),[size(combo_nan,1) size(combo_nan,2) 6]);
comboneigh2_nan = repmat(combo_nan(:,:,2),[1 1 6]);
adjacent_nan = sum(comboneigh1_nan==comboneigh2_nan,3);
nonadj_nan = adjacent_nan<1;


%maxnodes = zeros(size(fullmetric));
matlabpool open 4
parfor num = 1:size(fullmetric,2)
    num
    metric = fullmetric(:,num);
    maxnode = zeros(size(fullmetric,2),1);
    %For nodes with 6 neighbors
    metric_combo = metric(combo_nonan);
    init_compare = repmat(metric(neighbors_nonan(:,1)),[1 15 2]);
    max_combo = sum(init_compare > metric_combo,3)>1;
    max_nonadj = max_combo & nonadj_nonan;
    ismax_nonadj = single(sum(max_nonadj,2)>1);
    maxnode(ind_nonan) = ismax_nonadj;
  
    %For nodes with 5 neighbors
    metric_combo = metric(combo_nan);
    init_compare = repmat(metric(neighbors_nan(:,1)),[1 10 2]);
    max_combo = sum(init_compare > metric_combo,3)>1;
    max_nonadj = max_combo & nonadj_nan;
    ismax_nonadj = single(sum(max_nonadj,2)>1);
    maxnode(ind_nan) = ismax_nonadj;
    
    
    if noalone == 1
        alone_neigh_nonan = sum(maxnode(neighbors_nonan),2)<2;
        alone_ind_nonan = ind_nonan(alone_neigh_nonan);
        maxnode(alone_ind_nonan) = 0;
        
        alone_neigh_nan = sum(maxnode(neighbors_nan(:,1:6)),2)<2;
        alone_ind_nan = ind_nan(alone_neigh_nan);
        maxnode(alone_ind_nan) = 0;
    end
    maxnodes(:,num) = maxnode;    
    
end
matlabpool close
cd(startdir);   % output might have been specified as relative path
cd(outputdir)
outmaxnodes = mean(maxnodes,2);
%outmaxnodes_vals = mean(maxnodes_vals,2);
save(gifti(outmaxnodes),[outputname '.func.gii']);
