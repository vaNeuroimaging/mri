function [Y E diff_all_mean] = cmdscale_mat_polar(mat,groups,dist_type,varargin)
%
% Name:cmdscale_mat_polar.m
% $Revision: 1.1 $
% $Date: 2015/08/04 15:03:50 $
% This function performs multi-dimensional scaling on input matrices and
% displays a 2-dimensional plot of the relative positions of the input
% matrices colored by group. MDS is performed using euclidean distance.
% Inputs:
% mat - node x node x subject array of matrices from all subjects
% groups - identifies which subject belongs to which group
% dist_type - type of distance calculation (e.g., euclidan or hamming)
% varargin -
%   (1) group names can be specified by a third input, e.g.
%       {'group1';'group2'} and will be displayed in a legend.
%   (2) whether to plot both
%   (3) output directory
% Outputs:
% Y - configuration matrix
% E - eigenvalues of Y*Y'
%
% TOL, 09/14
% CG edits:
%   took out figure plotting command; do this outside this script
%   changed group indexing so it can deal with nonconsecutively numbered groups

numnodes = size(mat,1);
group_vals = unique(groups); % CG - added this so we can do non consecutive groups
numgroups = length(group_vals);
colors = distinguishable_colors(numgroups);
if nargin > 3
    groupnames = varargin{1};
else
    groupnames = {};
end
if length(varargin) > 1
    plot_best = varargin{2};
else
    plot_best = 0;
end
mask = ones(numnodes);
mask = triu(mask,1);

mat_col2 = zeros(size(find(mask),1),size(mat,3))*nan;
% ALSO REORDER SUBJECTS HERE SO PAIRED ARE ALWAYS SUCCESSIVE
inds_A = find(groups == group_vals(1));
inds_B = find(groups == group_vals(2));
for s = 1:size(mat,3)
    temp = mat(:,:,s);
    if groups(s) == group_vals(1) %if in the first group
        subi = find(inds_A == s);
        mat_col2(:,subi*2-1) = temp(logical(mask));
    elseif groups(s) == group_vals(2) % if in the second group
        subi = find(inds_B == s);
        mat_col2(:,subi*2) = temp(logical(mask));
    end
    mat_col(:,s) = temp(logical(mask));
end


% Multi-dimensional scaling
D = pdist(double(mat_col'),dist_type);
[Y E] = cmdscale(D);

D2 = pdist(double(mat_col2'),dist_type);
[Y2 E2] = cmdscale(D2);

% Display result
%figure('Color','white')    % CG - assumes figure has already been made outside


%Look at which dimension carries the most information about difference
%(consistently across subjects)
A_all = Y(inds_A,:);
B_all = Y(inds_B,:);
% A_all = Y2(1:2:end,:);
% B_all = Y2(2:2:end,:);
diff_all = A_all - B_all;
diff_all_mean = mean(diff_all,1); %average over subjects
[temp max_ind] = max(abs(diff_all_mean)); %find dimension w maximal difference
if max_ind > 2
    display(['Best dimension not in first two; it is ' num2str(max_ind)]);
end

if plot_best
    [b,i] = sort(abs(diff_all_mean));
    dim1 = i(end);
    dim2 = i(end -1);
    %figure;
    %plot(abs(diff_all_mean),'k'); hold on;
    %plot(dim1,abs(diff_all_mean(dim1)),'ro');
    %plot(dim2,abs(diff_all_mean(dim2)),'go');
    %xlabel('MDS dimensions');
    %ylabel('Abs distance between groups');
else
    dim1 = 1;
    dim2 = 2;
end


coords_A = [Y(inds_A,dim1), Y(inds_A,dim2)];
coords_B = [Y(inds_B,dim1), Y(inds_B,dim2)];
%coords_A = [Y2(1:2:end,1), Y2(1:2:end,2)];
%coords_B = [Y2(2:2:end,1), Y2(2:2:end,2)];
coords_diff = coords_B - coords_A;
compass(coords_diff(:,1),coords_diff(:,2),'r');
ticks_off = {'30','60','120','150','210','240','300','330'};
for t = ticks_off
    set(findall(gcf,'String',t{1}),'String','')
end


if length(groupnames) > 0
    legend(groupnames,'FontWeight','bold','FontSize',14);
end
