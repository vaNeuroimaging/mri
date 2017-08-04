function [g1g2_diff_mat group1_gstar group2_gstar pvalue] = bingat_glrtPvalue(graphs,groups,numbootstraps,stemname)
%
% Name:bingat_glrtPvalue_wrapper.m
% $Revision: 1.1 $
% $Date: 2015/08/03 20:32:36 $
% This is a wrapper for the bingat_glrtPvalue function run in R.
% Inputs:
% graphs - node x node x subject array of graphs from all subjects
% groups - vector of 1s and 2s identifying to which group a subject
% belongs. There are always only two groups.
% numbootstraps - number of bootstrap iterations used by glrtPvalue function
% stemname - filename prefix for output files
%
% Outputs:
% This function writes out a number of text files of results from the
% glrtPvalue function. It also outputs several matlab variables:
% g1g2_diff_mat - difference graph between the two groups, 1
% indicates an edge exists in group1 and not group2, -1
% indicates an edge exists in group2 and not group1.
% group1_gstar - gstar graph of group1
% group2_gstar - gstar graph of group2
% pvalue - significance of graph difference between groups
% TOL, 09/14

numnodes = size(graphs,1);
% Write out edges from upper triangle of input matrix
mask = ones(numnodes);
mask = triu(mask,1);
for s = 1:size(graphs,3)
    temp = graphs(:,:,s);
    graph_col(:,s) = temp(logical(mask));
end

dlmwrite([stemname '.txt'],graph_col,'\t')
group1_num = nnz(groups==1);
group2_num = nnz(groups==2);
sep = '\t';

% Write R script
fid = fopen([stemname '.r'],'wt');
fprintf(fid,'library(bingat)\n');
fprintf(fid,'data <- read.table("%s.txt")\n',stemname);
fprintf(fid,'groups <- c(rep(0,%f),rep(1,%f))\n',group1_num,group2_num);
fprintf(fid,'numbootstraps <- %f\n',numbootstraps);
fprintf(fid,'values <- glrtPvalue(data,"adjmatrix",groups,numbootstraps,FALSE)\n');
fprintf(fid,'write.table(values$b1.Differences,"%s_b1_diff.txt",sep="%s",col.names=FALSE,row.names=FALSE)\n',stemname,sep);
fprintf(fid,'write.table(values$b0.covs0,"%s_b0_covs0.txt",sep="%s",col.names=FALSE,row.names=FALSE)\n',stemname,sep);
fprintf(fid,'write.table(values$b0b1.covs1,"%s_b0b1_covs1.txt",sep="%s",col.names=FALSE,row.names=FALSE)\n',stemname,sep);
fprintf(fid,'write.table(values$pvalue,"%s_pvalue.txt",sep="%s",col.names=FALSE,row.names=FALSE)\n',stemname,sep);

% Run R script
evalc(['!chmod u+wrx ' stemname '.r']);
evalc(['!R CMD BATCH ' stemname '.r']);

% So parfor doesn't crap out
ispvalfile = exist([stemname '_pvalue.txt']);
while ispvalfile == 0
    pause(10)
    ispvalfile = exist([stemname '_pvalue.txt']);
end

mask = ones(numnodes);
mask = triu(mask,1);


b1_diff = load([stemname '_b1_diff.txt']);
matrix=zeros(numnodes);
matrix(logical(mask)) = b1_diff;
b1_diff_mat = matrix+matrix'; % Re-symmetrify matrix


b0_covs0 = load([stemname '_b0_covs0.txt']);
matrix=zeros(numnodes);
matrix(logical(mask)) = b0_covs0;
group1_gstar = matrix+matrix'; % Re-symmetrify matrix

b0b1_covs1 = load([stemname '_b0b1_covs1.txt']);
matrix=zeros(numnodes);
matrix(logical(mask)) = b0b1_covs1;
group2_gstar = matrix+matrix'; % Re-symmetrify matrix

pvalue = load([stemname '_pvalue.txt']);

g1g2_diff_mat = group1_gstar-group2_gstar; % Create difference of groups graph

evalc(['!rm ' stemname '.r']);
