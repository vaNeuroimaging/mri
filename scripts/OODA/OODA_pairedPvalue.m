function [g1g2_diff_mat group1_gstar group2_gstar pvalue] = OODA_pairedPvalue(graphs,groups,numperms,stemname)
%
% Name:OODA_pairedPvalue_wrapper.m
% $Revision: 1.1 $
% $Date: 2015/08/05 14:56:26 $
% This is a wrapper for the bingat_pairedPvalue function run in R.
% Inputs:
% graphs - node x node x subject array of graphs from all subjects
% groups - vector of 1s and 2s identifying to which group a subject
% belongs. There are always only two groups. 
%   numbers don't need to be 1 and 3 group1 will always be lower number and group2 will be higher - CG]
%   groups MUST be the same size and the position of the paired subject
%   MUST be the same in the graph order (e.g., 1st of group1 will
%   correspond to first of group2)
% numperms - number of permutation iterations to be used
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
% CG, 9/14
%
% adapted from TOL script: bingat_glrtPvalue_wrapper.m
% CG - must be run on Typhoon!

homedir = '/data/cn/data1/scripts/OODA';
numnodes = size(graphs,1);
% Write out edges from upper triangle of input matrix
mask = ones(numnodes);
mask = triu(mask,1);
for s = 1:size(graphs,3)
    temp = graphs(:,:,s);
    graph_col(:,s) = temp(logical(mask));
end

dlmwrite([stemname '.txt'],graph_col,'\t');

% number of subjects in each group (assumed to be tot/2)
nsubXgroup = size(graphs,3)/2;
if rem(size(graphs,3),2) ~= 0
    error('Must have an even number of subjects, paired in each group.')
end

%group1_num = nnz(groups==1); %CG
%group2_num = nnz(groups==2); %CG
%CG modified slightly so I don't need to use only groups 1 and 2 as labels
group_ids = unique(groups);
if length(group_ids) > 2 || length(group_ids) < 2
    error('incorrect number of groups');
end
group1_num = nnz(groups == group_ids(1));
group2_num = nnz(groups == group_ids(2));
sep = '\t';

% Write R script
fid = fopen([stemname '.r'],'wt');
fprintf(fid,'library(bingat)\n');
fprintf(fid,'source("%s/OODA_paired_permutation.R")\n',homedir);
fprintf(fid,'data <- read.table("%s.txt")\n',stemname);
fprintf(fid,'groups <- c(rep(0,%f),rep(1,%f))\n',group1_num,group2_num);
fprintf(fid,'numperms <- %f\n',numperms);
fprintf(fid,'groupsize <- %f\n',nsubXgroup);
%fprintf(fid,'gstarG1 <- estGStar(data[,1:%f])\n',nsubXgroup);
%fprintf(fid,'gstarG2 <- estGStar(data[,%f+1:2*%f])\n',nsubXgroup,nsubXgroup);
fprintf(fid,'bmp("%s_mdsplot.bmp")\n',stemname);
fprintf(fid,'values <- permutationPvalue(data,numperms,groupsize,plot=TRUE)\n');
fprintf(fid,'dev.off()\n');
fprintf(fid,'write.table(values$gstarG1,"%s_gstarG1.txt",sep="%s",col.names=FALSE,row.names=FALSE)\n',stemname,sep);
fprintf(fid,'write.table(values$gstarG2,"%s_gstarG2.txt",sep="%s",col.names=FALSE,row.names=FALSE)\n',stemname,sep);
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


G1 = load([stemname '_gstarG1.txt']);
matrix=zeros(numnodes);
matrix(logical(mask)) = G1;
group1_gstar = matrix+matrix'; % Re-symmetrify matrix

G2 = load([stemname '_gstarG2.txt']);
matrix=zeros(numnodes);
matrix(logical(mask)) = G2;
group2_gstar = matrix+matrix'; % Re-symmetrify matrix

pvalue = load([stemname '_pvalue.txt']);

g1g2_diff_mat = group1_gstar-group2_gstar; % Create difference of groups graph

evalc(['!rm ' stemname '.r'])
