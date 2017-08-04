function mat2pajek_mod_TL(mat,ind,roifile,outputname)
%
% Name:mat2pajek.m
% $Revision: 1.3 $
% $Date: 2011/03/08 20:30:07 $
%
% jdp 10/10/10
% 
% This takes a matrix, an roifile, and a name for the output, and writes
% the graph in the pajek format, which is used in Pajek (for windows only),
% and also for some compiled versions of things, like Infomap.
% 
% USAGE: mat2pajek(mat,roifile,outputname)
% USAGE: mat2pajek(mat,'rois.roi','/data/dolphins.net')
% 
% NOTE: pajek files should have .net extensions

% read in ROI names from roifile
[xyz name] = roifilereader(roifile);

% only the upper triangle
mat=triu(mat,1);

% get edges and values
[x y] = ind2sub(size(mat),ind);
z = mat(ind);

use = [x y z];
%%% make the input file %%%
[nodes nodes] = size(mat);
nodenum = 1:nodes;

c=clock; 
fprintf('\t%2.0f:%2.0f:%2.0f: mat2pajek: writing .net file, with %d vertices and %d edges\n',c(4),c(5),c(6),nodes,length(x));

fid=fopen(outputname,'W');
fprintf(fid,'*Vertices %d\n',nodes);
%for j=1:nodes
%    fprintf(fid,'%d "%s"\n',j,name{j,1});
%end
fprintf(fid,'%d "%d"\n',[nodenum; nodenum]);
fprintf(fid,'*Edges %d\n',length(x));

fprintf(fid,'%d %d %f\n',use');

fclose(fid);