function Infomap_new(rmat, dmatname, xdistance, thresholdarray, makebinary, thresholdtype, outdir)
%Infomap_new(rmat, dmatname, xdistance, thresholdarray, makebinary, thresholdtype, outdir)
%
% Run infomap on a matrix with a given distance exclusion, at various
% thresholds, and write the results from all thresholds into a single text
% file named "rawassn.txt". This can take a long time for large matrices.
% It will run up to eight infomaps simultaneously if the Parallel Computing
% Toolbox is installed.
%
% Inputs:
%
% rmat - a correlation matrix to be infomapped. Can be a numeric matrix or
%  a cifti file that will be loaded
% dmatname - a node-to-node distance matrix saved in a .mat file
% xdistance - the distance exclusion to apply, in mm (i.e., nodes closer
%  than xdistance are not allowed to be connected)
% thresholdarray - a vector of thresholds to apply to the matrix. Infomap
%  will be run separately for each threshold.
% makebinary - whether the matrix is binarized after thresholding. 1 = make
%  it binary; 0 = leave it weighted.
% thresholdtype - specify 'kden' for density thresholding or 'r' for
%  correlation thresholding
% outdir - the folder results will be written to. Will be created if it
%  doesn't exist.
%
%EMG 06/25/15


if ischar(rmat)
    rmat = ft_read_cifti_mod(rmat);
    rmat = rmat.data;
end

warning off
disp('Applying distance threshold')
tic
dmat = smartload(dmatname);

% apply a distance exclusion?
if ~isempty(xdistance)
    if isnumeric(xdistance) && (xdistance>=0)
        rmat(dmat < xdistance) = 0;
    else
        error('xdistance is not >=0 or is not numeric.\n');
    end
end

clear dmat
toc


numanalyses = length(thresholdarray);

if ~exist(outdir)
    mkdir(outdir);
end



fprintf('Form matrix and find r thresholds');
disp(' ')
%Remove diagonals
for i = 1:size(rmat,1)
    rmat(i,i) = 0;
end

tic
ind = matrix_thresholder_faster2(rmat,thresholdarray,thresholdtype);
toc

if makebinary
    rmat=ceil(rmat);
end

% Save out largest pajek file, i.e. top threshold
pajekfileorig = [ outdir '/pajek_col' num2str(length(thresholdarray)) '.net' ];
mat2pajek_mod_EG(rmat,ind,pajekfileorig)
numpossibleedges=(size(rmat,1)*(size(rmat,1)-1))/2;


% Write out other thresholds
for i = 1:numanalyses
    if i<numanalyses
        pajekfile = [outdir '/pajek_col' num2str(i) '.net'];
        edgesleft=round(thresholdarray(i)*numpossibleedges);
        numuse = edgesleft + size(rmat,1) + 2; % Number of edges plus number of nodes plus 2 lines for the headers in the pajek file
        evalc(['!head -n ' num2str(numuse) ' ' pajekfileorig ' >! ' pajekfile]);
    end
end
    


% do analyses at each threshold
matlabpool open 8
parfor i=1:numanalyses
    fprintf('Thr/box %d, pass %d\n',i);
    
    tic
    pajekfile = [ outdir '/pajek_col' num2str(i) '.net' ];
    rawclrs = run_infomap(pajekfile,100);
    dlmwrite([outdir '/rawassn_col' num2str(i) '.txt'],rawclrs,'\t')
    toc
end
matlabpool close

for i = 1:numanalyses
    pajekfile = [ outdir '/pajek_col' num2str(i) '.net' ];
    delete(pajekfile);
    delete([pajekfile(1:end-4) '.clu']);
end

clear rmat

for i = 1:numanalyses
    rawclrs_all(:,i) = load([outdir '/rawassn_col' num2str(i) '.txt']);
end
% write the raw assignments as .txt
dlmwrite([outdir '/rawassn.txt'],rawclrs_all,'\t');
delete([outdir '/rawassn_col*.txt'])
    
end

function [Ci] = run_infomap(pajekfilename,reps)
%[Ci] = run_infomap(pajekfilename,reps)
%
%
% This script runs infomap on a pajekfile with some number of
% repetitions. It then  returns the community assignments found.



% this will be the relevant output of infomap
[pathstr,cluname,ext] = filenamefinder(pajekfilename,'dotsout');
clufile = [ pathstr '/' cluname '.clu' ];

% obtain seed #
clear randnum;
randnum=ceil(rand*1000000);

% find out which computer we're on (infomap is compiled  - system specific
command= 'uname -m';
systemversion=evalc(['!' command]);
while length(systemversion)<1
    systemversion=evalc(['!' command]);
end
systemversion(end)=[]; % remove that carriage return
switch systemversion
    case 'i686'
        infomapcommand='infomap_i686';
    case 'x86_64'
        infomapcommand='infomap_x86_64';
    otherwise
        fprintf('Need to compile infomap yourself and add it to the infomap_wrapper switch possibilities.\n');
end

% run infomap
c=clock;
fprintf('\t%2.0f:%2.0f:%2.0f: infomap beginning\n',c(4),c(5),c(6));
command = ['!/data/cn4/evan/Infomap/Infomap-0.15.7/Infomap --clu -2 -s' num2str(randnum) ' -N' num2str(reps) ' ' pajekfilename ' ' pathstr];
evalc(command);
c=clock;
fprintf('\t%2.0f:%2.0f:%2.0f: infomap finished\n',c(4),c(5),c(6));



% So parfor doesn't crap out
isclufile = exist(clufile);
while isclufile == 0
    pause(60)
    isclufile = exist(clufile);
end

Ci = textread(clufile,'%d','headerlines',1);
end



