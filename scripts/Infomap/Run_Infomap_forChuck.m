function Run_Infomap(rmat, dmatname, xdistance, thresholdarray, makebinary, outdir, numreps, numpools)
%Run_Infomap(rmat, dmatname, xdistance, thresholdarray, makebinary, outdir, [numreps], [numpools])
%
% Run infomap on a matrix with a given distance exclusion, at various
% density thresholds, and write the results from all thresholds into a
% single text file named "rawassn.txt". This can take a long time for large
% matrices. It will run up to eight infomaps simultaneously if the Parallel
% Computing Toolbox is installed.
%
% Inputs:
%
% rmat - a correlation matrix to be infomapped. Can be a numeric matrix or
%  a cifti file that will be loaded
% dmatname - a .mat file containing a node-to-node distance matrix
% xdistance - the distance exclusion to apply, in mm (i.e., nodes closer
%  than xdistance are not allowed to be connected)
% thresholdarray - a vector of density thresholds to apply to the matrix.
%  Infomap will be run separately for each threshold.
% makebinary - whether the matrix is binarized after thresholding. 1 = make
%  it binary; 0 = leave it weighted.
% outdir - the folder results will be written to. Will be created if it
%  doesn't exist.
% numreps - the number of infomap "repetitions" (of the outside loop of the
%  algorithm. Increasing this number linearly increases the time needed to
%  run the algorithm but also increases the reliability of the output. Omit
%  or leave empty ([]) to use the default of 100.
% numpools - an OPTIONAL scalar input specifying the number of parallel
%  pools to use. Omit or leave empty ([]) to use the default of 8.
%
%
%
%
%EMG 06/25/15

if ~exist('numpools') || isempty(numpools)
    numpools = 8;
end

if ~exist('numreps') || isempty(numreps)
    numreps = 100;
end

dlmwrite([outdir '/thresholds.txt'],thresholdarray,'delimiter',' ')

[~,~] = system(['rm ' outdir '/pajek*']);

prevstring = [];

string = ['loading correlations...'];
fprintf(string);
%fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
%prevstring = string;

if ischar(rmat)
    rmat = ft_read_cifti_mod(rmat);
    rmat = rmat.data;
end
rmat = single(rmat);

warning off

string = ['applying distance exclusion...'];
fprintf(string);
%fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
%prevstring = string;
% tic

if isnumeric(dmatname)
    dmat = dmatname;
elseif strcmp(dmatname(end-9:end),'.dconn.nii')
    dmat = ft_read_cifti_mod(dmatname); dmat = dmat.data;
else
    dmat = smartload(dmatname);
end
clear dmatname
dmat = uint8(dmat);

% apply a distance exclusion?
if ~isempty(xdistance)
    if isnumeric(xdistance) && (xdistance>=0)
        rmat(dmat < xdistance) = 0;
    else
        error('xdistance is not >=0 or is not numeric.\n');
    end
end

clear dmat
%toc


numanalyses = length(thresholdarray);
numnodes = size(rmat,1);

if ~exist(outdir)
    mkdir(outdir);
end


string = ['finding r thresholds...'];
fprintf(string);
%fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
%prevstring = string;

%Remove diagonals
for i = 1:numnodes
    rmat(i,i) = 0;
end

%tic

ind = matrix_thresholder_simple(rmat,thresholdarray(end));
%toc

if makebinary
    rmat=ceil(rmat);
end

string = ['saving pajek files...'];
fprintf(string);
%fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
%prevstring = string;

% Save out largest pajek file, i.e. top threshold
pajekfileorig = [outdir '/pajek_col' num2str(length(thresholdarray)) '.net'];
mat2pajek_byindex(rmat,ind,pajekfileorig)
numpossibleedges=(numnodes*(numnodes-1))/2;

clear rmat

% Write out other thresholds
for i = 1:numanalyses
    if i<numanalyses
        pajekfile = [outdir '/pajek_col' num2str(i) '.net'];
        edgesleft=round(thresholdarray(i)*numpossibleedges);
        numuse = edgesleft + numnodes + 2; % Number of edges plus number of nodes plus 2 lines for the headers in the pajek file
        evalc(['!head -n ' num2str(numuse) ' ' pajekfileorig ' >! ' pajekfile]);
    end
end
    
string = ['running infomap'];
fprintf(string);
%fprintf([repmat('\b',1,length(prevstring)) '%s'],string);
disp(' ')

% do analyses at each threshold
processingpool = parpool(numpools);
parfor i=1:numanalyses
    %fprintf('Thr/box %d, pass %d\n',i);
    
    %tic
    pajekfile = [ outdir '/pajek_col' num2str(i) '.net' ];
    rawclrs = run_infomap_on_pajekfile(pajekfile,numreps);
    dlmwrite([outdir '/rawassn_col' num2str(i) '.txt'],rawclrs,'\t')
    %toc
end
delete(processingpool)

for i = 1:numanalyses
    pajekfile = [ outdir '/pajek_col' num2str(i) '.net' ];
    delete(pajekfile);
    delete([pajekfile(1:end-4) '.clu']);
end



for i = 1:numanalyses
    rawclrs_all(:,i) = load([outdir '/rawassn_col' num2str(i) '.txt']);
end
% write the raw assignments as .txt
dlmwrite([outdir '/rawassn.txt'],rawclrs_all,'\t');
delete([outdir '/rawassn_col*.txt'])
    
end


