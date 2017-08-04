function Infomap_new_structspecific_kden(rmatname, dmatname, xdistance, thresholdarray, makebinary, outdir,structinds)
%Infomap_new_structspecific_kden(rmatname, dmatname, xdistance, thresholdarray, makebinary, outdir,structinds)
% 
% rmatname is a gifti or cifti representing the Fisher-Transformed
%  all-to-all correlation matrix
% dmatname is a .mat file containing the distance matrix
% xdistance is the distance exclusion, in mm
% thresholdarray is an ordered vector of kden thresholds to run (smallest to largest)
% makebinary indicates whether the thresholded matrix should be binarized or not
% outdir is the directory the outputs will be written into
% structinds is a nodes X 1 array with a separate value for the indices of
%  each structure being seperately thresholded.

thresholdtype = 'kden';

if strcmp(rmatname(end-8:end),'.func.gii')
    rmat = gifti(rmatname); rmat = rmat.cdata;
elseif strcmp(rmatname(end-9:end),'.dconn.nii')
    rmat = cifti_read(rmatname);
else
    error('Needs a cifti or a gifti')
end

warning off
disp('Applying distance threshold')
tic
dmat = load(dmatname);
name = fieldnames(dmat);
name = char(name);
command = ['dmat.' name];
dmat = eval(command);

% apply a distance exclusion?
if ~isempty(xdistance)
    if isnumeric(xdistance) && (xdistance>=0)
        %dmat = load(dmatfile);
        %dmat = dmat.distmat_use;
        %dmat = dmat.distmat;
        %[rmat] = matrix_xdistance_wdistmat(rmat,dmat,xdistance);
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
%[rmat] = matrix_former_TL(matrix,subjectarray(1,2),subjectarray(1,2),'2D','diagout');
%Remove diagonals
for i = 1:size(rmat,1)
    rmat(i,i) = 0;
end

tic
ind = matrix_thresholder_faster2_structurespecific(rmat,thresholdarray,thresholdtype,structinds);
toc

if makebinary
    rmat=ceil(rmat);
end

% Save out largest pajek file, i.e. top threshold
pajekfileorig = [ outdir '/pajek_col' num2str(length(thresholdarray)) '.net' ];
mat2pajek_mod_EG(rmat,ind,pajekfileorig)
numnodes = size(rmat,1);
clear rmat ind
numpossibleedges=(numnodes*(numnodes-1))/2;

%divisions = 3;
%threshsperdivision = ceil(numanalyses/divisions);

% for division = 1:divisions
% thisdivision_threshs = [((threshsperdivision * (division-1)) +1) : min((threshsperdivision*division),numanalyses)];

% Write out other thresholds
for i = 1:numanalyses
    if i<numanalyses
        pajekfile = [outdir '/pajek_col' num2str(i) '.net'];
        edgesleft=round(thresholdarray(i)*numpossibleedges);
        numuse = edgesleft + numnodes + 2; % Number of edges plus number of nodes plus 2 lines for the headers in the pajek file
        evalc(['!head -n ' num2str(numuse) ' ' pajekfileorig ' >! ' pajekfile]);
    end
end
    


% do analyses at each threshold/boxcar
matlabpool open 10
parfor i=1:numanalyses
   % for j=1:finaldim
        fprintf('Thr/box %d, pass %d\n',i);
        
        tic
        pajekfile = [ outdir '/pajek_col' num2str(i) '.net' ];
        rawclrs = infomap_wrapper_mod_EG(pajekfile,100);
        dlmwrite([outdir '/rawassn_col' num2str(i) '.txt'],rawclrs,'\t')
        toc
end
matlabpool close

for i = 1:numanalyses
    pajekfile = [ outdir '/pajek_col' num2str(i) '.net' ];
    delete(pajekfile);
    delete([pajekfile(1:end-4) '.clu']);
end

%end
%clear rmat

    for i = 1:numanalyses
        rawclrs_all(:,i) = load([outdir '/rawassn_col' num2str(i) '.txt']);
    end
    % write the raw assignments as .txt and .tiff
    dlmwrite([outdir '/rawassn.txt'],rawclrs_all,'\t');
    delete([outdir '/rawassn_col*.txt'])
   % tiffmaker([outdir '/rawassn.tiff'],rawclrs);
    
   
    



