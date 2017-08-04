function graphcluster_Evan_surfacexd_confidence(prmfile,analysistype,xdistance,xdistancemat,savematfile,clustertype,std_water_corrmat,varargin)
%
% Name:graphcluster.m
% $Revision: 1.2 $
% $Date: 2011/03/08 20:29:54 $
%
% jdp 10/10/10
% 
% This script applies community detection algorithms to graphs. It can do
% mean matrices, or bootstrap analyses. Currently, the spectral modularity
% optimization algorithm of Newman 2006, and Infomap of Rosvall & Bergstrom
% 2008 are implemented. Look them up to understand them.
% 
% USAGE: graphcluster(prmfile,analysistype,xdistance,makebinary,savematfile,clustertype,*numbootstraps,bootstrapsamplesize*)
% USAGE: graphcluster('modbox.prm','thr',[],0,0,'modularity')
% USAGE: graphcluster('modbox.prm','thr',[],1,0,'infomap')
% USAGE: graphcluster('modbox.prm','thr',25,1,0,'infomap',1000,25)



%clf;

% read in the prmfile settings
[matfile roifile stemname subjectA subjectZ loend step hiend writepath threshold boxcarsize boxcarstep thresholdtype] = prmfilereader(prmfile);
if ~exist(writepath)
    mkdir(writepath);
end

% create filebase for output
[filestem] = filenameprep(subjectA,subjectZ,loend,step,hiend,threshold,boxcarsize,boxcarstep,analysistype,thresholdtype);

% determine the # thresholds, etc. for the analysis
[subjectarray thresholdarray numanalyses xarray] = matrix_parameter_setter(subjectA,subjectZ,loend,step,hiend,threshold,boxcarsize,boxcarstep,analysistype);

% load the matrix
[matrix nodes subjects] = matfile_loader(matfile);
clear matrix

if nodes<=10000
    issmallenough=1;
else
    issmallenough=0;
    
end

% apply a distance exclusion?
if ~isempty(xdistancemat)
    
    xd=dotremover(num2str(xdistance));
        filestem=[filestem '_surfxd' xd ];
        
end


% set variable dimensions according to bootstrapping needs
finaldim=1;
bootstrapping=0;
if ~isempty(varargin)
    bootstrapping=1;
    finaldim=varargin{1,1};
    if finaldim<1
        error('numbootstraps should be >=1');
    end
    bootsamplesize=varargin{1,1};
    if bootsamplesize<1
        error('bootsamplesize should be >=1');
    end
    filestem=[filestem '_boot' num2str(bootsamplesize) 'x' num2str(finaldim) ];
end


% create an output directory if needed
filestem=[stemname '_' filestem];
switch clustertype
    case 'modularity'
        outdir=[writepath '/' filestem '_MDLRTY' ];
    case 'infomap'
        outdir=[writepath '/' filestem '_INFMAP' ];
    otherwise
        error('Need to use ''modularity'' or ''infomap'' as clustering algorithms');
end
if ~exist(outdir)
    mkdir(outdir);
end

% check savematfile
switch savematfile
    case 0
    case 1
    otherwise
        error('savematfile should be ''0'' (no .mat saved) or ''1'' (.mat saved) ');
end

% initialize variables
rawclrs=zeros(nodes,numanalyses,finaldim,'single');
pattern=rawclrs;

Q=zeros(numanalyses,finaldim,'single');
foundr=Q;
foundkden=Q;

% do analyses at each threshold/boxcar
for i=1:numanalyses
    
    
    fprintf('Threshold %d, unaltered',i);
    
    [rmat nodes subjects] = matfile_loader(matfile);
    if ~isempty(xdistancemat)
        rmat = rmat .* xdistancemat;
    end
    
   
    
    % threshold the matrix
    [rmat foundr(i) foundkden(i)] = matrix_thresholder_Evan(rmat,thresholdarray(i,1),thresholdtype);
    
    
    
    % calculate the clusters for this network and store Q and partition
    switch clustertype
        case 'modularity'
            [realclrs(:,i) Q(i)]=olaf_modularity_und(rmat);
        case 'infomap'
            pajekfile = [ outdir '/' filestem '_col' num2str(i) '.net' ];
            %[rawclrs(:,i,j)] = infomap_wrapper(roifile,rmat,pajekfile,100,1);
            realclrs(:,i) = infomap_wrapper(roifile,rmat,pajekfile,100,1);
            %                 if issmallenough
            %                     [Q(i,j)] = M_calc_modularity(rawclrs(:,i,j),rmat);
            %                 end
    end
    
    realcolors = unique(realclrs(:,i)); realcolors(realcolors<1) = [];
    realsizes = zeros(length(realcolors),1);
    for clrnum=1:length(realcolors)
        realsizes(clrnum) = nnz(realclrs(:,i)==realcolors(clrnum));
    end
    [ign sorti] = sort(realsizes,1,'descend');
    sortrealcols = realcolors(sorti);
    
    
    randaddmat = zeros([size(rmat,1) size(rmat,2) 10000]);
    for nodei = 1:size(rmat,1);
        for nodej = 1:size(rmat,2)
            randaddmat(nodei,nodej,:) = normrnd(0,std_water_corrmat(nodei,nodej),10000,1);
        end
    end
    
    
    bootstrapclrs = zeros(size(realclrs,1),finaldim);
    for j=1:finaldim
        fprintf('Threshold %d, Boostrap %d',i,j);
        altered_rmat = (rmat + randaddmat(:,:,randi(10000,1))) .* logical(rmat);
        
        temppajekfile = [ outdir '/' filestem '_col' num2str(i) '_temp.net' ];
        iteration_clrs = infomap_wrapper(roifile,altered_rmat,temppajekfile,100,1);
        evalc('fclose all');
        
        regularized_iteration_clrs = zeros(size(iteration_clrs));
        
        iterationcolors = unique(iteration_clrs); iterationcolors(iterationcolors<1) = [];
        itercolorsavailable = iterationcolors;
        
        for clri = 1:length(sortrealcols)
            
            if ~isempty(itercolorsavailable)
                
                itersizes = zeros(length(itercolorsavailable),1);
                for clrj = 1:length(itercolorsavailable)
                    itersizes(clrj) = nnz(iteration_clrs(realclrs(:,i)==sortrealcols(clri))==itercolorsavailable(clrj));
                end
                [maxval maxi] = max(itersizes);
                if maxval>0
                    regularized_iteration_clrs(iteration_clrs==itercolorsavailable(maxi)) = sortrealcols(clri);
                    itercolorsavailable(maxi) = [];
                end
            end
        end
        for k = 1:length(itercolorsavailable)
            regularized_iteration_clrs(iteration_clrs==itercolorsavailable(k)) = 1000+k;
        end
        
        
%         previous_splitoffs = bootstrapclrs(:,1:(j-1)); previous_splitoffs(previous_splitoffs<1000) = 0;
%         previous_splitoffs = mode(previous_splitoffs,2);
%         previous_splitoff_colors = unique(previous_splitoffs(previous_splitoffs>1000));
%        
%         if (~isempty(previous_splitoff_colors)) && any(regularized_iteration_clrs>1000)
%             
%             maxnum = max(previous_splitoff_colors);
%             
%             thisiter_splitoff_colors = unique(regularized_iteration_clrs(regularized_iteration_clrs>1000));
%             clear thisiter_splitoff_color_sizes
%             for k = 1:length(thisiter_splitoff_colors)
%                 thisiter_splitoff_color_sizes(k,1) = nnz(regularized_iteration_clrs==thisiter_splitoff_colors(k));
%             end
%             [ign sorti] = sort(thisiter_splitoff_color_sizes,1,'descend');
%             thisiter_splitoff_colors_sorted = thisiter_splitoff_colors(sorti);
%             
%             
%             for thiscolor = thisiter_splitoff_colors_sorted(:)'
%                 clear numoverlap
%                 for l = 1:length(previous_splitoff_colors)
%                     numoverlap(l) = nnz((previous_splitoffs==previous_splitoff_colors(l)) .* (regularized_iteration_clrs==thiscolor));
%                 end
%                 [maxoverlap maxi] = max(numoverlap);
%                 if ((maxoverlap / nnz((previous_splitoffs==previous_splitoff_colors(maxi)))) > .5) && ((maxoverlap / nnz(regularized_iteration_clrs==thiscolor)) > .5)
%                     regularized_iteration_clrs(regularized_iteration_clrs==thiscolor) = previous_splitoff_colors(maxi);
%                     previous_splitoff_colors(maxi) = [];
%                 else
%                     maxnum = maxnum+1;
%                     regularized_iteration_clrs(regularized_iteration_clrs==thiscolor) = maxnum;
%                 end
%                 
%             end
%         end
        bootstrapclrs(:,j) = regularized_iteration_clrs;
        
    end
    
    for node = 1:size(realclrs,1)
        
        confidence(node,i) = sum(bootstrapclrs(node,:)==realclrs(node,i)) / size(bootstrapclrs,2);
        thisnode_otherids = bootstrapclrs(node,:); thisnode_otherids(thisnode_otherids==realclrs(node,i)) = []; thisnode_otherids(thisnode_otherids>1000) = [];
        secondid(node,i) = mode(thisnode_otherids);
        secondid_confidence(node,i) = sum(bootstrapclrs(node,:)==secondid(node,i)) / size(bootstrapclrs,2);
        thisnode_otherids(thisnode_otherids==secondid(node,i)) = [];
        thirdid(node,i) = mode(thisnode_otherids);
        thirdid_confidence(node,i) = sum(bootstrapclrs(node,:)==thirdid(node,i)) / size(bootstrapclrs,2);
    end
    
end

dlmwrite([outdir '/rawassn.txt'],realclrs,'\t');
dlmwrite([outdir '/confidence.txt'],confidence,'\t');
dlmwrite([outdir '/second_id.txt'],secondid,'\t');
dlmwrite([outdir '/secondid_confidence.txt'],secondid_confidence,'\t');
dlmwrite([outdir '/third_id.txt'],thirdid,'\t');
dlmwrite([outdir '/thirdid_confidence.txt'],thirdid_confidence,'\t');




