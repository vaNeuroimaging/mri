function graphcluster_wdistmat_new_mod(prmfile,analysistype,xdistance,makebinary,savematfile,clustertype,varargin)
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
% TOL modifed, 08/2014



clf;

% read in the prmfile settings
[matfile roifile dmatfile stemname subjectA subjectZ loend step hiend writepath threshold boxcarsize boxcarstep thresholdtype] = prmfilereader_wdistmat(prmfile);
if ~exist(writepath)
    mkdir(writepath);
end

% create filebase for output
[filestem] = filenameprep(subjectA,subjectZ,loend,step,hiend,threshold,boxcarsize,boxcarstep,analysistype,thresholdtype);

% determine the # thresholds, etc. for the analysis
[subjectarray thresholdarray numanalyses xarray] = matrix_parameter_setter(subjectA,subjectZ,loend,step,hiend,threshold,boxcarsize,boxcarstep,analysistype);

% load the matrix
[matrix nodes subjects] = matfile_loader(matfile);

if nodes<=10000
    issmallenough=1;
else
    issmallenough=0;
end

% apply a distance exclusion?
if ~isempty(xdistance)
    if isnumeric(xdistance) && (xdistance>=0)
        dmat = load(dmatfile);
        dmat = dmat.distmat_use;
        %dmat = dmat.distmat;
        [matrix] = matrix_xdistance_wdistmat(matrix,dmat,xdistance);
        xd=dotremover(num2str(xdistance));
        filestem=[filestem '_xd' xd ];
    else
        error('xdistance is not >=0 or is not numeric.\n');
    end
end

% clear matrix % For efficiency
% check makebinary
switch makebinary
    case 0
    case 1
        filestem=[filestem '_BI' ];
    otherwise
        error('makebinary should be ''0'' or ''1'' ');
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

fprintf('Form matrix and find r thresholds');
disp(' ')
[rmat] = matrix_former_TL(matrix,subjectarray(1,2),subjectarray(1,2),'2D','diagout');
tic
[ind foundr foundkden] = matrix_thresholder_faster2(rmat,thresholdarray,thresholdtype);
toc

if makebinary
    rmat=ceil(rmat);
end

% Save out largest pajek file, i.e. top threshold
pajekfileorig = [ outdir '/' filestem '_col' num2str(size(thresholdarray,1)) '.net' ];
mat2pajek_mod_TL(rmat,ind,roifile,pajekfileorig)
numpossibleedges=(size(rmat,1)*(size(rmat,1)-1))/2;

divisions = 3;
threshsperdivision = ceil(numanalyses/divisions);

for division = 1:divisions
thisdivision_threshs = [((threshsperdivision * (division-1)) +1) : min((threshsperdivision*division),numanalyses)];
    
% Write out other thresholds
for i = thisdivision_threshs
    if i<numanalyses
    pajekfile = [outdir '/' filestem '_col' num2str(i) '.net'];
    edgesleft=round(thresholdarray(i)*numpossibleedges);
    numuse = edgesleft + size(rmat,1) + 2; % Number of edges plus number of nodes plus 2 lines for the headers in the pajek file
    evalc(['!head -n ' num2str(numuse) ' ' pajekfileorig ' >! ' pajekfile]);
    end
end
    


% do analyses at each threshold/boxcar
matlabpool open 10
parfor i=thisdivision_threshs
   % for j=1:finaldim
        fprintf('Thr/box %d, pass %d\n',i);
        
   tic     
%         rmat(rmat<foundr(i)) = 0;
%         edgenum(i) = sum(sum(rmat>0,2),1);
%         calc_kden(i) = edgenum(i)/numpossibleedges;
        % if user wants binarized networks

        
        % calculate the clusters for this network and store Q and partition
        switch clustertype
          %  case 'modularity'
          %     [rawclrs(:,i,j) Q(i,j)]=olaf_modularity_und(rmat);
            case 'infomap'
                pajekfile = [ outdir '/' filestem '_col' num2str(i) '.net' ];
                rawclrs = infomap_wrapper_mod_TL(roifile,pajekfile,1,1);
                dlmwrite([outdir '/rawassn_col' num2str(i) '.txt'],rawclrs,'\t')
               % if issmallenough
     %           [Q(i,j)] = M_calc_modularity(rawclrs(:,i,j),rmat);
     %           dlmwrite([outdir '/Q_rawassn_col' num2str(i) '.txt'],Q(i,j),'\t') 
               % end
        end
  %  end
    toc
end
matlabpool close

for i = thisdivision_threshs
    pajekfile = [ outdir '/' filestem '_col' num2str(i) '.net' ];
    delete(pajekfile);
    delete([pajekfile(1:end-4) '.clu']);
end

end
clear rmat

% for simple analyses without any bootstrapping
if finaldim==1
    for i = 1:numanalyses
        rawclrs_all(:,i) = load([outdir '/rawassn_col' num2str(i) '.txt']);
    end
    % write the raw assignments as .txt and .tiff
    dlmwrite([outdir '/rawassn.txt'],rawclrs_all,'\t');
    delete([outdir '/rawassn_col*.txt'])
   % tiffmaker([outdir '/rawassn.tiff'],rawclrs);
    
    % write the thr/box and Q values out
    if issmallenough
        fid=fopen([outdir '/thrbox_r_kden.txt'],'w'); fprintf(fid,'Thr/box\tfoundr\tfoundkden\n'); fclose(fid);
        dlmwrite([outdir '/thrbox_r_kden.txt'],[xarray foundr' foundkden'],'delimiter','\t','-append');
 %       plot(xarray,foundr,'r.',xarray,foundkden,'b.',xarray,Q,'g.'); xlabel('Thr/box'); ylabel('red:threshold blue:kden green:Q value');
 %       saveas(gcf,[outdir '/thrbox_r_kden_Q.tiff'],'tiff');
    else
    %    fid=fopen([outdir '/thrbox_r_kden_Q.txt'],'w'); fprintf(fid,'Thr/box\tfoundr\tfoundkden\tQ\n'); fclose(fid);
    %    dlmwrite([outdir '/thrbox_r_kden_Q.txt'],[xarray foundr foundkden],'delimiter','\t','-append');
    end
    
    if issmallenough % don't attempt this with large networks
        % now sort the raw assignments into sensible patterns
        [pattern] =  rawoutput2clr(rawclrs_all);
        dlmwrite([outdir '/pattern.txt'],pattern,'\t');
  %      tiffmaker([outdir '/pattern.tiff'],pattern);
    end
    
    % save a matfile of these variables
    if savematfile
        save([outdir '/' filestem '.mat']);
    end
    
    if issmallenough
        M_visuals_wdistmat(prmfile,analysistype,[outdir '/pattern.txt'],0,'visuals',[],xdistance,makebinary,1)
    end
    
else % for analyses with bootstrapping
    
    % write the thr/box and Q values out
    fid=fopen([outdir '/thrbox_r_kden_Q.txt'],'w'); fprintf(fid,'Thr/box\tfoundr\tstdfoundr\tfoundkden\tstdfoundkden\tQ\tstdQ\n'); fclose(fid);
    dlmwrite([outdir '/thrbox_r_kden_Q.txt'],[xarray mean(foundr,2) std(foundr,[],2) mean(foundkden,2) std(foundkden,[],2) mean(Q,2) std(Q,[],2)],'delimiter','\t','-append');
    plot(xarray,mean(foundr,2),'r.',xarray,mean(foundkden,2),'b.',xarray,mean(Q,2),'g.'); xlabel('Thr/box'); ylabel('red:threshold blue:kden green:Q value');
    hold on; errorbar(xarray,mean(foundr,2),std(foundr,[],2),'r.'); hold off;
    hold on; errorbar(xarray,mean(foundkden,2),std(foundkden,[],2),'b.'); hold off;
    hold on; errorbar(xarray,mean(Q,2),std(Q,[],2),'g.'); hold off;
    saveas(gcf,[outdir '/thrbox_r_kden_Q.tiff'],'tiff');
        
    % calculate tallies of shared module assignments
    mastertally=zeros(nodes,nodes,numanalyses,'single');
    for i=1:numanalyses
        for j=1:finaldim
            temptally=zeros(nodes,nodes,'single');
            [temptally]=bootstrapper_tally(rawclrs_all(:,i,j));
            mastertally(:,:,i)=mastertally(:,:,i)+temptally;
        end
    end
    
    % normalize the mastertally matrix by the number of bootstraps
    normmastertally=mastertally/finaldim;
    save([outdir '/ntally.mat'],'normmastertally');
    
    % save a matfile of these variables
    if savematfile
        save([outdir '/' filestem '.mat']);
    end
    
    % get stuff for further analysis into the outdir folder
    copyfile(roifile,outdir);
    quickprmfile([outdir '/modbox.prm'],[ 'ntally.mat' ],roifile,[filestem],1,numanalyses,0,0.05,0.95,[outdir],.9,1,1,'r');
    
end


