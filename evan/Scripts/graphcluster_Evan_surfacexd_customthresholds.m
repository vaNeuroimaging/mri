function graphcluster_Evan_surfacexd_customthresholds(prmfile,analysistype,xdistance,xdistancemat,makebinary,savematfile,clustertype,thresholdarray,varargin)
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
%[subjectarray ign ign2 xarray] = matrix_parameter_setter(subjectA,subjectZ,loend,step,hiend,threshold,boxcarsize,boxcarstep,analysistype);

numanalyses = length(thresholdarray);

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

% do analyses at each threshold/boxcar
for i=1:numanalyses
    for j=1:finaldim
        fprintf('Thr/box %d, pass %d\n',i,j);
        
        [rmat nodes subjects] = matfile_loader(matfile);
        if ~isempty(xdistancemat)
            rmat = rmat .* xdistancemat;
        end
        
        % clear our a matrix
        %rmat=zeros(nodes,nodes,'single');
        
        % get the appropriate matrix
%         if bootstrapping
%             [rmat] = matrix_former(rmat,subjectarray(i,1),subjectarray(i,2),'2D','diagout',bootsamplesize);
%         else
%             [rmat] = matrix_former(rmat,subjectarray(i,1),subjectarray(i,2),'2D','diagout');
%         end
        
        % threshold the matrix
        [rmat foundr(i,j) foundkden(i,j)] = matrix_thresholder_Evan(rmat,thresholdarray(i,1),thresholdtype);
        
        % if user wants binarized networks
        if makebinary
            rmat=ceil(rmat);
        end
        
        % calculate the clusters for this network and store Q and partition
        switch clustertype
            case 'modularity'
                [rawclrs(:,i,j) Q(i,j)]=olaf_modularity_und(rmat);
            case 'infomap'
                pajekfile = [ outdir '/' filestem '_col' num2str(i) '.net' ];
                [rawclrs(:,i,j)] = infomap_wrapper_faster(roifile,rmat,pajekfile,100,1);
                %pause(.1)
%                 if issmallenough
%                     [Q(i,j)] = M_calc_modularity(rawclrs(:,i,j),rmat);
%                 end
        end
    end
end

% for simple analyses without any bootstrapping
    
    % write the raw assignments as .txt and .tiff
    dlmwrite([outdir '/rawassn.txt'],rawclrs,'\t');
    tiffmaker([outdir '/rawassn.tiff'],rawclrs);
    
%     % write the thr/box and Q values out
%     if issmallenough
%         fid=fopen([outdir '/thrbox_r_kden_Q.txt'],'w'); fprintf(fid,'Thr/box\tfoundr\tfoundkden\tQ\n'); fclose(fid);
%         dlmwrite([outdir '/thrbox_r_kden_Q.txt'],[xarray foundr foundkden Q],'delimiter','\t','-append');
%         plot(xarray,foundr,'r.',xarray,foundkden,'b.',xarray,Q,'g.'); xlabel('Thr/box'); ylabel('red:threshold blue:kden green:Q value');
%         saveas(gcf,[outdir '/thrbox_r_kden_Q.tiff'],'tiff');
%     else
%         fid=fopen([outdir '/thrbox_r_kden_Q.txt'],'w'); fprintf(fid,'Thr/box\tfoundr\tfoundkden\tQ\n'); fclose(fid);
%         dlmwrite([outdir '/thrbox_r_kden_Q.txt'],[xarray foundr foundkden],'delimiter','\t','-append');
%     end
%     
%     if issmallenough % don't attempt this with large networks
%         % now sort the raw assignments into sensible patterns
%         [pattern] =  rawoutput2clr(rawclrs);
%         dlmwrite([outdir '/pattern.txt'],pattern,'\t');
%         tiffmaker([outdir '/pattern.tiff'],pattern);
%     end
    
    % save a matfile of these variables
    if savematfile
        save([outdir '/' filestem '.mat']);
    end
    
%      if issmallenough
%          M_visuals(prmfile,analysistype,[outdir '/pattern.txt'],0,'visuals',[],xdistance,makebinary,1)
%      end
    



