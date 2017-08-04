function M_visuals_evan(prmfile,analysistype,clrfile,bycolumns,rgbvals,networknames,ROIreflectdistance,varargin)
%
% Name:M_visuals.m
% $Revision: 1.3 $
% $Date: 2011/03/08 20:29:28 $
%
% jdp 10/10/10
% 
% This script takes a network and makes a lot of visuals for it
% 
% There's no simple way to describe this script. It needs the full array of
% input used for most graphtools analyses, namely the prmfile and the
% analysistype. It then needs a set of assignments (or some parameter to
% map onto nodes), in the clrfile matrix. The values in the clrfile will be
% scale colors using either 1) the whole matrix, or 2) just the values in
% each column, controlled by the bycolumns switch. The colors that will map
% onto those values are controlled by rgbfile, which can be any MATLAB
% palette (e.g. 'jet','hsv',etc), or it can be a user-defined file of 3
% columns of 0-255 RGB values (e.g., 'myrgbfile.txt'), or even a single rgb
% color like [0 128 255]. ROIreflectdistance is used in the caret files,
% and will reflect rois within X mm of the midline to the other side, for
% when you just want to show a single hemisphere but want to see the
% midline nodes from either side. 
% 
% By default, no sonia files are written. SoNIA kinda sucks with larger
% files, so I don't recommend it for voxelwise analyses with lots of nodes,
% or if you're going to do a lot of thresholds. If you want SoNIA files
% written, then you supply additional variables, beginning with xdistance
% and makebinary (which eliminate edges of nodes within xdistance mm of one
% another, and either make the graph unweighted or leave it weighted).
% Additionall, one can pass in matrices or single entries for border,
% borderwidth, and nodesize (e.g. 'black',2,15, or matrices of these things
% that are [node x numanalyses] in if you want everything controlled in
% exact detail.
% 
% Finally, you can also pass in an edgergb and matbymat argument. Edgergb
% specifies how to color edges. This can be 'black' or 'red/blue', in which
% case edges are black, or red=positive/blue=negative. You can also supply
% any MATLAB colormap ('jet','hsv'), and matbymat will apply this colormap
% scaled over the entire matrix (matbymat=0), or scaled to each 3rd
% dimension matrix (matbymat=1). Alternatively, if you have a [node x node
% x numanalyses x 3] RGB matrix of 0-255 values, you can supply this
% instead.
% 
% prmfile: typical graphtools
% analysistype: typical graphtools 'thr' or 'box'
% clrfile: a matrix of integers or values of some sort
% bycolumns: whether to scale coloring over the whole matrix (0) or within
%   each column (1)
% rgbfile: the colormap or 3-column 0-255 RGB text file
% ROIreflectdistance: reflect ROIs within X mm of midline
%   use [] to avoid any reflecting
% 
% *optional variables, for sonia output*
% xdistance: typical graphtools, exclude edges of nodes <xdistance apart
% makebinary: 0-leave matrix weighted, 1-binarize matrix
% border: a single entry, or cell matrix [node x numanalyses] of strings
%   for SoNIA node borders.
%   standard sonia colors: 'blue','red',etc are also possibilities
% borderwidth: same as border. DEFAULT=2
% nodesize: scaling factor, same as border. DEFAULT=15
% edgergb: can be 'black','red/blue', or a MATLAB colormap, or a [node node
%   numanalyses 3] 0-255 RGB array. DEFAULT='black'
% matbymat: whether to scale edge coloring by the entire matrix (0), or by 
%   the matrix at each position in the 3rd dimension (1). DEFAULT=0
% metanodes: combined like nodes into supernodes, and write sonia files
%   for these supernodes, as well as the supernodes with edges normalized by
%   the number of nodes within the supernodes (possible edges between nodes
%   = n1 x n2)
% 
% USAGE: M_visuals(prmfile,analysistype,clrfile,colbycol,rgbfile,ROIreflectdistance,*xdistance,makebinary,metanodes*,**border,borderwidth,nodesize**,***edgergb,matbymat***)
% USAGE: M_visuals('modbox.prm','thr','clr.txt',0,'jet',[])                     nodes colored by 'jet', no reflection, all clrs scaled at once
% USAGE: M_visuals('modbox.prm','thr','clr.txt',1,'hsv',5,[],0,1)               nodes colored by 'hsv', scaled column by column, reflected within 5mm of midline, make sonia using prmfile settins with xdist=[], makebinary=0, make metanodes
% USAGE: M_visuals('modbox.prm','thr','clr.txt',0,'cool',5,15,1,'black',2,15,0)   nodes colored by 'cool', reflect 5mm midline, 15mm xdistance and make binary, make edges black with borderwidth 2 and nodesize 15
% USAGE: M_visuals('modbox.prm','thr','clr.txt',0,'gray',[],15,1,0,'black',2,15,'red/blue',1)
% USAGE: M_visuals('modbox.prm','thr','clr.txt',0,'jet',[],15,1,0,'black',2,15,'hot',0)      
% USAGE: M_visuals('modbox.prm','thr','clr.txt',0,'jet',[],15,1,0,'black',2,15,'jet',0)        color edges with 'jet', whole matrix at once, and write out normalized versions of everything
% 
% NOTES: 11/10/10 fixed 'border' with char instead of cell array






%%%%% FIRST, A LOT OF DATA NEEDS READING IN AND MASSAGING %%%%%


% read in the prmfile settings
[matfile roifile stemname subjectA subjectZ loend step hiend writepath threshold boxcarsize boxcarstep thresholdtype] = prmfilereader(prmfile);

% load in xyz coordinates from roifile
[xyz name]=roifilereader(roifile);

% create filebase for output
[filestem] = filenameprep(subjectA,subjectZ,loend,step,hiend,threshold,boxcarsize,boxcarstep,analysistype,thresholdtype);

% determine the # thresholds, etc. for the analysis
[subjectarray thresholdarray numanalyses xarray] = matrix_parameter_setter(subjectA,subjectZ,loend,step,hiend,threshold,boxcarsize,boxcarstep,analysistype);

% make a preliminary stemname for all output
[pth, clrname, ext] = filenamefinder(clrfile,'dotsout');
if ~isnumeric(rgbvals)
    [trash, rgbname, ext] = filenamefinder(rgbvals,'dotsout');
else
    rgbname=['customrgb'];
end
filestem = [ clrname '_' rgbname ];

% load in color matrix
clrmatrix=load(clrfile);
d=size(clrmatrix);
if d(1)~=size(xyz,1)
    error('Clrfile doesn''t have the same number of nodes as the roifile.');
end

% obtain the rgb and caret/sonia parameters for this color matrix
[rgb shape border modulecolors moduleborders] = rgbmapper(clrmatrix,bycolumns,'hot');

% reflect nodes around midline if desired
reflected=0;
startnumnodes=size(xyz,1);
if ~isempty(ROIreflectdistance)
    reflectedones=abs(xyz(:,1))<=ROIreflectdistance;
    if ~isempty(reflectedones)
        
        % indicate that we reflected
        reflected=1;
        
        % initialize reflected matrices
        refxyz=xyz;
        refshape=shape;
        refrgb=rgb;
        refname=name;
        
        % fill in reflected matrices
        refxyz(startnumnodes+1:startnumnodes+nnz(reflectedones),1)=-xyz(reflectedones,1);
        refxyz(startnumnodes+1:startnumnodes+nnz(reflectedones),2:3)=xyz(reflectedones,2:3);
        refshape(startnumnodes+1:startnumnodes+nnz(reflectedones),:)=shape(reflectedones,:);
        refrgb(startnumnodes+1:startnumnodes+nnz(reflectedones),:,:)=rgb(reflectedones,:,:);
        for i=startnumnodes+1:startnumnodes+nnz(reflectedones)
            refname{i,1}=[num2str(refxyz(i,1)) '_' num2str(refxyz(i,2)) '_' num2str(refxyz(i,3)) ];
        end
        
    end
    filestem = [ filestem '_RF' num2str(ROIreflectdistance) ];
    fprintf('Reflecting nodes within %4.2f mm of midline increased the number of nodes from %d to %d.\n',ROIreflectdistance,startnumnodes,size(refxyz,1));
end

% make a path for output data
outdir = [ pth '/' filestem ];
if ~exist(outdir)
    mkdir(outdir);
end



%%%%% NOW START WRITING SOME OUTPUT, STARTING WITH CARET FILES %%%%%



% save the RGB colors that are used
dlmwrite([outdir '/' filestem '_rgbmat.txt'],modulecolors,'\t');

% make a .tiff of the clr assignments
imwrite(imresize(double(rgb/255),10,'nearest'),[outdir '/' filestem '.tiff'],'tiff');

% write the caret files
if reflected
    for i=1:numanalyses
        M_caretfilemaker(refname,refxyz,refrgb(:,i,:),refshape(:,i),[outdir '/' filestem '.foci'],[outdir '/' filestem '_col' num2str(i) '.focicolor']);
    end
else
    for i=1:numanalyses
        M_caretfilemaker(name,xyz,rgb(:,i,:),shape(:,i),[outdir '/' filestem '.foci'],[outdir '/' filestem '_col' num2str(i) '.focicolor']);
    end
end

% write list of caret files for making movies
fid=fopen([ outdir '/' filestem '_focilist.txt' ],'w');
fprintf(fid,'%s\n',[filestem '.foci']);
for i=1:numanalyses
    fprintf(fid,'%s\n',[filestem '_col' num2str(i) '.focicolor']);
end
fclose(fid);


%%%%% NOW MAKE SONIA FILES IF DESIRED %%%%%

% the file sizes can get quite large here if you don't watch out

if ~isempty(varargin)
    
    %%% SET SONIA PARAMETERS %%%
    %%% soniawriter is completely customizable, but varying levels of
    %%% detail will be needed depending on the analysis. set defaults but
    %%% let users override them
    
    nvargin=size(varargin,2);
    
    % set defaults (to be overridden as needed)
    xdistance=[];
    makebinary=0;
    donormalized=1;
    borderwidth=2;
    nodesize=15;
    edgergb='black';
    matbymat=0;

    
    % needed to properly form the matrix
    xdistance=varargin{1,1};
    makebinary=varargin{1,2};
    donormalized=varargin{1,3};
      
    % to make nodes
    if ((nvargin==6) || (nvargin==8));
        border=varargin{1,4};
        borderwidth=varargin{1,5};
        nodesize=varargin{1,6};
    end
    
    % to color edges in a particular way
    if (nvargin==8);
        edgergb=varargin{1,7};
        matbymat=varargin{1,8};
    end

    
    %%% MAKE A 3D MATRIX CORRESPONDING TO THE ANALYSIS %%%
    
    % load the matrix
    [matrix nodes subjects] = matfile_loader(matfile);
    
    % apply a distance exclusion?
    if ~isempty(xdistance)
        if isnumeric(xdistance) && (xdistance>=0)
            [matrix] = matrix_xdistance(matrix,roifile,xdistance);
            xd=dotremover(num2str(xdistance));
            filestem=[filestem '_xd' xd ];
        else
            error('xdistance is not >=0 or is not numeric.\n');
        end
    end
    
    % check makebinary
    switch makebinary
        case 0
        case 1
            filestem=[filestem '_BI' ];
        otherwise
            error('makebinary should be ''0'' or ''1'' ');
    end
    
    % make matrices
    rmat=zeros(nodes,nodes,numanalyses,'single');
    for i=1:numanalyses
        % pull the matrix out of the master matrix
        [rmat(:,:,i)] = matrix_former(matrix,subjectarray(i,1),subjectarray(i,2),'2D','diagout');
        
        % threshold the matrix
        [rmat(:,:,i) foundr(i,1) foundkden(i,1)] = matrix_thresholder(rmat(:,:,i),thresholdarray(i,1),thresholdtype);      
    end
    clear matrix;
    
    % if user wants binarized networks
    if makebinary
        rmat=ceil(rmat);
    end
    
    %%% NOW THE (bigass) MATRIX IS PROPERLY FORMED %%%
    %%% NOW ADJUST SONIA PARAMETERS TO MATCH THE MATRIX %%%
    
    % adjust border to [nodes numanalyses]
    if ischar(border)
        newborder={border};
        border=repmat(newborder,[nodes numanalyses]);
    elseif isequal(size(border),[nodes numanalyses])
    else
        error('Border size should be either a single cell entry, or a cell matrix of [nodes x numanalyses].');
    end
    
    % adjust borderwidth to [nodes numanalyses]
    if ~isequal(size(borderwidth),[nodes numanalyses]) && isequal(size(borderwidth),[1 1])
        oldborderwidth=single(borderwidth);
        borderwidth = repmat(oldborderwidth,[nodes numanalyses]);
    else
        error('Borderwidth should be either a single entry, or a matrix of [nodes x numanalyses].');
    end
    
    % adjust nodesize to [nodes numanalyses]
    if ~isequal(size(nodesize),[nodes numanalyses]) && isequal(size(nodesize),[1 1])
        oldnodesize=single(nodesize);
        nodesize = repmat(oldnodesize,[nodes numanalyses]);
    else
        error('Nodesize should be either a single entry, or a matrix of [nodes x numanalyses].');
    end
    
    % adjust edgergb to [nodes nodes numanalyses 3]
    if ~isnumeric(edgergb)
        switch edgergb
            case {'black','red/blue'}
            case {'hsv','hot','gray','bone','copper','pink','white','flag','lines','colorcube','vga','jet','prism','cool','autumn','spring','winter','summer'}
                [edgergb]=matrixrgbmapper(rmat,edgergb,matbymat);
            otherwise
                error('If not passing a rgb matrix in, need to use ''black'',''red/blue'', or a MATLAB colormap, e.g. ''jet'' ');
        end
    else
        if ~isequal((size(edgergb)),[nodes nodes numanalyses 3])
            error('edgergb needs to be a matrix of [nodes nodes numanalyses 3] of 0-255 values.');
        end
    end
    
    %%% NOW EVERYTHING SHOULD MATCH UP AND WE CAN WRITE THE SONIA FILE %%%

    soniawriter(rmat,networknames,rgbvals,border,nodesize,borderwidth,edgergb,[ outdir '/' filestem '.son' ])
    
    
    %%% WE CAN ALSO CONSOLIDATE NODES INTO MODULENODES %%%
    
    if donormalized    
        
        % i'm really presuming these are modularity assignments. it's meaningless otherwise.
        [MDCX nMDCX modulenodesize modulename modulenumber firstexamplex firstexampley] = node2module_matrix(clrmatrix,rmat);
        [n1 n2]=size(modulenodesize);
        nmodulenodesize=repmat(single(15),[n1 n2]);
        moduleborderwidth=repmat(single(2),[n1 n2]);
        moduleedgergb='black';
        fprintf('Setting colormap to ''black'' for edges of the modulenodes.son files\n');
        fprintf('Setting borderwidth to 2 for modulenodes.son and nmodulenodes.son files\n');
        fprintf('Setting nodesize to 15 for nmodulenodes.son files\n');
        
        % set border, borderwidth, rgb stuff
        moduleborder=cell(n1,n2);
        modulergb=zeros(n1,n2,3,'single');
        for i=1:n1
            for j=1:n2
                moduleborder(i,j)=border(firstexamplex(i),firstexampley(i));
                modulergb(i,j,1)=rgb(firstexamplex(i),firstexampley(i),1);
                modulergb(i,j,2)=rgb(firstexamplex(i),firstexampley(i),2);
                modulergb(i,j,3)=rgb(firstexamplex(i),firstexampley(i),3);
            end
        end

        soniawriter(MDCX,modulename,modulergb,moduleborder,modulenodesize,moduleborderwidth,moduleedgergb,[ outdir '/' filestem '_modulenodes.son' ]);
        soniawriter(nMDCX,modulename,modulergb,moduleborder,nmodulenodesize,moduleborderwidth,moduleedgergb,[ outdir '/' filestem '_modulenodes_norm.son' ]);
    end
    
end







