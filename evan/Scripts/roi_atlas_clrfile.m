function roi_atlas_clrfile(roilist,clrfile,roirelationships,rgbfile,bycolumns,varargin)

% jdp 9/25/10
%
% This script accepts an roilist, a clrfile, and *an optional rgbfile*
% It first makes an roi atlas from the roilist, assigning voxels to nodes
% Then, using the clrfile, and some rgb palette (a default is set but the
% user can specify a custom palette), the colors are written to the voxel
% positions using .paint and .areacolor files for Caret.
%
% All of this is written to a /clrfile_paint4dfps/ directory
%
% USAGE: roi_atlas_clrfile(roilist,clrfile,roirelationships,rgbfile,bycolumns)

volvoxels=147456; % assume 333 space
switch roirelationships
    case 'between'
        [roireference roireferencefilename] = roi_atlas('between',roilist); % make the roi atlas
    case 'within'
        [roireference roireferencefilename] = roi_atlas('within',roilist);
    otherwise
        fprintf('Please choose ''between'' or ''within'' for the roirelationship\n');
end

% load modularity assignments
clrs=load(clrfile);
[c d]=size(clrs);
[pthstr namestem ext]=filenamefinder(clrfile,'dotsout');
[pthstr rgbstem ext]=filenamefinder(rgbfile,'dotsout');
writedir = [ pthstr '/' namestem '_' rgbstem '_paints/'];
if ~exist(writedir)
    mkdir(writedir);
end
allclrs=unique(clrs);
numclrs=size(allclrs,1);


% set a default palette that can be overriden with varargin
[clrrgbs shape border modulecolors moduleborders] = rgbmapper(clrs,bycolumns,rgbfile);

% go through all assignment columns
for k=1:d
    
    % load in the threshold assignments and start an output matrix
    clrcolumn=clrs(:,k);
    module4dfp=zeros(volvoxels,1);
    
    for i=1:c % cycle through each position of the clrcolumn
        clear a b;
        [a]=ismember(roireference,i); % find voxels for ROI #i
        module4dfp(a)=clrs(i,k); % set those voxels to the value of clrs(ROI#i,k)
    end
    % write the 4dfp.img and ifh file
    outputimg{k,1}=[ writedir '/' namestem '_col' num2str(k) '.4dfp.img' ];
    outputifhstem{k,1}=[ namestem '_col' num2str(k) '.4dfp.ifh' ];
    write_4dfpimg(module4dfp,outputimg{k,1},'bigendian');
    write_4dfpifh(outputimg{k,1},1,'bigendian');
end

% now we have a 4dfp with the assignment numbers at every voxel. we need to
% write a file that says the color of those assignment numbers
if ~isempty(varargin)
    % write areacolorfiles
    for k=1:d
        areacolorname{k,1}=[ writedir '/' namestem '_col' num2str(k) '.areacolor' ];
        areacolornamestem{k,1}=[ namestem '_col' num2str(k) '.areacolor' ];
        fid=fopen(areacolorname{k,1},'w');
        fprintf(fid,'CSVF-FILE,0,,,,,,,\n');
        fprintf(fid,'csvf-section-start,header,2,,,,,,\n');
        fprintf(fid,'tag,value,,,,,,,\n');
        fprintf(fid,'Caret-Version,5.614,,,,,,,\n');
        fprintf(fid,'Date,2010-06-15T13:57:00,,,,,,,\n');
        fprintf(fid,'comment,,,,,,,,\n');
        fprintf(fid,'encoding,COMMA_SEPARATED_VALUE_FILE,,,,,,,\n');
        fprintf(fid,'pubmed_id,,,,,,,,\n');
        fprintf(fid,'csvf-section-end,header,,,,,,,\n');
        fprintf(fid,'csvf-section-start,Colors,9,,,,,,\n');
        fprintf(fid,'Name,Red,Green,Blue,Alpha,Point-Size,Line-Size,Symbol,SuMSColorID\n');
        fprintf(fid,'???,170,170,170,255,2.000000,1.000000,POINT,\n');
        
        tempclrs=unique(clrs(:,k)); % find the colors that appear in a column
        temprgb=zeros(size(tempclrs,1),3); % these will be the RGB values used in this column
        for j=1:size(tempclrs,1) % write the rgb values for those colors to the areacolor file
            clear firstnode;
            firstnode=find((clrs(:,k)==tempclrs(j,1)),1); % find the first node with this clrfile value
            temprgb(j,:)=[clrrgbs(firstnode,k,1) clrrgbs(firstnode,k,2) clrrgbs(firstnode,k,3) ];
            fprintf(fid,'Unknown_name_%d,%d,%d,%d,255,2.000000,1.000000,POINT,\n',tempclrs(j,1),temprgb(j,1),temprgb(j,2),temprgb(j,3));
        end
        fprintf(fid,'csvf-section-end,Colors,,,,,,,,\n');
        fclose(fid);
        
        % write a 2d image of the file % gotta feed this the temprgb!!
        temprgb=zeroinsert(temprgb,tempclrs);
        fcimage_2dviewer(outputimg{k,1},temprgb,[ writedir '/' namestem '_col' num2str(k) '.tiff' ]);
        
    end
    
    areacolorlist=[ writedir '/' namestem '_caretpainterlist.txt' ];
    fid=fopen(areacolorlist,'w');
    for k=1:d
        fprintf(fid,'%s\t%s\n',areacolornamestem{k,1},outputifhstem{k,1});
    end
    fclose(fid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [temprgb] = zeroinsert(temprgb,clrs)

% this script is intended to map out 4dfps. real modules should start at 1,
% and the junk module is set to -1 currently. so there will be N modules,
% but the 4dfp also has values of zero everywhere where there is no module,
% and so the rgb mapping will be off by one since rgbmapper doesn't know to
% ignore those zero values. this finds where zero sits in the module
% hierarchy, and inserts a white color at that position.

zerorgb=[255 255 255];

if nnz(clrs==0)==0 % if zero is not in the present clr list
    insertposition=nnz(clrs<0);
    if insertposition==0 % if no current clrs are below zero it goes first
        temprgb=[zerorgb;temprgb];
    elseif insertposition==(size(temprgb,1)) % if all current clrs are below zero it goes last
        temprgb=[temprgb;zerorgb];
    else % if zerorgb goes in the middle somewhere
        temprgb=[temprgb(1:insertposition,:) ; zerorgb ; temprgb((insertposition+1):end,:)];
    end
end









