function soniawriter_Evan(mat,roinames,rgb,border,nodesize,borderwidth,edgergb,outputname,varargin)

% jdp 10/10/10
%
% This script takes networks and turns them into SoNIA files to see spring
% embedded layouts. It needs
% 
% mat: a 3D or 2D matrix of the network
% roinames: [nodes x 1] names of the nodes
% rgb: [node x analyses x 3] 0-255 rgb matrix
% border: [node x analyses] a cell array of node border colors
% nodesize: [node x analyses ] node sizes (15 is standard)
% borderwidth: [node x analyses] border widths (2 is standard)
% arccolors: a switch for coloring purposes
%   black: all edges black
%   red/blue: positive red, negative blue
%   custom: probably jet, but some color scale has been used
% outputname: the filename to write the data to 

if ~isempty(varargin)
    edgergb=varargin{1,1};
end

% get matrix dimensions
d=size(mat);

% presuming variables are correctly formed
numanalyses=size(borderwidth,2);

% presuming the user passes in a 0-255 RGB matrix
rgb=rgb/255;

% open the sonia file, write node information
fid=fopen(outputname,'w');
fprintf(fid,'AlphaId\tLabel\tRedRGB\tGreenRGB\tBlueRGB\tStartTime\tEndTime\tBorderColor\tBorderWidth\tNodeSize\n');
for i=1:d(1)
    for j=1:numanalyses
        starttime=j-1; endtime=j;
        fprintf(fid,'%d\t%s\t%g\t%g\t%g\t%g\t%g\t%g\t%g\t%g\\t%g\t%g\n',i,roinames{i},rgb(i,j,1),rgb(i,j,2),rgb(i,j,3),starttime,endtime,border(i,j,1),border(i,j,3),border(i,j,3),borderwidth(i,j),nodesize(i,j));
    end
end

% now write the edge information
fprintf(fid,'FromId\tToId\tStartTime\tEndTime\tArcWeight\tRedRGB\tGreenRGB\tBlueRGB\n');
for k=1:numanalyses    
    for i=1:d(1)
        for j=i+1:d(1)
            if abs(mat(i,j,k))>0.001
                if ~isnumeric(edgergb)
                    switch edgergb
                        case 'black'
                            arccolor=[0 0 0];
                        case 'red/blue'
                            if mat(i,j,k)>0;
                                arccolor=[1 0 0];
                            else
                                arccolor=[0 0 1];
                            end
                    end
                else
                    arccolor(1,1)=edgergb(i,j,k,1);
                    arccolor(1,2)=edgergb(i,j,k,2);
                    arccolor(1,3)=edgergb(i,j,k,3);
                    arccolor=arccolor/255;
                end
                fprintf(fid,'%d\t%d\t%g\t%g\t%g\t%f\t%f\t%f\n',i,j,k-1,k,mat(i,j,k),arccolor(1,1),arccolor(1,2),arccolor(1,3));
            end
        end
    end
end

% close the file
fclose(fid);