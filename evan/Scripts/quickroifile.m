function quickroifile(numrois,varargin)

% jdp 10/10/10
% This is a versatile script that makes roifiles. It works in several ways, illustrated below:
% 
% The roifile is modeled after the standard SoNIA roifile (the header is ignored):
% 
% ID X Y Z NAME CATA COLORA CATB COLORB
% 1 3 7 6 PCC Brain Black Brain Black
% 2 6 2 2 ACC Brain Black Brain Black 
% 3 9 1 3 pFC Brain Black Brain Black
% ...
%
% USAGE: quickroifile(numrois/roifilename/roiarray,filename)
% USAGE: quickroifile(23) - makes 23.roi, an roifile with 23 dummy positions
% USAGE: quickroifile(23,'mine') - makes mine.roi, an roifile with 23 dummy positions
% USAGE: quickroifile('xyzlist.txt') - makes xyzlist.roi, with rois corresponding to the coordinates
% USAGE: quickroifile('xyzlist.txt','newroifile') - makes newroifile.roi, with rois corresponding to the coordinates
% USAGE: quickroifile(xyz) - makes 44.roi, for an xyz array of 44 coordinates
% USAGE: quickroifile(xyz,'mynewrois') - makes mynewrois.roi, with rois cooresponding to the xyz coordinates
% 
% NOTES:
% 10/13/10 updated the header for graphtools

if isnumeric(numrois)
    if size(numrois,2)==3
        xyz=numrois;
    else
        xyz=1:numrois;
        xyz=xyz';
        xyz=repmat(xyz,[1 3]);
    end
    filename = [ pwd '/' num2str(size(xyz,1)) '.roi' ];
elseif exist(numrois)==2 % if you have an xyz list you want make into an roifile
    xyz=load(numrois);
    [pth fname ext]=filenamefinder(numrois,'dotsout');
    filename=[ pth '/' fname '.roi' ];
end


if ~isempty(varargin)
    filename=varargin{1,1};
end
    

fid=fopen(filename,'w');
fprintf(fid,'ROI\tX\tY\tZ\tName\tAnatomy\tAcolor\tBname\tBcolor\n');
for i=1:size(xyz,1);
    clear tempname;
    tempname=[ num2str(i) ];
    fprintf(fid,'%d\t%d\t%d\t%d\t%s\tBrain\tBlack\tBrain\tBlack\n',i,xyz(i,1),xyz(i,2),xyz(i,3),tempname);
end
fclose(fid);

    
