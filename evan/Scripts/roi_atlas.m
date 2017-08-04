function [roireference roireferencename] = roi_atlas(atlastype,roilist)

% jdp 9/25/10
% This script assumes 333 data. It also assumes that ROIs do not overlap
%
% You enter 2 things: atlastype and roilist. Atlastype instructs the script
% to make a particular form of the atlas, and roilist is a list of hard
% paths to ROI 4dfps.
%
% If 'between' is used as the atlastype, then all voxels in roi#x will
% receive the value x in the atlas. This can then be used to index which
% colors to assign voxels using caret_painter and the results of
% graphtools.
%
% If 'within' is used, the rois are read in, and all voxels in the rois are
% added to the atlas in the order of the rois. Each voxel will have some
% integer value, and these voxels will correspond to the results of
% graphtools output as well.
%
% Output is written to /roilist_atlastype_roiatlas.4dfp.img
%
% USAGE: [roireference reference4dfp] = roi_atlas(atlastype,roilist)
% USAGE: [roireference reference4dfp] = roi_atlas('between','roilist.txt')
% USAGE: [roireference reference4dfp] = roi_atlas('within','roilist.txt')

% assume 333
volumevoxels=147456;

%read in roilist
[rois x y z] = textread(roilist,'%s%f%f%f');
numrois=size(rois,1);

switch atlastype
    case 'between'
        % voxels in each ROI are assigned the position of the ROI in the
        % list. for example, voxels in ROI#3=3;
        roireference=zeros(volumevoxels,1);
        for j=1:numrois
            roidata=zeros(volumevoxels,1);
            [roidata frames voxels]=read_4dfpimg(rois{j,1});
            roidata=double(logical(roidata));
            roidata=roidata*j;
            roireference=roireference+roidata;
        end
        
        % this will tag the name of the roiatlas
        atlastypeinsert='between';
        
    case 'within'
        roireference=zeros(volumevoxels,1);
        collectedvoxels=0;
        % each voxel receives a unique number corresponding to its position
        % in the correlation matrix.
        for j=1:numrois
            roidata=zeros(volumevoxels,1); % clear and read the roidata
            [roidata frames voxels]=read_4dfpimg(rois{j,1});
            roidata=single(logical(roidata)); % set roi voxels to 1 or 0
            numvoxels=nnz(roidata); % how many voxels will we add
            clear srt ind; [srt ind]=sort(roidata); % sort adds an index position to each value (voxel) in the matrix
            maskvoxels=srt.*ind; % clear out the index values at voxels outside the mask
            maskvoxels(maskvoxels==0)=[]; % ditch those voxels from the array
            clear srt ind; [srt ind]=sort(maskvoxels); % this provides a new index to add to the already collected voxels
            ind=ind+collectedvoxels; % shift the index up by the number of voxels already collected
            roidata(maskvoxels)=ind; % now put those index numbers at the voxel positions of the roi
            collectedvoxels=collectedvoxels+numvoxels; % bump up the number of collected voxels
            roireference=roireference+roidata; % now add the newly indexed voxels to the master atlas
        end
      
        % this will tag the name of the roiatlas
        atlastypeinsert='within';
        
    otherwise
        fprintf('atlastype should be ''between'' or ''within''.\n');
end

[pth fname ext] = filenamefinder(roilist,'dotsout');
roireferencename = [ pth '/' fname '_' atlastypeinsert '_roiatlas.4dfp.img' ];
write_4dfpimg(roireference,roireferencename,'bigendian');
write_4dfpifh(roireferencename,1,'bigendian');

