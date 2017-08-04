function [datamat frames voxelsize] = read_4dfpimg(imgname)


% get necessary info from the ifhfile
[pth fname ext] = filenamefinder(imgname,'dotsin');
ifh = [ pth '/' fname '.ifh' ];
[voxelsize frames etype] = read_4dfpifh(ifh);
switch etype
    case 'littleendian'
        etypespec=['ieee-le'];
        [echeck] = endian_checker(imgname,'little');
    case 'bigendian'
        etypespec=['ieee-be'];
        [echeck] = endian_checker(imgname,'big');
    otherwise
        error('Endian type selected was neither big nor little..');
end

% read in the 4dfp
fid=fopen(imgname,'r',etypespec);
datamat = single(fread(fid,'float'));
fclose(fid);

% get the dimensions
d=size(datamat);

switch voxelsize
    case 3
        volumevoxels=147456;
        vols=d(1)/volumevoxels;
        if ~isequal(vols,frames)
            error('ifh frames and img volumes not matching');
        end
        datamat=reshape(datamat,[volumevoxels vols]);
    otherwise
        error('not set up for anything besides 333 yet');
end

