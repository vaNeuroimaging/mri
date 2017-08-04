function [datamat frames voxelsize I J K] = read_4dfpimg_HCP_noendian(imgname)


result = evalc(['!endian_4dfp ' imgname]);
databegins = strfind(result,imgname) + length(imgname);
for i = databegins + 10 : length(result)
    if isempty(str2num(result(databegins : i)))
        dataends = i-1;
        break
    end
end

data = str2num(result(databegins : dataends));
I = data(1); J = data(2); K = data(3); frames = data(4);

databegins = dataends + 1;
for i = databegins+10 : length(result)
    if isempty(str2num(result(databegins : i)))
        dataends = i-1;
        break
    end
end
data = str2num(result(databegins : dataends));
voxelsize = data(1);


etype = endian_checker(imgname);
if strcmp(etype,'big')
    etypespec=['ieee-be'];
elseif strcmp(etype,'little')
    etypespec=['ieee-le'];
else
    error('Endian type selected was neither big nor little..');
end


% % get necessary info from the ifhfile
% [pth fname ext] = filenamefinder(imgname,'dotsin');
% ifh = [ pth '/' fname '.ifh' ];
% [voxelsize frames I J K etype] = read_4dfpifh_HCP(ifh);
% switch etype
%     case 'littleendian'
%         etypespec=['ieee-le'];
%         [echeck] = endian_checker(imgname,'little');
%     case 'bigendian'
%         etypespec=['ieee-be'];
%         [echeck] = endian_checker(imgname,'big');
%     otherwise
%         error('Endian type selected was neither big nor little..');
% end

% read in the 4dfp
fid=fopen(imgname,'r',etypespec);
datamat = single(fread(fid,'float'));
fclose(fid);

% get the dimensions
d=size(datamat);

% for i=1:length(I)
%     if length(str2num(I(1:i))) == 1
%         numI = str2num(I(1:i));
%     end
% end
% 
% for i=1:length(J)
%     if length(str2num(J(1:i))) == 1
%         numJ = str2num(J(1:i));
%     end
% end
% 
% for i=1:length(K)
%     if length(str2num(K(1:i))) == 1
%         numK = str2num(K(1:i));
%     end
% end
% 
% volumevoxels = numI*numJ*numK;
volumevoxels = I*J*K;
vols = d(1)/volumevoxels;
datamat = reshape(datamat,[volumevoxels vols]);