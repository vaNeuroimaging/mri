function [voxelsize frames etype] = read_4dfpifh(ifh)

commands=[ 'grep "matrix.*\[4\]" ' ifh ' | awk -F ":= " ''{print $2}''' ];
[trash,frames] = system(commands);
commands=[ 'grep "scaling.*\[1\]" ' ifh ' | awk -F ":= " ''{print $2}''' ];
[trash,voxelsize] = system(commands);
commands=[ 'grep "imagedata" ' ifh ' | awk -F ":= " ''{print $2}''' ];
[trash,etype] = system(commands);
etype=etype(1:end-1); % removes the carriage return

voxelsize=str2num(voxelsize);
frames=str2num(frames);