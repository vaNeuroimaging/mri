function write_4dfpifh(imgname,frames,etype,varargin)

% create ifhname
[pth fname ext] = filenamefinder(imgname,'dotsin');
ifhname=[ pth '/' fname '.ifh' ];

% write ifh
fidifh=fopen(ifhname,'wt');
fprintf(fidifh,'INTERFILE       :=\n');
fprintf(fidifh,'version of keys := 3.3\n');
fprintf(fidifh,'number format           := float\n');
fprintf(fidifh,'conversion program      := write_ifh.m\n');
fprintf(fidifh,'name of data file       := %s\n',fname);
fprintf(fidifh,'number of bytes per pixel       := 4\n');
fprintf(fidifh,'imagedata byte order    := %s\n',etype);
fprintf(fidifh,'orientation             := 2\n');
fprintf(fidifh,'number of dimensions    := 4\n');
fprintf(fidifh,'matrix size [1] := 48\n');
fprintf(fidifh,'matrix size [2] := 64\n');
fprintf(fidifh,'matrix size [3] := 48\n');
fprintf(fidifh,'matrix size [4] := %d\n',frames);
fprintf(fidifh,'scaling factor (mm/pixel) [1]   := 3.000000\n');
fprintf(fidifh,'scaling factor (mm/pixel) [2]   := 3.000000\n');
fprintf(fidifh,'scaling factor (mm/pixel) [3]   := 3.000000\n');
fprintf(fidifh,'mmppix  :=   3.000000 -3.000000 -3.000000\n');
fprintf(fidifh,'center  :=    73.5000  -87.0000  -84.0000\n');

if ~isempty(varargin) %%% THIS IS NOT FINISHED HERE %%%
    numvoxels=varargin{1,1};
    fprintf(fidifh,'region names := 0 region %d',numvoxels);
end

fclose(fidifh);