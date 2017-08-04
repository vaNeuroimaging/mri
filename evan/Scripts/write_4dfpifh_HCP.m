function write_4dfpifh_HCP(imgname,frames,etype,x,y,z,dx,dy,dz)
%TOL modified 04/09/12
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
fprintf(fidifh,['matrix size [1] := ' num2str(x) '\n']);
fprintf(fidifh,['matrix size [2] := ' num2str(y) '\n']);
fprintf(fidifh,['matrix size [3] := ' num2str(z) '\n']);
fprintf(fidifh,'matrix size [4] := %d\n',frames);
fprintf(fidifh,['scaling factor (mm/pixel) [1]   := ' num2str(dx) '\n']);
fprintf(fidifh,['scaling factor (mm/pixel) [2]   := ' num2str(dy) '\n']);
fprintf(fidifh,['scaling factor (mm/pixel) [3]   := ' num2str(dz) '\n']);

fclose(fidifh);