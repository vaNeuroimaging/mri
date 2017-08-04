function caret_painter(areacolorlist)

% jdp 10/10/10
% 
% This script takes the output of roi_atlas_clrfile (a series of .imgs and
% corresponding .areacolor files), which is all stored in a list called the
% areacolorlist. This reads in the .img and .areacolor filenames, and then
% calls caret_4dfp_to_paint.csh, which takes the image volumes, projects
% them to surfaces, and then colors them according to the areacolor
% descriptions. The surfaces are snapshotted on fiducial, inflated, and
% very inflated surfaces, and .jpegs are saved accordingly.
% 
% USAGE: caret_painter(areacolorlist)


[areacolors images]=textread(areacolorlist,'%s%s');
[numfiles b]=size(areacolors);
[c d]=size(images);
if numfiles~=c
    error('Number of images and areacolor files are not the same!');
end


for i=1:numfiles
    thisimage=images{i,1};
    thisacfile=areacolors{i,1};
    commands=[ 'caret_4dfp_to_paint.csh ' thisimage ' ' thisacfile ];
    system(commands); 

    
    [pth,fname,ext]=filenamefinder(thisacfile,'dotsin');
    
    % load and resize images
    file1 = [ fname '_L.paint_LATERAL.jpg'];
    file2 = [ fname '_R.paint_LATERAL.jpg'];
    file3 = [ fname '_L.paint_MEDIAL.jpg'];
    file4 = [ fname '_R.paint_MEDIAL.jpg'];
    frame1=imread(file1); frame2=imread(file2); frame3=imread(file3); frame4=imread(file4);
    frame1b=frame1(160:870,400:1370,:); frame2b=frame2(160:870,220:1180,:); frame3b=frame3(160:870,220:1180,:); frame4b=frame4(160:870,400:1370,:);
    
    % combine images into one master image
    frame=[frame1b frame2b; frame3b frame4b];
    frameb=frame(1:2:end,1:2:end,:);
    jpegname=[ fname '.paint.jpeg' ];
    imwrite(frameb,jpegname,'jpeg');
    
    commands = ['rm ' file1 ' ' file2 ' ' file3 ' ' file4 ];
    system(commands);
    
    % load and resize images
    file1 = [ fname '_L.paint_LATERAL_INFLATED.jpg'];
    file2 = [ fname '_R.paint_LATERAL_INFLATED.jpg'];
    file3 = [ fname '_L.paint_MEDIAL_INFLATED.jpg'];
    file4 = [ fname '_R.paint_MEDIAL_INFLATED.jpg'];
    frame1=imread(file1); frame2=imread(file2); frame3=imread(file3); frame4=imread(file4);
    frame1b=frame1(200:1000,225:1335,:); frame2b=frame2(210:1010,255:1375,:); frame3b=frame3(220:990,265:1375,:); frame4b=frame4(230:1000,230:1350,:);
    
    % combine images into one master image
    frame=[frame1b frame2b; frame3b frame4b];
    frameb=frame(1:2:end,1:2:end,:);
    jpegname=[ fname '.paint_INFLATED.jpeg' ];
    imwrite(frame,jpegname,'jpeg');

    commands = ['rm ' file1 ' ' file2 ' ' file3 ' ' file4 ];
    system(commands);


    % load and resize images
    file1 = [ fname '_L.paint_LATERAL_VERY_INFLATED.jpg'];
    file2 = [ fname '_R.paint_LATERAL_VERY_INFLATED.jpg'];
    file3 = [ fname '_L.paint_MEDIAL_VERY_INFLATED.jpg'];
    file4 = [ fname '_R.paint_MEDIAL_VERY_INFLATED.jpg'];
    frame1=imread(file1); frame2=imread(file2); frame3=imread(file3); frame4=imread(file4);
    frame1b=frame1(200:1000,225:1335,:); frame2b=frame2(210:1010,255:1375,:); frame3b=frame3(220:990,265:1375,:); frame4b=frame4(230:1000,230:1350,:);
    
    % combine images into one master image
    frame=[frame1b frame2b; frame3b frame4b];
    frameb=frame(1:2:end,1:2:end,:);
    jpegname=[ fname '.paint_VERY_INFLATED.jpeg' ];
    imwrite(frame,jpegname,'jpeg');

    commands = ['rm ' file1 ' ' file2 ' ' file3 ' ' file4 ];
    system(commands);

     
end


    