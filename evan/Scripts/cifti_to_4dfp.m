function cifti_to_4dfp(ciftifile,manualmask)
%cifti_to_4dfp(ciftifile,[mask]);

slashloc = strfind(ciftifile,'/');
if isempty(slashloc)
    slashloc = 0;
end

ciftifile_noext = ciftifile(slashloc(end)+1:end-13);
result = evalc(['!wb_command -cifti-convert -to-gifti-ext ' ciftifile ' ' ciftifile_noext '.func.gii']);
if (~strcmp(result,'')) && all(result ~=0)
    error(result)
end
data = gifti([ciftifile_noext '.func.gii']);
data = data.cdata;
delete([ciftifile_noext '.func.gii'])
delete([ciftifile_noext '.func.gii.data'])

%[maskdata frames voxelsize] = read_4dfpimg_HCP(maskfile);
%[voxelsize frames x y z etype] = read_4dfpifh_HCP([maskfile(1:end-3) 'ifh']);

enoughvoxels = 0;
resolution = 3;
while enoughvoxels == 0
    
    if exist('manualmask')
        maskfile = manualmask;
    else
        maskfile = ['/home/usr/fidl/lib/glm_atlas_mask_' num2str(resolution) num2str(resolution) num2str(resolution) '.4dfp.img'];
    end
    
    [maskdata, frames, voxelsize, x, y, z] = read_4dfpimg_HCP_noendian(maskfile);
    maskinds = find(maskdata);
    
    if size(data,1) > length(maskinds)
        if (resolution==1) || (exist('manualmask'))
            error('Insufficient data points inside mask to contain all the cifti data!')
        else
            resolution = resolution - 1;
        end
    else
        enoughvoxels = 1;
    end
end

voldata = rand(size(maskdata,1),size(data,2));
voldata(maskinds(1:size(data,1)),:) = data;

write_4dfpimg(voldata,[ciftifile_noext '_' num2str(resolution) num2str(resolution) num2str(resolution) '.4dfp.img'],'bigendian')
write_4dfpifh_HCP([ciftifile_noext '_' num2str(resolution) num2str(resolution) num2str(resolution) '.4dfp.img'],size(voldata,2),'bigendian',x,y,z,voxelsize,voxelsize,voxelsize)


fid = fopen([ciftifile_noext '.txt'],'at'); %open the output file for writing
fprintf(fid,'%s\n',['Mask ' maskfile]); %
fprintf(fid,'%s',['Cifti ' ciftifile]); %write the output file header
fclose(fid);

