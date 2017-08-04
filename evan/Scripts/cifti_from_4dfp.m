function cifti_from_4dfp(volfile,textfile)
%cifti_from_4dfp(volfile,textfile);


[fields files] = textread(textfile,'%s%s');

for i = 1:length(fields)
    if strcmp(fields{i},'Mask')
        maskdata = read_4dfpimg_HCP_noendian(files{i});
    elseif strcmp(fields{i},'Cifti')
        ciftifile = files{i};
    end
end

[data frames voxelsize] = read_4dfpimg_HCP(volfile);

data = data(logical(maskdata),:);

evalc(['!wb_command -cifti-convert -to-gifti-ext ' ciftifile ' ' ciftifile(1:end-13) '.func.gii']);
temp = gifti([ciftifile(1:end-13) '.func.gii']);
cifti_write_wHDR(data(1:size(temp.cdata,1),:),[ciftifile(1:end-13) '.func.gii'],[volfile(1:end-9) '.dtseries.nii'])
delete([ciftifile(1:end-13) '.func.gii*'])
