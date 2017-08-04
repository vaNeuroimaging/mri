files = dir('./*.zip');
for filenum = 1:length(files)
     filenames = unzip(files(filenum).name,'Temp');
     for i = 1:length(filenames)
         if strcmp(filenames{i}(end-11:end),'.dscalar.nii')
             fileparts = tokenize(filenames{i},'/');
             copyfile(filenames{i},fileparts{end})
         end
     end
     rmdir('Temp','s')
end