function Batch_workbench_pics(scene,filenames,timepoints)

bufsize = 524288;
scenetext = textread(scene,'%s','delimiter','\r','bufsize',bufsize);

for timepoint = timepoints(:)'
    
    delete('Temp_forpics.scene')
    fid = fopen('Temp_forpics.scene','at'); %open the output file for writing
    fclose(fid);

    overwrite = 0;
    dontwrite = 0;
for row = 1:length(scenetext)
    
    for filenum = 1:length(filenames)
        if strcmp(scenetext{row},['File: ' filenames{filenum}])
            overwrite = row+1;
        end
    end
    if strcmp('<Image Encoding="Base64"',scenetext{row})
        dontwrite = [row:(row+2)];
    end
    
    
    if row==overwrite
        dlmwrite('Temp_forpics.scene',['Map Index: ' num2str(timepoint)],'-append','delimiter','');
    elseif ~any(row==dontwrite)
        dlmwrite('Temp_forpics.scene',scenetext{row},'-append','delimiter','');
    end
    
    
end

system(['wb_command -show-scene Temp_forpics.scene 1 outpic_' num2str(timepoint) '.png 1900 500'])
end
    
            