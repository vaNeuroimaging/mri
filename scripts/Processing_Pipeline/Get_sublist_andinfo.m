projectstring = 'MAV';
scanstring = 'RSFC';
subfolder = '/home/data/subjects/';

outputfilename = [subfolder 'MAV_subs_info.txt'];
delete([outputfilename]);
fid = fopen([outputfilename],'at'); %open the output file for writing
fprintf(fid,'%s\t\%s\t\%s\t\%s\t\%s\n\r\','Subject','Runs','Sessions','meanFD','PctRetained'); %write the output file header
fclose(fid);
dlmwrite([outputfilename],' ','-append');

cd(subfolder);
allsubjects = dir([projectstring '*']);

for subnum = 1:length(allsubjects)
    runsessfile = [subfolder allsubjects(subnum).name '/fc_processed/' scanstring '_runs_sessions.txt'];
    if exist(runsessfile)
        [runs, sessions] = textread(runsessfile,'%f%f');
        tmask = load([subfolder allsubjects(subnum).name '/fc_processed/' scanstring '_all_tmask.txt']);
        FD = load([subfolder allsubjects(subnum).name '/fc_processed/' scanstring '_all_FD.txt']);
        
        outstring = [allsubjects(subnum).name ' ' num2str(nnz(unique(runs))) ' ' num2str(nnz(unique(sessions))) ' ' num2str(mean(FD)) ' ' num2str(nnz(tmask) / numel(tmask))];
        dlmwrite([outputfilename],outstring,'-append','delimiter','');
    end
end

