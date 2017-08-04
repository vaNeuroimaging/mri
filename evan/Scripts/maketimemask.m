subjects = {'vc34096' 'vc33457' 'vc34125' 'vc34126' 'vc34128' 'vc34140' 'vc34141' 'vc34198' 'vc34199' 'vc34200' 'vc34201' 'vc34220' 'vc33378' 'vc35175' 'vc34252' 'vc33775_2' 'vc35469' 'vc34306' 'vc34307' 'vc34308' 'vc34330' 'vc34331' 'vc34401' 'vc34402' 'vc34403' 'vc34404' 'vc33769' 'vc34408'};

directory = [];
outputfilename = 'Tmask.txt';

delete([directory outputfilename]);
fid = fopen([directory outputfilename],'at'); %open the output file for writing
fprintf(fid,'%s\t\%s\t\%s\n\r\','Subject','Timepoint','Value'); %write the output file header
fclose(fid);
dlmwrite([directory outputfilename],' ','-append');

for subject = 1:length(subjects)
     mask = dlmread(['/data/cn4/evan/RestingState/28subjects_FDp2_00_DVinit_20_00/tmask/28subjects_FDp2_00_DVinit_20_00_' subjects{subject} '_tmask.txt']);
     
     for timepoint = 1:length(mask)
    dlmwrite([directory outputfilename],[subjects{subject} '  ' num2str(timepoint) '   ' num2str(mask(timepoint))],'-append','delimiter','');%write the data to the output file
     end

end

