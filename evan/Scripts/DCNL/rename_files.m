% correctname = '256';
% wrongsubname = '296';

files_directory = '/fmri/data3/Evan/Gene-Rest-Nback/Data/307/avw_wrongname/';
files_destination = '/fmri/data3/Evan/Gene-Rest-Nback/Data/307/avw/';

files = dir(files_directory);
files = files(3:end);

for filenum = 1:length(files)
    oldname = files(filenum).name;
%     namelocation = strfind(oldname,wrongsubname);
    %newname = [oldname(1:namelocation-1) correctname oldname(namelocation+3:end)];
    newname = ['DCNL-' oldname];
    copyfile([files_directory oldname], [files_destination newname]);
end