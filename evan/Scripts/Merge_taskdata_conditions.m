datafolder = '/data/cn3/steven/Task_fcMapping/INDIVIDUALS/';

files = dir([datafolder '*.img']);
files2 = struct2cell(files);
condnames = [];
for i=1:length(files2)
    if strcmp(files2{1,i}(1:3),'_vc')
    condnames{end+1} = files2{1,i}(end-70:end);
    end
end
conditions = unique(condnames);

finaldata = [];

for condition = 1:length(conditions)
    
%     underscoreindices = strfind(conditions{condition},'_');
%     
%     if length(underscoreindices)>1
%         shortsubname = conditions{condition}(2:underscoreindices(2)-1);
%     else
%         shortsubname = conditions{condition}(2:end);
%     end
    
    %allshortsubnames{subject,1} = shortsubname;
    
    files = dir([datafolder '*' conditions{condition}]);
    
    for file = 1:length(files)
        data(:,:,file) = read_4dfpimg([datafolder files(file).name]);
    end
    
    finaldata = [finaldata mean(data,3)];
    
    
    
end

write_4dfpimg(finaldata,'/data/cn4/evan/Task_parcellation/SteveSubjects/AvgNewtaskTimecourses.4dfp.img','bigendian');
write_4dfpifh('/data/cn4/evan/Task_parcellation/SteveSubjects/AvgNewtaskTimecourses.4dfp.img',size(finaldata,2),'bigendian');

map_vol_to_surface('/data/cn4/evan/Task_parcellation/SteveSubjects/AvgNewtaskTimecourses.4dfp.img','L')
    
