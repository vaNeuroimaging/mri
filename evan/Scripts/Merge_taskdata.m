datafolder = '/data/cn3/steven/Task_fcMapping/';

files = dir([datafolder '*.img']);
files2 = struct2cell(files);
for i=1:368
subnames{i} = files2{1,i}(1:8);
end
subjects = unique(subnames);

for subject = 1:length(subjects)
    
    underscoreindices = strfind(subjects{subject},'_');
    
    if length(underscoreindices)>1
        shortsubname = subjects{subject}(2:underscoreindices(2)-1);
    else
        shortsubname = subjects{subject}(2:end);
    end
    
    %allshortsubnames{subject,1} = shortsubname;
    
    files = dir([datafolder subjects{subject} '*.img']);
    
    %string = 'paste_4dfp';
    
    listname = '/data/cn4/evan/Task_parcellation/SteveSubjects/Templist.txt';
    delete([listname]);
    fid = fopen([listname],'at'); %open the output file for writing
    fclose(fid);
    dlmwrite([listname],' ','-append');
    
    for file = 1:length(files)
        
        %string = [string ' ' datafolder files(file).name];
        string = [datafolder files(file).name];
        dlmwrite([listname],string,'-append','delimiter','');
        
    end
    
    string = ['paste_4dfp ' listname ' /data/cn4/evan/Task_parcellation/SteveSubjects/' shortsubname '_Mergefile.4dfp.img -a -p16'];
    
    system(string)
    
    map_vol_to_surface(['/data/cn4/evan/Task_parcellation/SteveSubjects/' shortsubname '_Mergefile.4dfp.img'],'L')
    
end