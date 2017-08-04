excludefile = '/home/data/subjects/QC/MAV_QC_031517.xlsx';
[NUM,TXT,RAW]=xlsread(excludefile);
subinds = [];
for i = 1:length(TXT)
    if ~isempty(TXT{i})
        subinds = [subinds i];
    end
end
%subinds = find(~isempty(TXT(:,1)));

runtrackerfile = '/home/data/subjects/QC/Allruns_toexclude.txt';
delete(runtrackerfile);
fid = fopen(runtrackerfile,'at');
fclose(fid);

for subnum = 1:length(subinds)
    subject = TXT{subinds(subnum),1};
    if subnum<length(subinds)
        excludenums = NUM(subinds(subnum) : (subinds(subnum+1)-1));
    else
        excludenums = NUM(subinds(subnum) : end);
    end
    
    [filenames, boldrunnames] = textread(['/home/data/subjects/' subject '/preprocessed/BOLD_run_tracker.txt'],'%s%s');
    
    for excludenum = excludenums(:)'
        excludefile = ['/home/data/subjects/' subject '/preprocessed//RSFC_' num2str(excludenum) '.nii.gz'];
        origfile_toexclude = filenames{strcmp(excludefile,boldrunnames)};
        dlmwrite(runtrackerfile,origfile_toexclude,'-append','delimiter','')
    end
end