% studies = {'MAV','MBG','ROB'};
% potentialsubjects = dir('/home/data/subjects/');
% subjects = cell(0,1);
% for p = 1:length(potentialsubjects)
%     if (length(potentialsubjects(p).name) > 2) && (any(strcmp(potentialsubjects(p).name(1:3),studies))) && exist(['/home/data/subjects/' potentialsubjects(p).name '/freesurfer/label'])
%         subjects(end+1,1) = {potentialsubjects(p).name};
%         %mkdir(['/home/data/Analysis/ENIGMA/' potentialsubjects(p).name '_freesurfer'])
%         system(['ln -s /home/data/subjects/' potentialsubjects(p).name '/freesurfer /home/data/Analysis/ENIGMA/' potentialsubjects(p).name '_freesurfer'])
%     end
% end

subjects = '/home/data/subjects/processing_list_060717.txt';
if iscell(subjects)

elseif exist(subjects,'file')
    
    subjects = textread(subjects,'%s');
    
elseif ischar(subjects)
    
    subjects = {subjects};
    
end

for s = 1:length(subjects)
    
    system(['ln -s /home/data/subjects/' subjects{s} '/freesurfer /home/data/Analysis/ENIGMA/' subjects{s}])
    
end
    


