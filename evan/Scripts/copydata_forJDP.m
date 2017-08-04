subjects = textread('/home/data/subjects/list_forJDP.txt','%s');

for s = 1:length(subjects)
    subject = subjects{s};
    mkdir(['/home/data/evan/for_JDP/' subject])
    for i=1:100
        file = ['/home/data/subjects/' subject '/preprocessed/restingstate_' num2str(i) '.nii.gz'];
        if exist(file)
            copyfile(file,['/home/data/evan/for_JDP/' subject '/'])
        end
    end
    
    T1s = dir(['/home/data/subjects/' subject '/raw/' subject '-1/oMPRageT1*.nii.gz']);
    for i = 1:length(T1s)
        file = ['/home/data/subjects/' subject '/raw/' subject '-1/' T1s(i).name];
        copyfile(file,['/home/data/evan/for_JDP/' subject '/'])
    end
end