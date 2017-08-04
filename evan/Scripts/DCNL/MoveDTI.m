subs = {'172','181','182','187','189','199','202','207','211','214','215','221','225','229','232','233','242','250','254','255','269','272','274'};
%'101','102','110','113','118','120','122','125','127','132','138','147','150','151','154','156','159','160','161','162',
%'166',
rootdir = '/fmri/data3/Evan/Gene-Rest-Nback/Data/';

for subject = 1:length(subs)
    
    subj = subs{subject};
    disp(subj)
    cd([rootdir subj]);
     mkdir('DTI');
     mkdir('DTI/run1');
     mkdir('DTI/run2');
     
     filestomove = dir(['/fmri/data3/Evan/DTI/Converted/*CV-' subj '_*']);
     
     copyfile(['/fmri/data3/Evan/DTI/Converted/' filestomove(1).name], [rootdir subj '/DTI/run1']);
     copyfile(['/fmri/data3/Evan/DTI/Converted/' filestomove(2).name], [rootdir subj '/DTI/run2']);
     
     copyfile(['/fmri/data3/phil/connectivity/dti/control/023/run1/bvals'],[rootdir subj '/DTI/run1']);
     copyfile(['/fmri/data3/phil/connectivity/dti/control/023/run1/bvecs'],[rootdir subj '/DTI/run1']);

     copyfile(['/fmri/data3/phil/connectivity/dti/control/023/run1/bvals'],[rootdir subj '/DTI/run2']);
     copyfile(['/fmri/data3/phil/connectivity/dti/control/023/run1/bvecs'],[rootdir subj '/DTI/run2']);

end