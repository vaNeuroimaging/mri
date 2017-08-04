function CreateDicomDict(filename)
% This function CreateDicomDict, will read a dicom-dict.txt file to
% generate a Matlab dicom-tag lookup file.
%  
%  CreateDicomDict(filename);
%
% example,
%   filename='dicom-dict-large.txt';
%   CreateDicomDict(filename);
%
%
% Function is written by D.Kroon University of Twente (October 2010)

dcmdic=struct();
tagdata=struct();
[dcm_tag dcm_type dcm_name dcm_length]=textread(filename,'%s%s%s%s');
for num=1:length(dcm_tag);
    tag=dcm_tag{num}; 
    if((tag(1)=='(')&&(length(tag)>9)),
        tag=['tag' tag(2:5) '' tag(7:10)];
        tagdata.type=dcm_type{num};
        tagdata.name=dcm_name{num};
        tagdata.length=dcm_length{num};
        dcmdic.(tag)=tagdata;
    end
end
save('DicomTagDictionary.mat','dcmdic');