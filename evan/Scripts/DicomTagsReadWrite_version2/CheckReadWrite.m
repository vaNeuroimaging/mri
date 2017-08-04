function err=CheckReadWrite(Filename)
% This function CheckReadWrite, first reads the tags of a dicom file into 
% a struct with ReadDicomElementList, and then write them back with
% WriteDicomElementList to a test-file, looking for differences.
%
% E = CheckReadWrite(filename)
%
% inputs,
%   filename : The dicom file name
% 
% outputs,
%   E : True if there are differences
%
%
% example,
%
%   CheckReadWrite('images\example2.dcm');
%
% Function is written by D.Kroon University of Twente (October 2010)

if(nargin<1)
    [fn, dn] = uigetfile('.dcm', 'Select a dicom file'); Filename=[dn fn];
end

FilenameTest='testreadwrite.dcm';
Elements = ReadDicomElementList(Filename);
WriteDicomElementList(Elements,FilenameTest);

fid = fopen(Filename, 'r'); bd1= fread(fid, inf,'uint8')'; fclose(fid);
fid = fopen(FilenameTest, 'r'); bd2= fread(fid, inf,'uint8')'; fclose(fid);
delete(FilenameTest);

if(length(bd1)>length(bd2)), bd2(end+1:length(bd1))=0; elseif(length(bd2)>length(bd1)), bd1(end+1:length(bd2))=0; end
bytediff=bd1~=bd2;
err=nnz(bytediff)>0;
if(err)
    disp('Failure : There are differences between the original and test file')
else
    disp('Succes : Original and test file are identical')
end


 
