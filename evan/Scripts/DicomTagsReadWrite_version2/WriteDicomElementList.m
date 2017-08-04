function WriteDicomElementList(Elements,fname)
% This function WriteDicomElementList will write a struct with (raw) tags
% to a new dicom file.
%
%   WriteDicomElementList(Elements,filename)
%
% inputs,
%   Elements : A struct with raw tags from ReadDicomElementList, used are:
%             Elements(*).name
%             Elements(*).data
%             Elements(*).group
%             Elements(*).number
%             Elements(*).type
%             Elements(*).length
%             Elements(*).explicit
%
% example, Read&Write
%   Elements = ReadDicomElementList('images\example.dcm');
%   WriteDicomElementList(Elements,'output.dcm');
%
% example, Replace Value
%
%   Elements=ReadDicomElementList('images\example3.dcm');
%   i=structfind(Elements,'name','PatientsName');
%   Elements(i).data = ['Peter^Pan^-^-'];
%   Elements(i).length= length(Elements(i).data );
%   WriteDicomElementList(Elements,'output.dcm');
%   info=dicominfo('images\example3.dcm'); info.PatientName
%   info=dicominfo('output.dcm'); info.PatientName
%
% Function is written by D.Kroon University of Twente (October 2010)

% Display a file choose box if not provide as function input
if(nargin<2)
    [fn, dn] = uiputfile('.dcm', 'Choose a dicom file name'); fname=[dn fn];
end

% Open the dicom file
f=fopen(fname,'w', 'ieee-le');
if(f<0), 
    error('WriteDicomElementList:input',['could not open file' fname]);
end

if(strcmp(Elements(1).group,'0002'))
    % Write the Dicom header Prefix
    fwrite(f,zeros(1,128,'uint8'),'uint8');

    % Write a dicom identifier
    fwrite(f, 'DICM', 'char*1');
end

% Load Dicom Tag Library
functionname='WriteDicomElementList.m';
functiondir=which(functionname);
functiondir=functiondir(1:end-length(functionname));
load([functiondir 'Dictonary/DicomTagDictionary.mat']); 

% Write all tags
for i=1:length(Elements)
    WriteDicomElement(f,Elements(i),dcmdic);
end
fclose('all');

function WriteDicomElement(f,element,dcmdic)
fwrite(f,hex2dec(element.group(3:4)));
fwrite(f,hex2dec(element.group(1:2)));

fwrite(f,hex2dec(element.number(3:4)));
fwrite(f,hex2dec(element.number(1:2)));

if(element.length<0), element.length=2^32-1; end
        
if(element.explicit>0)
    fwrite(f,element.type,'char');
    if(element.explicit>1)
        switch(element.type)
        case {'OB','OW','UN','SQ'}
           fwrite(f, 0, 'uint16');
        end    
        fwrite(f,element.length,'uint32');
    else
        fwrite(f,element.length,'uint16');
    end
else
    fwrite(f,element.length,'uint32');
end

switch(element.type(1:2))
    case 'OB' %   OB -     |single trailing 0x00 to make even number of bytes. Transfer Syntax determines len
        fwrite(f, element.data, 'uint8');
    case 'SH' %   SH 16    |Short String. may be padded
        fwrite(f, element.data, 'char');
    case 'SQ' %   SQ  -    |Sequence of zero or more items
        fwrite(f, element.data, 'uint8');
    case 'UI' %   UI 64    |Unique Identifier (delimiter = .) 0-9 only, trailing space to make even #
        fwrite(f, element.data, 'char');
    case 'UL' %   UL 4     |Unsigned long integer
        fwrite(f,element.data,'ulong');
    case 'US' %   US 2     |Unsigned short integer (word)
        fwrite(f,element.data,'ushort');
    case 'SS' %   SS 2     |signed short integer (word)
        fwrite(f,element.data,'short');
    case 'SL' %   SL 4     |signed long integer
        fwrite(f,element.data,'long');
    case 'FL' %   FL 4     |Single precision floating pt number (float)
        fwrite(f,element.data,'float');
    case 'FD' %   FD 16    |Double precision floating pt number (double)
        fwrite(f,element.data,'double');
    case 'AE' %   AE 16    |Application Name
        fwrite(f, element.data, 'char');  
    case 'OF' %   OF -     |Other Float String. floats
        fwrite(f, element.data, 'char');  
    case 'OW' %   OW -     |Other Word String. words
        fwrite(f, element.data, 'uint8');  
    case 'PN' %   PN -     |Person's Name 64byte max per component. 5 components. delimiter = ^
        fwrite(f, element.data, 'char');  
    case 'ST' %   ST 1024  |Short Text of chars
        fwrite(f, element.data, 'char');  
    case 'TM' %   TM 16    |Time hhmmss.frac (or older format: hh:mm:ss.frac)
        fwrite(f, element.data, 'char');  
    case 'UT' %   UT -     |Unlimited Text. trailing spaces ignored
        fwrite(f, element.data, 'char');  
    case 'AS' %   AS 4     |Age String: nnnW or nnnM or nnnY
        fwrite(f, element.data, 'char');  
    case 'AT' %   AT 4     |Attribute Tag gggg,eeee
        fwrite(f, element.data, 'char');  
    case 'CS' %   CS 16    |Code String
        fwrite(f, element.data, 'char');  
    case 'DA' %   DA 8     |Date yyyymmdd (check for yyyy.mm.dd also and convert)
        fwrite(f, element.data, 'char');
    case 'DS' %   DS 16    |Decimal String may start with + or - and may be padded with l or t space
        fwrite(f, element.data, 'char');        
    case 'DT' %   DT 26    |Date Time YYYYMMDDHHMMSS.FFFFFF&ZZZZ (&ZZZ is optional & = + or -)
        fwrite(f, element.data, 'char');        
    case 'IS' %   IS 12    |Integer encoded as string. may be padded
        fwrite(f, element.data, 'char');       
    case 'LO' %   LO 64    |Character string. can be padded. cannot contain \ or any control chars except ESC
        fwrite(f, element.data, 'char');       
    case 'LT' %   LT 10240 |Long Text. Leading spaces are significant. trailing spaces arent
        fwrite(f, element.data, 'char');            
    otherwise
        fwrite(f, element.data, 'char');
end
