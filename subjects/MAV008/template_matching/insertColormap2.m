function [ output_args ] = insertColormap2( filename )
%INSERTCOLORMAPINNIFTI Summary of this function goes here
%   Detailed explanation goes here

if(~exist('filename'))
    filename = 'restingstate_templatematch_dice_kden0.05.dtseries.nii';
end

sampleTarget= sprintf('<Matrix>\n');
insertionFilename = 'insertion.txt';
if(~exist(insertionFilename))    
    sampleFilename = '/home/data/subjects/MAV006/template_matching/restingstate_templatematch_dice_kden0.05.dtseries.nii';
    samplefileId = fopen(sampleFilename);
    samplefileBytes = fread(samplefileId);
    fclose(samplefileId);
    samplefileChar = char(fileBytes)';
    sampleStart = strfind(samplefileChar, sampleTarget) + length(sampleTarget);
    sampleTarget2 = sprintf('</MetaData>\n');
    sampleEnd = strfind(samplefileChar, sampleTarget2) + length(sampleTarget2);
    insertionText = samplefileChar(sampleStart:sampleEnd);
    fileId = fopen(insertionFilename, 'w');
    fprintf(fileId, insertionText);
else
    fileId = fopen(insertionFilename);
    insertionBytes = fread(fileId);
    fclose(fileId);
    insertionText = char(insertionBytes)';
end

fileId = fopen(filename);
fileBytes = fread(fileId, inf);
fclose(fileId);
fileChar = char(fileBytes)';

sampleStart = strfind(fileChar, sampleTarget) + length(sampleTarget);
pre = fileChar(1:sampleStart);
post = fileChar(sampleStart:end);
newfileChar = strcat(pre, insertionText, post);


%insertstring = sprintf('\n\t<MetaData>\n\t\t<MD>\n\t\t\t<Name>PaletteColorMapping</Name>\n\t\t\t<Value>&lt;PaletteColorMapping Version=&quot;1&quot;&gt;\n\t&lt;ScaleMode&gt;MODE_USER_SCALE&lt;/ScaleMode&gt;\n\t&lt;AutoScalePercentageValues&gt;98.000000 2.000000 2.000000 98.000000&lt;/AutoScalePercentageValues&gt;\n\t&lt;AutoScaleAbsolutePercentageValues&gt;2.000000 98.000000&lt;/AutoScaleAbsolutePercentageValues&gt;\n\t&lt;UserScaleValues&gt;-100.000000 0.000000 1.000000 18.000000&lt;/UserScaleValues&gt;\n\t&lt;PaletteName&gt;power_surf&lt;/PaletteName&gt;\n\t&lt;InterpolatePalette&gt;true&lt;/InterpolatePalette&gt;\n\t&lt;DisplayPositiveData&gt;true&lt;/DisplayPositiveData&gt;\n\t&lt;DisplayZeroData&gt;false&lt;/DisplayZeroData&gt;\n\t&lt;DisplayNegativeData&gt;false&lt;/DisplayNegativeData&gt;\n\t&lt;ThresholdTest&gt;THRESHOLD_TEST_SHOW_OUTSIDE&lt;/ThresholdTest&gt;\n\t&lt;ThresholdType&gt;THRESHOLD_TYPE_OFF&lt;/ThresholdType&gt;\n\t&lt;ThresholdFailureInGreen&gt;false&lt;/ThresholdFailureInGreen&gt;\n\t&lt;ThresholdNormalValues&gt;-1.000000 1.000000&lt;/ThresholdNormalValues&gt;\n\t&lt;ThresholdMappedValues&gt;-1.000000 1.000000&lt;/ThresholdMappedValues&gt;\n\t&lt;ThresholdMappedAvgAreaValues&gt;-1.000000 1.000000&lt;/ThresholdMappedAvgAreaValues&gt;\n\t&lt;ThresholdDataName&gt;&lt;/ThresholdDataName&gt;\n\t&lt;ThresholdRangeMode&gt;PALETTE_THRESHOLD_RANGE_MODE_MAP&lt;/ThresholdRangeMode&gt;\n\t&lt;ThresholdLowHighLinked&gt;false&lt;/ThresholdLowHighLinked&gt;\n&lt;/PaletteColorMapping&gt;\n</Value>\n\t\t</MD>\n\t</MetaData>');
% insertstring = strrep(insertstring, '&lt;', '<');
% insertstring = strrep(insertstring, '&gt;', '>');
% insertstring = strrep(insertstring, '&quot;', '"');

%findstring = '<Matrix>';
% findstringlocation = strfind(fileChar,findstring);
% 
% newfileChar = [fileChar(1:(findstringlocation+length(findstring))) insertstring fileChar((findstringlocation+length(findstring)+1):end)];

%newfileChar = strrep(fileChar, findstring, [findstring insertstring]);


%newfileChar = fileChar;
writeBytes = single(newfileChar');

testName = strcat('autoColor2',filename);
fileId = fopen(testName,'w');
fwrite(fileId, writeBytes);
fclose(fileId);
end

