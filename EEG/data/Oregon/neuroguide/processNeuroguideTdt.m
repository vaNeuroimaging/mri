function [ output_args ] = processNeuroguideTdt( outputVariableNames )
%PROCESSNEUROGUIDETDT Summary of this function goes here
%   Detailed explanation goes here

if(~exist('outputVariableNames'))
    outputVariableNames = {'CAPS', 'BDI1'};
end

load('neuroguideTdtANT.mat');
m = table2array(tab(:,2:end));
variables = tab{:,1};
subjectIds = tab.Properties.VariableNames;
subjectIds(1) = [];
for i = 1:length(subjectIds)
    id = subjectIds{i};
    id(1) = [];
    id = strrep(id, 'POINT', '.');
    subjectIds{i} = id;
end
for i = size(m,1):-1:1
    if(any(isnan(m(i,:))))
        m(i,:) = [];
        variables(i) = [];
    end
end

%load('/Users/Geoff/Documents/MATLAB/EEG/wahbehVariables.mat');
load('/home/gmay/Documents/MATLAB/wahbehVariables.mat');
cleanData = vetmindData;
for i = size(cleanData,2):-1:1
    dataPoint = cleanData{1,i};
    if(~isnumeric(dataPoint))
        cleanData(:,i) = [];
    end
end

targetIndices = NaN(1,length(outputVariableNames));
for i = 1:length(targetIndices)
    for j=1:length(cleanData.Properties.VariableNames)
        varName = cleanData.Properties.VariableNames{j};
        if(strcmp(outputVariableNames{i}, varName))
            targetIndices(i) = j;
        end
    end
end

cleanData = cleanData{:,targetIndices};

outputData = NaN(size(cleanData,2),size(m,2));
dataOwner = cell(1, length(subjectIds));
for i = 1:length(subjectIds)
    id1 = subjectIds{i};
    id1 = id1(1:3);
    for j = 1:size(cleanData,1)
        id2 = vetmindData{j,1};
        id2 = id2{1};
        id2 = id2(3:5);
        if(strcmp(id1,id2))
           outputData(:,i) = cleanData(j,:); 
           dataOwner{i} = id2;
        end
    end
end

inputData = m';
maxInput = max(max(inputData));
minInput = min(min(inputData));
maxOutput = max(max(outputData));
minOutput = min(min(outputData));

inputData = (inputData - minInput) ./ (maxInput-minInput);
outputData = (outputData - minOutput) ./ (maxOutput-minOutput);

outputData = outputData';
save('neuroguideNeuralTraining.mat', 'inputData', 'outputData');

[traIn, traOut, tesIn, tesOut, valIn, valOut] = assignToGroups(inputData,outputData, dataOwner, 0.1, 0.1);
clear inputData;

varNames = {'valOut', 'tesOut', 'traOut', 'valIn', 'tesIn', 'traIn'};
csvStamp = cell(1,length(varNames)); csvStamp(1:end) = {'.csv'};
npyStamp = cell(1,length(varNames)); npyStamp(1:end) = {'.npy'};
csvNames = strcat(varNames,csvStamp);
npyNames = strcat(varNames,npyStamp);

for i = 1:length(varNames)
  csvwrite(csvNames{i}, eval(varNames{i}));
  command = '';
  command = sprintf('%secho "import numpy as np";',command);
  command = sprintf('%secho "array = np.genfromtxt(''%s'',delimiter='','')";',command, csvNames{i});
  command = sprintf('%secho "print(np.shape(array))";',command);
  command = sprintf('%secho "np.save(''%s'', arr=array)";',command, npyNames{i});
  command = sprintf('(%s) | python', command);
  unix(command);
end


if(exist('outputVariableNames'))
    outputVars = outputVariableNames;
end
notes = sprintf('summary of neuroguide tdts %d\n   created %s\npatient ids: ', date);
for i = 1:length(subjectIds)
  %[folder file extension] = fileparts(allEvents(targetFileIndices(i)).filename);
  notes = [notes subjectIds{i} '; '];
end
notes = [notes sprintf('\noutputVariables: ')];
for i = 1:length(outputVars)
  notes = [notes outputVars{i} '; '];
end

notes = [notes sprintf('\nmin input: %f\n max input: %f\n min output: %f\nmax output:', ...
    minInput, maxInput, minOutput, maxOutput)];
% notes = [notes sprintf('\nvariable min Values: ')];
% for i = 1:length(minVals)
%   notes = [notes num2str(minVals(i)) '; '];
% end
% notes = [notes sprintf('\nvariable max Values: ')];
% for i = 1:length(maxVals)
%   notes = [notes num2str(maxVals(i)) '; '];
% end

channels = wahbehChannelLocs;
notes = [notes sprintf('\nEEG channel sequence: ')];
for i = 1:length(channels)
  notes = [notes channels{i} '; '];
end
notes = [notes sprintf('\ninput data are from each channel listed above, with the first %d set of numbers corresponding to the first sample from those channels', length(channels))];
notesFilename = 'generationNotes.txt';
fileId = fopen(notesFilename, 'w');
fprintf(fileId, notes);
fclose(fileId);


end

