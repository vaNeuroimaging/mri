function [ outputData, inputData ] = getPtsdOutputData( inputTable )
%GETFLANKEROUTPUTDATA Summary of this function goes here
%   Pass in eegdata as an input table.  Returns corresponding output
%   as well as modified inputTable (discarding EEG data with no
%   corresponding psychological variables)

%folder = fileparts(which('concatenateEegAndPsych'));
load('wahbehVariables.mat');

%track first or second session in t1ort2; rename Eeg data columns to match the
%vetmindData table
wahbehIds = ptsdData{:,1};
tableIds=inputTable.Properties.VariableNames(2:end);
t1ort2 = zeros(length(tableIds),1);
 for i = 1:length(tableIds)
     temp = tableIds{i};
     temp = temp(3:end);
%     temp = strrep(temp, 'POINT', '.');
%     t1ort2(i) = str2num(temp(5));
%     temp = temp(1:3);
%     temp = strcat('VM', temp);
     tableIds{i} = str2num(temp) - 100;
end

inputData = [];
outputData = table;
inputTemplate = table2array(inputTable(:,2:end))';
matches = cell(0,0);
match = zeros(1, length(tableIds));
for i=1:length(wahbehIds)
    for j = 1:length(tableIds)
        if(wahbehIds(i) == tableIds{j})
            matches{end+1} = wahbehIds(i);
            outputData(end+1,1:size(vetmindData,2))=vetmindData(i,:);
            pclScore = vetmindData{i,'PCL1'};
            if(t1ort2(i) == 2)
                pclScore = vetmindData{i,'PCL2'};
            end
            if(iscell(pclScore))
                pclScore = NaN;
            end
            outputData{end, size(vetmindData,2)+1} = pclScore;
            if(i == 1)
                outputData.Properties.VariableNames(1:size(vetmindData,2)) = vetmindData.Properties.VariableNames;
                outputData.Properties.VariableNames{size(vetmindData,2)+1} = 'PclCombined';
            end
            inputData(end+1,:)=inputTemplate(j,:);
            match(j) = 1;
        end
    end
end

end

