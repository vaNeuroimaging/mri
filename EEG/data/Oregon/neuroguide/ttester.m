function [ output_args ] = ttester( input_args )
%TTESTER Summary of this function goes here
%   Detailed explanation goes here

measureName = 'PclCombined';

%flanker = load('/Users/Geoff/Documents/MATLAB/neuralNetPrep/neuroguide/neuroguideFlankerTraining.mat');
%wahbeh = load('correctedWahbehNeuroguide');
%load('neuroguideTdtPtsd.mat');
if(false)%PTSD
    [outputTable, inputData] = getPtsdOutputData(data);
    outputData = outputTable{:,measureName};
else%flanker
    %neural = load('/Users/Geoff/Documents/MATLAB/neuralNetPrep/neuroguide/neuroguideNeuralTraining.mat');
    neural = load('artifactRejectedFlankerEeg.mat');
%     flankerData = table;
%     tonesData = table;
%     flankerData(:,1) = neural.data(:,1);
%     tonesData(:,1) = neural.data(:,1);
%     
%     for i = 2:size(neural.data,2)
%         columnName = neural.data.Properties.VariableNames{i};
%         if(strcmp('ANT', columnName(end-2:end)))
%             flankerData{:,end+1}=neural.data{:,i};
%             flankerData.Properties.VariableNames{end} = neural.data.Properties.VariableNames{i};
%         else
%             tonesData{:,end+1}=neural.data{:,i};
%             tonesData.Properties.VariableNames{end} = neural.data.Properties.VariableNames{i};
%         end
%     end
    
    %wahbeh.data = neural.data;
    wahbeh.data = neural.data(:, 49:94);
    
    
    [outputTable, inputData] = getFlankerOutputData(wahbeh.data);
end
% data = tab;
% outputData = outputData(:,1);
% data = inputData;




outputVariableNames = outputTable.Properties.VariableNames;
% inputMatrix = table2array(wahbeh.data(:,2:end));
% inputData= inputMatrix';
allCorrelations = [];
for outputVariableNumber = 1:length(outputVariableNames)
    inputVariableNames = wahbeh.data{:,1};
    rhoArray = NaN(size(inputData, 2),1);
    pArray = NaN(size(inputData, 2),1);
    inputIndices = NaN(size(inputData,2),1);
    tic;
    measureName = outputVariableNames{outputVariableNumber};
    disp(sprintf('%d of %d (%s)', outputVariableNumber, length(outputVariableNames), measureName));
    output = outputTable{:,measureName};
    if(isnumeric(output))

        
        for i = size(inputData,2):-1:1
            %for each eeg value, find the correlation with CAPS
            matrix = [inputData(:,i), output];
            remove = any(isnan(matrix),2);
            matrix(remove,:) = [];
            if(length(matrix) == 0)
                inputVariableNames(i) = [];
                rhoArray(i) = [];
                pArray(i) = [];
                inputIndices(i) = [];
            else
                [rho, p] = corr(matrix);
                rhoArray(i) = rho(2,1);
                pArray(i) = p(2,1);
                inputIndices(i) = i;
            end
            %    [rhoArray(i), pArray(i)] = corr(inputData(:,i), outputData(:,1));
        end
        %correctedPArray = pArray .* i;
        correctedPArray = NaN(1,length(pArray));
        [sortedPValues, pIndexes] = sort(pArray);
        for i = 1:length(pArray)
            gentleBonferroni = length(pArray) - i + 1;
            index = pIndexes(i);
            correctedPArray(index) = pArray(index) * gentleBonferroni;
        end
        stillSignificant = correctedPArray < 0.05;
        sigIndexes = find(stillSignificant);
        
        alphaCorrected = pArray * 171;
        kindOfSignificant = alphaCorrected < 0.05;
        kindofSignIndexes = find(kindOfSignificant);
        
        if(length(sigIndexes > 0))
            correlation.sigIndexes = sigIndexes;
            correlation.sigVariables = inputVariableNames(sigIndexes);
            disp(inputVariableNames(sigIndexes));
            correlation.pValues = pArray(sigIndexes);
            correlation.rValues = rhoArray(sigIndexes);
            correlation.output = outputTable{:,measureName};
            correlation.input = inputData(:,sigIndexes);
            correlation.measureName = measureName;
            if(length(allCorrelations) > 0)
                allCorrelations(end+1) = correlation;
            else
                allCorrelations = correlation;
            end
            
        end
    end
    
    if(false)%display sum
        close all;
        sumVector = zeros(size(inputData,1), 1);
        for i = 1:length(sigIndexes)
            index = sigIndexes(i);
            matrix = [inputData(:,inputIndices(index)), outputTable{:,measureName}];
            normalizedInput= inputData(:,inputIndices(index));
            normalizedInput = normalizedInput - mean(normalizedInput);
            normalizedInput = normalizedInput ./ std(normalizedInput);
            sumVector = sumVector + normalizedInput;
            figure;
            scatter(matrix(:,1),matrix(:,2));
            title(inputVariableNames(index));
        end
        %     matrix = [sumVector, outputTable{:,measureName}];
        %     remove = any(isnan(matrix),2);
        %     matrix(remove,:) = [];
        %     [rho, p] = corr(matrix);
        %     figure;
        %     scatter(matrix(:,1),matrix(:,2));
        %     title('composite');
        %
        %     tilefigs;
        %     end
        % rhoArray(sigIndexes)
    end
toc;
%
% if(false)
%     %remove outliers
%     matrix = [sumVector(:), outputData(:,1)];
%     remove = any(isnan(matrix),2);
%     matrix(remove,:) = [];
%     %tempInput = inputData(:,i);
%     numberToRemove = 5;
%     patientRemovalList = NaN(1,numberToRemove);
%     for i = 1:numberToRemove
%         removeIndex = find(matrix(:,1)==min(matrix(:,1)));
%         patientRemovalList(i) = removeIndex;
%         matrix(removeIndex, :) = [];
%         %     tempInput(removeIndex) = [];
%     end
% end
% 
% close all;
% outlierMatrix = [];
% for i = 1:length(sigIndexes)
% %     patientRemovalList = NaN(1,numberToRemove);
% %     for j = 1:numberToRemove
% %         removeIndex = find(matrix(:,1)==min(matrix(:,1)));
% %         patientRemovalList(j) = removeIndex;
% %         disp(removeIndex);
% %         matrix(removeIndex, :) = [];
% %         %     tempInput(removeIndex) = [];
% %     end
% %     outlierMatrix(end+1,:) = patientRemovalList;
% %     %     index = sigIndexes(i);
% %     %     matrix = [inputData(:,inputIndices(index)), outputData(:,1)];
% %     %     remove = any(isnan(matrix),2);
% %     %     matrix(remove,:) = [];
% %     %     matrix(patientRemovalList,:) = [];
% %     %
%     figure;
%     scatter(matrix(:,1),matrix(:,2));
%     title(variableNames(index));
% end
% tilefigs;
% 
% 
% 
% 
% [rho, p] = corr(matrix);
end

for i = 1:length(allCorrelations)
    correlation = allCorrelations(i);
    for j= 1:length(allCorrelations(i).sigIndexes)
        
    end
end

save('artifactRejectedTruncatedFlankerCorrelationResultsBlock2.mat', 'allCorrelations');

