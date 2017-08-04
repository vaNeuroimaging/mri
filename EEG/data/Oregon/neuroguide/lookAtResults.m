function [ output_args ] = lookAtResults( input_args )
%LOOKATRESULTS Summary of this function goes here
%   Detailed explanation goes here
load('artifactRejectedFlankerCorrelationResults.mat');
index = find(strcmp({allCorrelations.measureName}, 'CAPS'));
capsData = allCorrelations(index);

scatter(capsData.input(:,2), capsData.output);

end

