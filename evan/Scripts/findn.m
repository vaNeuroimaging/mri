function [sub v] = findn(A)
% DESCRIPTION:
%       "find" for n-D matrix
% INPUT:
%       A: n-D matrix
% OUTPUT: 
%       sub: the k x n index matrix (k is the number of the nonzeros)
%       v: a column or row vector v of the nonzero entries in A
% USAGE:
%       [sub ] = findn(A)
% AUTHOR:
%       Brian H. Hui (brianhui@alumni.usc.edu) 
%       PhD candidate, Electrical Engineering, University of Southern California
% HISTORY: 
%       03/18/09:  Created
%       
% TODO: 
%       
% SEE ALSO: 
%       FIND, IND2SUB
% TESTING CODE:
%{
if 0 % for Matlab 6.5 and below wich doesnot support %{ %} comments
    A = zeros(4,3,5,7);
    A(4,2,1,4)=1;
    A(2,3,4,6)=2;
    A(1,2,3,4)=3;
    tic;
    findn(A)
    toc;
end
%}
sA = size(A);
dA = length(sA);
B = reshape(A,prod(sA(1:end-1)),sA(end));

[sub(:,dA-1),sub(:,dA),v] = find(B);

for k = dA-1:-1:2
    [sub(:,k-1) sub(:,k)] = ind2sub([prod(sA(1:k-1)) sA(k)],sub(:,k));
end