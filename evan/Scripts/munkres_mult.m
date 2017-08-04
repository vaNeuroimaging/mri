function [assignment1,assignment2,cost1, cost2] = munkres_mult(costMat,max_mult_assigns)

bigcostMat = repmat(costMat,1,max_mult_assigns);

[assignment,cost1] = munkres(bigcostMat);
assignment1 = mod(assignment,size(costMat,2));
assignment1(logical((mod(assignment,size(costMat,2))==0) .* (assignment>0))) = size(costMat,2);
%assignment1 = assignment1(1:size(costMat,1));

bigcostMat = repmat(costMat,max_mult_assigns,1);

[assignment,cost2] = munkres(bigcostMat');
assignment2 = mod(assignment,size(costMat,1));
assignment2(logical((mod(assignment,size(costMat,1))==0) .* (assignment>0))) = size(costMat,1);
%assignment2 = assignment2(1:size(costMat,2));

