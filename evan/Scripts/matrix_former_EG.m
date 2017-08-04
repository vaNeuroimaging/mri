function [rmat] = matrix_former_TL(matrix,subjectA,subjectZ,dims_to_return,diagonalstatus,varargin)
%
% Name:matrix_former.m
% $Revision: 1.2 $
% $Date: 2011/03/08 20:30:12 $
%
% jdp 10/10/10
%
% This script returns processed matrices in service of graphtools. 
% 
% A matrix is passed in, and can be returned as a mean (2D) or 3D matrix,
% depending on the dimensions to return switch. Additionally, the diagonal
% will be taken out if takediagonalout says to.
% 
% USAGE: [rmat] = matrix_former(matrix,subjectA,subjectZ,dims_to_return,diagonalstatus,*#bootstrapsamples*)
% USAGE: [rmat] = matrix_former(mat,1,30,'3D','diagout',32)
% USAGE: [rmat] = matrix_former(mat,20,34,'2D','diagout')
% 
% SubjectA and SubjectZ are the bounds of the 3rd matrix dimension to
% consider. Set these to 1 and 1 if using a 2D matrix
% 
% dims_to_return can be '2D' or '3D'
% diagonalstatus can be 'diagout' or 'diagin'
% #bootstrapsamples is optional, for when you're bootstrapping
% 

d=size(matrix);
matrix=single(matrix);
if isempty(varargin) % no bootstrapping
    switch dims_to_return
        case '2D'
            rmat=mean(matrix(:,:,subjectA:subjectZ),3); % save RAM-make 2D
        case '3D'
            rmat=matrix(:,:,subjectA:subjectZ); % make 3D
        otherwise
            fprintf('Please select ''2D'' or ''3D'' as the dims_to_return\n');
    end 
else % if bootstrapping
    numsamples=varargin{1,1};
    subjectrange=subjectZ-subjectA+1;
    samples=ceil(rand(numsamples,1)*subjectrange); % pick samples with resampling
    samples=samples+subjectA-1; % shift to the appropriate position of the matrix
    switch dims_to_return
        case '2D'
            rmat=mean(matrix(:,:,samples),3); % save RAM-make 2D
        case '3D'
            rmat=matrix(:,:,samples); % make 3D
        otherwise
            fprintf('Please select ''2D'' or ''3D'' as the dims_to_return\n');
    end
end

% now remove the diagonal values, if desired
switch diagonalstatus
    case 'diagin'
    case 'diagout'
        for i = 1:size(rmat,1)
            rmat(i,i,:) = 0;
        end
%         for i=1:size(rmat,3)
%             tempmat=rmat(:,:,i);
%             tempmat(logical(eye(size(tempmat,1))))=0;
%             rmat(:,:,i)=tempmat;
%         end
    otherwise
        fprintf('Please use ''diagin'' or ''diagout'' for diagonalstatus\n');
end
