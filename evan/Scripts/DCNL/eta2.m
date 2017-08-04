function [ e2 ] = eta2(a,b)
% Compute the eta-squared coefficient between two signals.
% Eta-squared is the fraction of the variance in one signal
% accounted for by variance in a second signal.
%
% FORMAT [ e2 ] = eta2(a, b)
%
% a  - First signal, any dimensionality.
% b  - Second signal, must be of same dimensionality as 'a'.
% e2 - Eta-squared coefficient between 'a' and 'b'.
%
% If an error is generated (e.g., a and b are not of the same dimensionality), this
% function returns NaN.  Otherwise, returns the eta-squared coefficient between 'a' and 'b'.
% Reference: see equation 1 in Cohen et al., Neuroimage 41(1):45-57 (2008)
% Variable names in this script follow the naming conventions in the paper.
% 

% Check dimensionality.
e2    = NaN; % Initialize
sizeA = size(a);
sizeB = size(b);
if ( length(sizeA) ~= length(sizeB) )
    return
end

% Dimensionality the same, now check dimensions.
sameDims = ( sizeA == sizeB );
if ( ~ all(sameDims) )
    return
end

% Compute within-element mean 'm' and grand mean 'M'.
% The colon operator converts arrays of arbitrary dimension to column vectors.
m = ( a(:) + b(:) ) / 2;
M = mean(m);

% Compute numerator.
% The colon operator converts arrays of arbitrary dimension to column vectors.
aTerm     = a(:) - m(:);
aTerm     = aTerm' * aTerm;
bTerm     = b(:) - m(:);
bTerm     = bTerm' * bTerm;
SS_Within = sum(aTerm(:) + bTerm(:));

% Compute denominator.
aTerm    = a(:) - M;
aTerm    = aTerm' * aTerm;
bTerm    = b(:) - M;
bTerm    = bTerm' * bTerm;
SS_Total = sum(aTerm(:) + bTerm(:));

% Compute eta-squared.
if ( SS_Total == 0 )
    e2 = Inf;
    return
else
    e2 = 1 - ( SS_Within / SS_Total );
end
