function [coef,pval] = partialcorr(x,z,varargin)
%PARTIALCORR Linear or rank partial correlation coefficients.
%   RHO = PARTIALCORR(X,Z) returns the sample linear partial correlation
%   coefficients between pairs of variables in X, controlling for the
%   variables in Z.  X is an N-by-P matrix, and Z an N-by-Q matrix, with rows
%   corresponding to observations, and columns corresponding to variables.  R
%   is a symmetric P-by-P matrix.
%
%   RHO = PARTIALCORR(X,Y,Z) returns the sample linear partial correlation
%   coefficients between pairs of variables between X and Y, controlling for
%   the variables in Z.  X is an N-by-P1 matrix, Y an N-by-P2 matrix, and Z an
%   N-by-Q matrix, with rows corresponding to observations, and columns
%   corresponding to variables.  R is a P1-by-P2 matrix, where the (I,J)-th
%   entry is the sample linear partial correlation between the I-th column in
%   X and the J-th column in Y.
%
%   If the covariance matrix of [X,Z] is S = [S11 S12; S12' S22], then the
%   partial correlation matrix of X, controlling for Z, can be defined
%   formally as a normalized version of the covariance matrix
%
%      S_XZ = S11 - S12*inv(S22)*S12'.
%
%   [RHO,PVAL] = PARTIALCORR(...) also returns PVAL, a matrix of p-values for
%   testing the hypothesis of no partial correlation against the alternative
%   that there is a non-zero partial correlation.  Each element of PVAL is the
%   p-value for the corresponding element of RHO.  If PVAL(i,j) is small, say
%   less than 0.05, then the partial correlation RHO(i,j) is significantly
%   different from zero.
%
%   [...] = PARTIALCORR(...,'PARAM1',VAL1,'PARAM2',VAL2,...) specifies
%   additional parameters and their values.  Valid parameters are the
%   following:
%
%        Parameter  Value
%         'type'    'Pearson' (the default) to compute Pearson (linear)
%                   partial correlations or 'Spearman' to compute Spearman
%                   (rank) partial correlations.
%         'rows'    'all' (default) to use all rows regardless of missing
%                   values (NaNs), 'complete' to use only rows with no
%                   missing values, or 'pairwise' to compute RHO(i,j) using
%                   rows with no missing values in column i or j.
%         'tail'    The alternative hypothesis against which to compute
%                   p-values for testing the hypothesis of no partial
%                   correlation.  Choices are:
%                      TAIL         Alternative Hypothesis
%                   ---------------------------------------------------
%                     'both'     correlation is not zero (the default)
%                     'right'    correlation is greater than zero
%                     'left'     correlation is less than zero
%
%   The 'pairwise' option for the 'rows' parameter can produce RHO that is not
%   positive definite.  The 'complete' option always produces a positive
%   definite RHO, but when data are missing, the estimates will in generally
%   be based on fewer observations.
%
%   PARTIALCORR computes p-values for linear and rank partial correlations
%   using a Student's t distribution for a transformation of the correlation.
%   This is exact for linear partial correlation when X and Z are normal, but
%   is a large-sample approximation otherwise.
%
%   See also CORR, CORRCOEF, TIEDRANK.

%   References:

%   Copyright 1984-2009 The MathWorks, Inc.
%   $Revision: 1.1.6.2 $  $Date: 2009/05/07 18:31:56 $

%   Partial correlation for X, controlling for Z, can be computed by
%   normalizing the full covariance matrix S_XZ = S11 - S12*inv(S22)*S12'.
%   However, PARTIALCORR instead computes it as the correlation of the
%   residuals from a regression of X on Z for linear partial correlation, or
%   from a regression of the ranks of X on the ranks of Z for rank partial
%   correlation.
%
%   An equivalent recursive definition in terms of the individual full
%   and partial correlation coefficients is
%      rxy_z = (rxy - rxz*ryx) / sqrt(1-rxz^2)*(1-ryz^2))
%      rxy_zw = (rxy_z - rxw_z*ryw_z) / sqrt(1-rxw_z^2)*(1-ryw_z^2))
%   etc.

if nargin<2
    error('stats:partialcorr:NotEnoughInputs', 'Not enough input arguments.');
end
if length(varargin) > 0 && isnumeric(varargin{1})
    % Two separate x and y matrices, make into a single matrix
    y = z;
    if size(x,1)~=size(y,1)
        error('stats:partialcorr:InputSizeMismatch', 'X and Y must have the same number of rows.');
    elseif ndims(x) > 2 || ndims(y) > 2
        error('stats:partialcorr:InputsMustBeMatrices', 'X and Y must be matrices.');
    end
    crossCorr = true;
    dx = size(x,2); dy = size(y,2);
    x = [x y];
    
    z = varargin{1};
    varargin(1) = [];
    
else
    if ndims(x) > 2
        error('stats:partialcorr:InputsMustBeMatrices', 'X must be a matrix.');
    end
    crossCorr = false;
    dx = size(x,2);
end
[n,d] = size(x);

if ndims(z) > 2
    error('stats:partialcorr:InputsMustBeMatrices', 'Z must be a matrix.');
end
if size(z,1) ~= n
    error('stats:partialcorr:InputSizeMismatch', 'X and Z must have the same number of rows.');
end
[n,dz] = size(z);

outClass = superiorfloat(x,z);

pnames = {'type'  'rows' 'tail'};
dflts  = {'p'     'a'    'both'};
[errid,errmsg,type,rows,tail] = internal.stats.getargs(pnames,dflts,varargin{:});
if ~isempty(errid)
    error(sprintf('stats:partialcorr:%s',errid),errmsg);
end

% Validate the rows parameter.
rowsChoices = {'all' 'complete' 'pairwise'};
if ischar(rows)
    i = strmatch(lower(rows),rowsChoices);
    if isscalar(i)
        rows = rowsChoices{i}(1);
    elseif isempty(i)
        error('stats:partialcorr:UnknownRows', ...
              'The ''rows'' parameter value must be ''all'', ''complete'', or ''pairwise''.');
    end
else
    error('stats:partialcorr:UnknownRows', ...
          'The ''rows'' parameter value must be ''all'', ''complete'', or ''pairwise''.');
end

% Validate the type parameter.
typeChoices = {'pearson' 'kendall' 'spearman'};
if ischar(type)
    i = strmatch(lower(type),typeChoices);
    if isscalar(i)
        type = typeChoices{i}(1);
    elseif isempty(i)
        error('stats:partialcorr:UnknownType', ...
              'The ''type'' parameter value must be ''Pearson'' or ''Spearman''.');
    end
    if type == 'k'
        error('stats:partialcorr:Kendall', ...
              'Cannot compute Kendall''s partial rank correlation.');
    end
else
    error('stats:partialcorr:UnknownType', ...
          'The ''type'' parameter value must be ''Pearson'' or ''Spearman''.');
end

% Validate the tail parameter.
tailChoices = {'left','both','right'};
if ischar(tail) && (size(tail,1)==1)
    i = find(strncmpi(tail,tailChoices,length(tail)));
    if isempty(i)
        i = find(strncmpi(tail,{'lt','ne','gt'},length(tail)));
    end
    if isscalar(i)
        tail = tailChoices{i}(1);
    elseif isempty(i)
        error('stats:partialcorr:UnknownTail', ...
              'TAIL must be one of the strings ''both'', ''right'', or ''left''.');
    end
else
    error('stats:partialcorr:UnknownTail', ...
          'TAIL must be one of the strings ''both'', ''right'', or ''left''.');
end

% Regress X on Z, and compute the correlation of the residuals.  Works
% even when the full (unconditional) correlation matrix of [X Z] would be
% indefinite due to pairwise missing data removal.
switch rows
case {'a' 'c'} % 'all' 'complete'
    if rows == 'c'
        notnans = ~any(isnan(x),2) & ~any(isnan(z),2);
        x = x(notnans,:);
        z = z(notnans,:);
        n = size(x,1);
    end
    if type == 's'
        x = tiedrank(x);
        z = tiedrank(z);
    end
    z1 = [ones(n,1) z];
 
    gx = gsingle(x); gz1 = gsingle(z1); gn = gsingle(n); gdz = gsingle(dz);
    resid = gzeros(size(z1,1),size(x,2));
    resid = gx - gz1*(gz1 \ gx);
    tol = max(gn,gdz)*eps(class(x))*sqrt(sum(abs(gx).^2,1));
    resid(:,sqrt(sum(abs(resid).^2,1)) < tol) = 0;
    resid = double(resid);
    clear gx gz1 gn gdz tol;
    
    if crossCorr

        first = gdouble(resid(:,1:dx));
        all = gdouble(resid(:,dx+1:dx+dy));
        coef = gzeros(1,size(all,2));
        
        gfor k = 1:length(coef)
            coef(k) = corr2(first, all(:,k));
        gend
        coef = double(coef);
        clear first all
        
        %coef = corr(resid(:,1:dx),resid(:,dx+1:dx+dy),'type','pearson');
    else
        
        coef = corr(resid,'type','pearson');
    end

case 'p' % 'pairwise'
    % Compute correlations on each pair of columns in X, with complete
    % row removal.
    znotnans = ~any(isnan(z),2);
    x = x(znotnans,:);
    z = z(znotnans,:);
    coef = zeros(d,outClass);
    n = zeros(size(coef));
    isnanx = isnan(x);
    for i = 1:dx
        % For cross correlation, only do the x:y cross terms. For
        % autocorrelation, only do the lower triangle and diagonal.
        if crossCorr
            j0 = dx+1; j1 = dx+dy;
        else
            j0 = 1; j1 = i;
        end
        for j = j0:j1
            notnans = ~(isnanx(:,i) | isnanx(:,j));
            xx = x(notnans,[i j]);
            zz = z(notnans,:);
            nn = size(xx,1);
            if type == 's'
                xx = tiedrank(xx);
                zz = tiedrank(zz);
            end
            z1 = [ones(nn,1) zz];
            resid = xx - z1*(z1 \ xx);

            tol = max(nn,dz)*eps(class(xx))*sqrt(sum(abs(xx).^2,1));
            resid(:,sqrt(sum(abs(resid).^2,1)) < tol) = 0;

            coef(i,j) = sum(prod(resid,2)) ./ prod(sqrt(sum(abs(resid).^2,1)),2);
            n(i,j) = nn;
        end
    end
    
    if ~crossCorr
        % Force a one on the diagonal, but preserve NaNs.
        ii = find(~isnan(diag(coef)));
        coef((ii-1)*d+ii) = 1;
        
        % Reflect to lower triangle.
        coef = tril(coef) + tril(coef,-1)';
        n = tril(n) + tril(n,-1)';
    else
        % Keep only the x:y cross term elements.
        coef = coef(1:dx,dx+1:dx+dy);
        n = n(1:dx,dx+1:dx+dy);
    end
end

if nargout > 1
    df = max(n - dz - 2,0); % a matrix for 'pairwise'
    t = sign(coef) .* Inf;
    k = (abs(coef) < 1);
    t(k) = coef(k) ./ sqrt(1-coef(k).^2);
    t = sqrt(df).*t;
    switch tail
    case 'b' % 'both or 'ne'
        pval = 2*tcdf(-abs(t),df);
    case 'r' % 'right' or 'gt'
        pval = tcdf(-t,df);
    case 'l' % 'left or 'lt'
        pval = tcdf(t,df);
    end
end

