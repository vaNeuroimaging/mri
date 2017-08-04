% POSTPAD.M
% %POSTPAD   Pads or truncates a vector x to a specified length L.
% %
% %   This is a simpler version of the file distributed with octave.
% %   The function does not check input parameters.
% 
% ---------------------------------------------------------------
%
% COPYRIGHT : (c) NUHAG, Dept.Math., University of Vienna, AUSTRIA
%             http://nuhag.eu/
%             Permission is granted to modify and re-distribute this
%             code in any manner as long as this notice is preserved.
%             All standard disclaimers apply.
%
function x = postpad (x, L)
%POSTPAD   Pads or truncates a vector x to a specified length L.
%
%   This is a simpler version of the file distributed with octave.
%   The function does not check input parameters.


xl=size(x,1);

if xl==L
  % Do nothing, bails out as quickly as possible.
  return;
end;

xw=size(x,2);

if ndims(x)>2
  error('Postpad of multidim not done yet.');
end;

if xl<L
  x=[x; zeros(L-xl,xw)];
else
  x=x(1:L,:);
end;
