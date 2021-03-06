%% Copyright (C) 1999 Paul Kienzle
%%
%% This program is free software; you can redistribute it and/or modify
%% it under the terms of the GNU General Public License as published by
%% the Free Software Foundation; either version 2 of the License, or
%% (at your option) any later version.
%%
%% This program is distributed in the hope that it will be useful,
%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%% GNU General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with this program; if not, write to the Free Software
%% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
%
%% Generate a butterworth filter.
%% 
%% [b,a] = butter(n, Wc)
%%    low pass filter with cutoff pi*Wc radians
%%
%% [b,a] = butter(n, Wc, 'high')
%%    high pass filter with cutoff pi*Wc radians
%%
%% [b,a] = butter(n, [Wl, Wh])
%%    band pass filter with edges pi*Wl and pi*Wh radians
%%
%% [b,a] = butter(n, [Wl, Wh], 'stop')
%%    band reject filter with edges pi*Wl and pi*Wh radians
%%
%% [z,p,g] = butter(...)
%%    return filter as zero-pole-gain rather than coefficients of the
%%    numerator and denominator polynomials.
%% 
%% References: 
%%
%% Proakis & Manolakis (1992). Digital Signal Processing. New York:
%% Macmillan Publishing Company.

%% Author: pkienzle@cs.indiana.edu

function [Zz, Zp, Zg] = butter(n, W, stype)

  if (nargin>3 || nargin<2) || (nargout>3 || nargout<2)
    usage ('[b, a] or [z, p, g] = butter (n, W, [, ''ftype''])');
  end

  %% interpret the input parameters
  if ~((length(n)==1 && n == round(n) && n > 0))
    error ('butter: filter order n must be a positive integer');
  end

  stop = nargin==3;
  if stop && (strcmp(stype, 'high') || strcmp(stype, 'stop'))
    error ('butter: ftype must be ''high'' or ''stop');
  end

  [r, c]=size(W);
  if ~((length(W)<=2 && (r==1 || c==1)))
    error ('butter: frequency must be given as w0 or [w0, w1]');
  elseif ~(all(W >= 0 & W <= 1))
    error ('butter: critical frequencies must be in (0 1)');
  elseif ~((length(W)==1 || length(W) == 2))
    error ('butter: only one filter band allowed');
  elseif ~(length(W)==2 && (W(1) < W(2)))
    error ('butter: first band edge must be smaller than second');
  end

  %% Prewarp to the band edges to s plane
  T = 2;       % sampling frequency of 2 Hz
  Ws = 2/T*tan(pi*W/T);

  %% Generate splane poles for the prototype butterworth filter
  %% source: Kuc
  C = 1; % default cutoff frequency
  Sp = C*exp(1i*pi*(2*[1:n] + n - 1)/(2*n));
  Sz = [];
  Sg = C^n;

  %% splane frequency transform
  [Sz, Sp, Sg] = sftrans(Sz, Sp, Sg, Ws, stop);

  %% Use bilinear transform to convert poles to the z plane
  [Zz, Zp, Zg] = bilinear(Sz, Sp, Sg, T);

  if nargout==2, [Zz, Zp] = zp2tf(Zz, Zp, Zg); end

