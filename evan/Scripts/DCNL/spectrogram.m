function X = spectrogram(x,nfft,fs,window,noverlap,doplot,dbdown);

%SPECTROGRAM Calculate spectrogram from signal.
%   B = SPECTROGRAM(A,NFFT,Fs,WINDOW,NOVERLAP) calculates the 
%   spectrogram for the signal in vector A.  
% 
%   NFFT is the FFT size used for each frame of A.  It should be a
%   power of 2 for fastest computation of the spectrogram.
%
%   Fs is the sampling frequency. Since all processing parameters are
%   in units of samples, Fs does not effect the spectrogram itself,
%   but it is used for axis scaling in the plot produced when
%   SPECTROGRAM is called with no output argument (see below).
%
%   WINDOW is the length M window function applied, IN ZERO-PHASE
%   FORM, to each frame of A.  M cannot exceed NFFT.  For M<NFFT,
%   NFFT-M zeros are inserted in the FFT buffer (for interpolated
%   zero-phase processing).  The window should be supplied in CAUSAL
%   FORM.
%
%   NOVERLAP is the number of samples the sections of A overlap, if
%   nonnegative.  If negative, -NOVERLAP is the "hop size", \ie, the
%   number of samples to advance successive windows.  (The overlap is
%   the window length minus the hop size.)  The hop size is called
%   NHOP below.  NOVERLAP must be less than M.
%
%   Thus, SPECTROGRAM splits the signal into overlapping segments of
%   length M, windows each segment with the length M WINDOW vector, in
%   zero-phase form, and forms the columns of B with their
%   zero-padded, length NFFT discrete Fourier transforms.
%
%   With no output argument B, SPECTROGRAM plots the dB magnitude of
%   the spectrogram in the current figure, using
%   IMAGESC(T,F,20*log10(ABS(B))), AXIS XY, COLORMAP(JET) so the low
%   frequency content of the first portion of the signal is displayed
%   in the lower left corner of the axes.
%
%   Each column of B contains an estimate of the short-term,
%   time-localized frequency content of the signal A.  Time increases
%   linearly across the columns of B, from left to right.  Frequency
%   increases linearly down the rows, starting at 0.
%
%   If A is a length NX complex signal, B is returned as a complex
%   matrix with NFFT rows and
%        k = floor((NX-NOVERLAP)/(length(WINDOW)-NOVERLAP)) 
%          = floor((NX-NOVERLAP)/NHOP)
%   columns.  When A is real, only the NFFT/2+1 rows are needed when
%   NFFT even, and the first (NFFT+1)/2 rows are sufficient for
%   inversion when NFFT is odd.
%
%   If DOPLOT is nonzero, the spectrogram will be plotted in any case.
%   DBDOWN (default=100) sets the clip level below the maximum in dB.
%
%   See also: Matlab's SPECGRAM function.
%   02/04/02/jos: Created
%   02/12/04/jos: Added dbdown

if nargin<7, dbdown=100; end
if nargin<6, doplot=0; end
if nargin<5, noverlap=256; end
if nargin<4, window=hamming(512); end
if nargin<3, fs=1; end
if nargin<2, nfft=2048; end

x=x(:);

nfft
M = length(window)
if (M<2) error(...
  'spectrogram: Expect complete window, not just its length'); 
end;
Modd = mod(M,2); % 0 if M even, 1 if odd
Mo2 = (M-Modd)/2;
w = window(:); % Make sure its a column
zp = zeros(nfft-M,1);
wzp = [w(Mo2+1:M);zp;w(1:Mo2)];

if noverlap<0
  nhop = - noverlap;
  noverlap = M-nhop;
else
  nhop = M-noverlap;
end

nx = length(x);
nframes = 1+floor((nx-noverlap)/nhop);

X = zeros(nfft,nframes);
x = x(:); % make sure it's a column
xoff = 0;
for m=1:nframes-1
  xframe = x(xoff+1:xoff+M); % extract frame of input data
  xoff = xoff + nhop;   % advance in-pointer by hop size
  xzp = [xframe(Mo2+1:M);zp;xframe(1:Mo2)];
  xw = wzp .* xzp;
  X(:,m) = fft(xw);
end

if (nargout==0) | doplot
  t = (0:nframes-1)*nhop/fs;
  f = 0.001*(0:nfft-1)*fs/nfft;
  Xdb = 20*log10(abs(X));
  Xmax = max(max(Xdb));
  % Clip lower limit to -dbdown dB so nulls don't dominate:
  clipvals = [Xmax-dbdown,Xmax];
  imagesc(t,f,Xdb,clipvals);
  % grid;
  axis xy;
  colormap(jet);
  xlabel('Time (sec)');
  ylabel('Freq (kHz)');
end
