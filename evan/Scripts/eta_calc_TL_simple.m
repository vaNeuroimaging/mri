function eta=eta_calc_TL_simple(vec1,vec2)
%tic


% mean correlation value over all locations in both images
Mgrand  = (mean(vec1) + mean(vec2))/2;
%
% mean value matrix for each location in the 2 images
Mwithin = (vec1+vec2)/2;
SSwithin = sum((vec1-Mwithin).^2) + sum((vec2-Mwithin).^2);
SStot    = sum((vec1-Mgrand).^2) + sum((vec2-Mgrand).^2);
%
% N.B. SStot = SSwithin + SSbetween so eta can also be written as SSbetween/SStot
eta = 1 - SSwithin/SStot;
