function [Ci] = infomap_wrapper_mod_EG(pajekfilename,reps)
%[Ci] = infomap_wrapper_mod_EG(pajekfilename,reps)
%
%
% This script runs infomap on a pajekfile with some number of
% repetitions. It then  returns the community assignments found.



% this will be the relevant output of infomap
[pathstr,cluname,ext] = filenamefinder(pajekfilename,'dotsout');
clufile = [ pathstr '/' cluname '.clu' ];

% obtain seed #
clear randnum;
randnum=ceil(rand*1000000);

% find out which computer we're on (infomap is compiled  - system specific
command= 'uname -m';
systemversion=evalc(['!' command]);
while length(systemversion)<1
    systemversion=evalc(['!' command]);
end
systemversion(end)=[]; % remove that carriage return
switch systemversion
    case 'i686'
        infomapcommand='infomap_i686';
    case 'x86_64'
        infomapcommand='infomap_x86_64';
    otherwise
        fprintf('Need to compile infomap yourself and add it to the infomap_wrapper switch possibilities.\n');
end

% run infomap
c=clock;
fprintf('\t%2.0f:%2.0f:%2.0f: infomap beginning\n',c(4),c(5),c(6));
command = ['!/data/cn4/evan/Infomap/Infomap-0.15.7/Infomap --clu -2 -s' num2str(randnum) ' -N' num2str(reps) ' ' pajekfilename ' ' pathstr];
evalc(command);
c=clock;
fprintf('\t%2.0f:%2.0f:%2.0f: infomap finished\n',c(4),c(5),c(6));



% So parfor doesn't crap out
isclufile = exist(clufile);
while isclufile == 0
    pause(60)
    isclufile = exist(clufile);
end

Ci = textread(clufile,'%d','headerlines',1);