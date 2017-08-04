function [Ci] = infomap_wrapper_faster(roifile,rmat,pajekfilename,reps,deleteit)
%
% Name:infomap_wrapper.m
% $Revision: 1.3 $
% $Date: 2011/03/08 20:30:00 $
%
% jdp 10/10/10
% 
% This script takes an roifile and a matrix, a name to write a pajekfile
% to, and then runs infomap on the pajekfile with some number of
% repetitions. It then does or does not delete the output, depending on
% deleteit status, and it returns the community assignments found.
% 
% USAGE: [modules] = infomap_wrapper(roifile,rmat,pajekfilename,reps,deleteit)


% create the pajekfile
%mat2list_TL(rmat,pajekfilename);
mat2pajek_TL(rmat,roifile,pajekfilename);

% this will be the relevant output of infomap
[pathstr,cluname,ext] = filenamefinder(pajekfilename,'dotsout');
clufile = [ pathstr '/' cluname '.clu' ];

% obtain seed #
clear randnum;
randnum=ceil(rand*1000000);

% find out which computer we're on (infomap is compiled  - system specific
% command= 'uname -m';
% systemversion=evalc(['!' command]);
% while length(systemversion)<1
%     systemversion=evalc(['!' command]);
% end
% systemversion(end)=[]; % remove that carriage return
% switch systemversion
%     case 'i686'
%         infomapcommand='infomap_i686';
%     case 'x86_64'
%         infomapcommand='infomap_x86_64';
%     otherwise
%         fprintf('Need to compile infomap yourself and add it to the infomap_wrapper switch possibilities.\n');
% end

% run infomap
c=clock;
fprintf('\t%2.0f:%2.0f:%2.0f: infomap beginning\n',c(4),c(5),c(6));
command = ['!/data/cn4/evan/Infomap/Infomap-0.15.7/Infomap --clu -2 -s' num2str(randnum) ' -N' num2str(reps) ' ' pajekfilename ' ' pathstr];
%command = [ infomapcommand ' ' num2str(randnum) ' ' pajekfilename ' ' num2str(reps) ' > junk.txt' ];
evalc(command);
c=clock;
fprintf('\t%2.0f:%2.0f:%2.0f: infomap finished\n',c(4),c(5),c(6));

% suppress and remove screen output
%command = [ 'rm junk.txt' ];
%system(command);

% get the column of assignments
fid=fopen(clufile,'r');
tsc=textscan(fid,'%d','HeaderLines',1);
Ci=double(tsc{1});

if deleteit
    command = [ 'rm ' pathstr '/' cluname '*' ];
    system(command);
end
fclose('all');

