function [Ci community_nums] = infomap_wrapper_v157(roifile,rmat,pajekfilename,reps,deleteit)
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
mat2pajek(rmat,roifile,pajekfilename);


% this will be the relevant output of infomap
[pathstr,cluname,ext] = filenamefinder(pajekfilename,'dotsout');
clufile = [ pathstr '/' cluname '.clu' ];
treefile = [ pathstr '/' cluname '.tree' ];

% obtain seed #
clear randnum;
randnum=ceil(rand*1000000);

% find out which computer we're on (infomap is compiled  - system specific
command=['uname -m'];
[ trash systemversion]=system(command);
systemversion(end)=[]; % remove that carriage return
switch systemversion
    case 'i686'
        infomapcommand='infomap_i686';
    case 'x86_64'
        infomapcommand='/data/cn4/evan/Infomap/Infomap-0.15.7/Infomap --clu --tree';
    otherwise
        fprintf('Need to compile infomap yourself and add it to the infomap_wrapper switch possibilities.\n');
end



% run infomap
c=clock;
fprintf('\t%2.0f:%2.0f:%2.0f: infomap beginning\n',c(4),c(5),c(6));
%command = [infomapcommand ' -s' num2str(randnum) ' -N' num2str(reps) ' ' pajekfilename ' ' pathstr ' > junk.txt'];
command = [infomapcommand ' -s' num2str(randnum) ' ' pajekfilename ' ' pathstr ' > junk.txt'];
system(command);
c=clock;
fprintf('\t%2.0f:%2.0f:%2.0f: infomap finished\n',c(4),c(5),c(6));

% suppress and remove screen output
command = [ 'rm junk.txt' ];
system(command);

% get the column of assignments
fid=fopen(clufile,'r');
tsc=textscan(fid,'%d','HeaderLines',1);
Ci=double(tsc{1});

[communities walkers nodes] = textread(treefile,'%s%f%s','delimiter',' ','headerlines',1);
maxlevels = 0;

for nodenum = 1:length(nodes)
    
    nodeIDs(nodenum) = str2num(nodes{nodenum}(2:end-1));
    
    coloninds = strfind(communities{nodenum},':');
    communities{nodenum}(coloninds) = ' ';
    communities{nodenum} = str2num(communities{nodenum});
    communities{nodenum} = communities{nodenum}(1:end-1);
    if maxlevels < length(communities{nodenum})
        maxlevels = length(communities{nodenum});
    end
end

communities(nodeIDs) = communities;

community_nums = zeros(length(nodes),maxlevels);

for nodenum = 1:length(nodes)
    community_nums(nodenum,1) = communities{nodenum}(1);
    for levelnum = 2:maxlevels
        higherlevelval = community_nums(nodenum,levelnum-1) * 100;
        if length(communities{nodenum}) >= levelnum
            community_nums(nodenum,levelnum) = higherlevelval + communities{nodenum}(levelnum);
        else
            community_nums(nodenum,levelnum) = higherlevelval;
        end
    end
end

for levelnum = 2:maxlevels
    uniquevals = unique(community_nums(:,levelnum));
    temp = community_nums(:,levelnum);
    for valnum = 1:length(uniquevals)
         community_nums(temp==uniquevals(valnum),levelnum) = valnum;
    end
end

if deleteit
    %system(['rm ' pathstr '/' cluname '*' ]);
    
end

