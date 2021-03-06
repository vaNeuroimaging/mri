function [Ci] = run_infomap_on_pajekfile(pajekfilename,reps)
%[Ci] = run_infomap_on_pajekfile(pajekfilename,reps)
%
%
% This script runs infomap on a pajekfile with some number of
% repetitions. It then returns the community assignments found.



% this will be the relevant output of infomap
[pathstr,cluname,ext] = filenamefinder(pajekfilename,'dotsout');
clufile = [ pathstr '/' cluname '.clu' ];

% obtain seed #
clear randnum;
randnum=ceil(rand*1000000);


[~,~] = system(['/home/data/scripts/Infomap/Infomap-0.15.7/Infomap --clu -2 -s' num2str(randnum) ' -N' num2str(reps) ' ' pajekfilename ' ' pathstr]);


% So parfor doesn't crap out
isclufile = exist(clufile);
while isclufile == 0
    pause(60)
    isclufile = exist(clufile);
end

Ci = textread(clufile,'%d','headerlines',1);
end
