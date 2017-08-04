function [Q] = M_calc_modularity(groups,rmat)

% jdp 10/10/10
%
% This script calculates the modularity (Newman, 2004) of a matrix given a set of assignments
%
% USAGE: [Q] = M_calc_modularity(groups,rmat)


% find the unique groups
ugroups=unique(groups);

% initialize a matrix of groups x groups
g=size(ugroups,1);
e=zeros(g);

% calculate the strengths of edges between and within groups
for x=1:g
    for y=1:g
        e(x,y)=sum(sum(rmat(find(groups==ugroups(x)),find(groups==ugroups(y)))));
    end
end
e=e/sum(sum(e));
tr=trace(e);
e2=sum((sum(e)).^2);
Q=tr-e2;

