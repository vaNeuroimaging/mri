assignmentsfile = 'rawassn.txt';

minsize = 5;

mincol = 35;

recolor = [1 15;2 10;3 5;4 7;5 -1;6 1;7 -1;8 2;9 9;10 -1;11 -1;12 8;13 11;22 3;70 16;27 2.5;56 11.5;42 -1];
%recolor = [];

% %% Size threshold
% 
% simplified = modify_clrfile('simplify',assignmentsfile,minsize);
% 
% 
% %% Regularize
% 
% regularized = rawoutput2clr(simplified);
% regularized(regularized<=1) = 0; regularized = regularized-1;
% dlmwrite([assignmentsfile(1:end-4) '_minsize' num2str(minsize) '_regularized.txt'],regularized,'delimiter','\t');


%% Create initial consensus by accepting all assignments at the mincol threshold and assigning unassigned nodes to their higher threshold assignments
regularized = load([assignmentsfile(1:end-4) '_minsize' num2str(minsize) '_regularized.txt']);

consensusmap = regularized(:,mincol);

unassigned = find(consensusmap<1);
for unassignedindex = unassigned'
    thisassignments = regularized(unassignedindex,mincol:end);
    thisassignments(thisassignments<1) = [];
    if ~isempty(thisassignments)
        consensusmap(unassignedindex) = thisassignments(1);
    end
end


%% Reapply size threshold

consensusvals = unique(consensusmap);
for val = consensusvals(:)'
    if nnz(consensusmap==val) < minsize
        toosmallinds = find(consensusmap==val);
        for smallindex = toosmallinds(:)'
            thisassignments = regularized(unassignedindex,mincol:end);
            thisassignments(thisassignments<1) = [];
            thisassignments(thisassignments==val) = [];
            if ~isempty(thisassignments)
                consensusmap(smallindex) = thisassignments(1);
            end
        end
        
    end
end



dlmwrite([assignmentsfile(1:end-4) '_minsize' num2str(minsize) '_regularized_consensus.txt'],consensusmap)

%% Recolor

if ~isempty(recolor)
    recolored = zeros(size(consensusmap));
    for i =1:size(recolor,1)
        recolored(consensusmap==recolor(i,1)) = recolor(i,2);
    end
    dlmwrite([assignmentsfile(1:end-4) '_minsize' num2str(minsize) '_regularized_consensus_recolored.txt'],recolored)
end
