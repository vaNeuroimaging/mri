assignmentsfile = 'rawassn.txt';

minsize = 5;

mincol = 3;

%% Size threshold

simplified = modify_clrfile('simplify',assignmentsfile,minsize);


%% Regularize

regularized = rawoutput2clr(simplified);
regularized(regularized<=1) = 0; regularized = regularized-1;
dlmwrite([assignmentsfile(1:end-4) '_minsize' num2str(minsize) '_regularized.txt'],regularized,'delimiter','\t');


%% Create initial consensus by accepting all assignments at the mincol threshold and assigning unassigned nodes to their higher threshold assignments

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
