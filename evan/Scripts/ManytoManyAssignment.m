function [grp1assigns , grp2assigns] = ManytoManyAssignment(distance_metrics,maxnummatches)

if ~exist('maxnummatches')
    maxnummatches = 7;
end

[g1tog2_assign,g2tog1_assign,subtogroup_cost, grouptosub_cost] = munkres_mult(distance_metrics,maxnummatches);

unassigned_g1 = [1:length(g1tog2_assign)];
unassigned_g2 = [1:length(g2tog1_assign)];

matrix = zeros(length(g1tog2_assign),length(g2tog1_assign));
for i = 1:length(g1tog2_assign)
    if g1tog2_assign(i) > 0
        matrix(i,g1tog2_assign(i)) = 1;
    end
end
for i = 1:length(g2tog1_assign)
    if g2tog1_assign(i) > 0
        matrix(g2tog1_assign(i),i) = matrix(g2tog1_assign(i),i) + 1;
    end
end


grp1assigns = zeros(length(g1tog2_assign),1);
grp2assigns = zeros(length(g2tog1_assign),1);

unassigned_g1(logical(sum(matrix,2)==0)) = [];
unassigned_g2(logical(sum(matrix,1)==0)) = [];

[g1inds g2inds] = find(matrix==2);

for i = 1:length(g1inds)
    grp1assigns(g1inds(i)) = i;
    grp2assigns(g2inds(i)) = i;
    unassigned_g1(unassigned_g1==g1inds(i)) = [];
    unassigned_g2(unassigned_g2==g2inds(i)) = [];
end

while (numel(unassigned_g1) > 0) || (numel(unassigned_g2) > 0)

    for g1ind = unassigned_g1(:)'
        theseassigns = find(matrix(g1ind,:));
        theseassigns(grp2assigns(theseassigns) == 0) = [];
        if ~isempty(theseassigns)
            [mindist mini] = min(distance_metrics(g1ind,theseassigns));
            grp1assigns(g1ind) = grp2assigns(theseassigns(mini));
            matrix(g1ind,theseassigns(mini)) = 2;
            unassigned_g1(unassigned_g1==g1ind) = [];
        end
    end
    
    for g2ind = unassigned_g2(:)'
        theseassigns = find(matrix(:,g2ind));
        theseassigns(grp1assigns(theseassigns) == 0) = [];
        if ~isempty(theseassigns)
            [mindist mini] = min(distance_metrics(theseassigns,g2ind));
            grp2assigns(g2ind) = grp1assigns(theseassigns(mini));
            matrix(theseassigns(mini),g2ind) = 2;
            unassigned_g2(unassigned_g2==g2ind) = [];
        end
    end


end
