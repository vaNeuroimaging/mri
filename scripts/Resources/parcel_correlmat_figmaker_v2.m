function parcel_correlmat_figmaker_v2(corrmat,assignmentsfile,limits,titlename)
%parcel_correlmat_figmaker_v2(corrmat,assignmentsfile,limits,titlename)

networklabels = {'Unassigned','Default','Visual','FrontoPar','PrimaryVisual','DorsalAttn','Premotor','VentAttn','Salience','CingOperc','MortorHand','MotorMouth','Auditory','MTL1','MTL2','MedPar','ParOccip','MotorFoot','Unknown','Unknown'};
colors = [.5 .5 .5;1 0 0;0 0 .6;.9 .9 0;1 .7 .4;0 .8 0;1 .6 1;0 .6 .6;0 0 0;.3 0 .6;.2 1 1;1 .5 0;.6 .2 1;0 .2 .4;.2 1 .2;0 0 1;1 1 1;.0 .4 0;.8 .35 .5;.5 .75 .2];

nodes = size(corrmat,1);

yaxis_label_position_params = [1.4 (nodes/100) * -1]; %horiz vert
xaxis_label_position_params = [(nodes/100) (nodes/100)*1.75 + .6]; %horiz vert

corrmat(logical(diag(ones(size(corrmat,1),1)))) = 0;

if ischar(assignmentsfile)
    assignments = load(assignmentsfile);
else
    assignments = assignmentsfile;
end
assignments(assignments<0) = 0;
IDs = unique(assignments);
networklabels = networklabels(IDs+1);
colors = colors(IDs+1,:);

[communities, sorti] = sort(assignments);

sorted_corrmat = corrmat(sorti,sorti);

transitions = find(communities(1:end-1) ~= communities(2:end));
transitions_plusends = [1 transitions(:)' length(communities)];
centers = transitions_plusends(1:end-1) + ((transitions_plusends(2:end) - transitions_plusends(1:end-1))/2);

figure;
imagesc(sorted_corrmat,limits)


axis(gca,'off')


set(gca,'FontSize',20)

set(gcf,'Position',[813 30 1102 805])

colormap hot; hotmap = colormap; coolmap = [hotmap(:,3),hotmap(:,2),hotmap(:,1)];
hotmap = hotmap(1:end-3, :); coolmap = coolmap(1:end-3, :);
combined = [flipdim(coolmap,1); zeros(10,3); hotmap];
colormap(combined);
colorbar

boxwidth = round(size(sorted_corrmat,1) ./ 20);
xlim([-boxwidth size(sorted_corrmat,1)])
ylim([-boxwidth size(sorted_corrmat,1)])

boxedges = [0; transitions;  size(sorted_corrmat,1)];
for i = 1:(length(boxedges)-1)
    rectangle('Position',[-boxwidth+.1 boxedges(i) boxwidth (boxedges(i+1) - boxedges(i))],'FaceColor',colors(i,:))
    rectangle('Position',[boxedges(i) -boxwidth+.1 (boxedges(i+1) - boxedges(i)) boxwidth],'FaceColor',colors(i,:))
end



        


if exist('titlename')
    title(titlename);
end

if ~isempty(transitions)
    hline(transitions,'g')%hline(transitions+.5,'g')
    vline(transitions,'g')%vline(transitions+.5,'g')
end

set(gcf,'Color',[1 1 1])