function make_hierarchical_dendrogram(recolored_networks_bythresh_file,thresholds)
%recolored_networks_bythresh_file = 'MSC02_infomap_allcolumns_recoloredv3.dscalar.nii';
%thresholds = [.003 .004 .005:.005:.05];

power_surf_colormap = [1 0 0;0 0 .8;1 1 0;1 .8 .6;0 1 0;1 .6 1;0 .6 .6;0 0 0;.35 0 .65;.2 1 1;1 .5 0;.65 .25 1;0 .25 .6;.6 1 .6;.2 .3 1;1 1 .8;0 .4 0];

apriori_merge_structure = [8 15;
    17 10;
    6 10;
    7 1;
    4 2;
    14 1;
    11 10;
    16 1;
    15 1;
    5 3;
    12 10;
    3 1;
    13 1;
    10 2;
    9 2;
    2 1];

data = ft_read_cifti_mod(recolored_networks_bythresh_file); data = data.data;
data(data > 17) = 0;
data = data(1:59412,:);
% temp = data(:,4:end);
% temp(temp==7) = 0;
% data(:,4:end) = temp;


colors = unique(data); colors(colors<1) = [];

Zplus = zeros(length(colors)-1,3);
mergecounter = 0;

for col = 2:size(data,2)
    colors_prevcol = unique(data(:,col-1)); colors_prevcol(colors_prevcol<1) = [];
    colors_thiscol = unique(data(:,col)); colors_thiscol(colors_thiscol<1) = [];
    colorslost = setdiff(colors_prevcol,colors_thiscol);
    
    for lostcolor = colorslost(:)'
        if nnz(data(:,col+1:end)==lostcolor)==0
            allvalues_inlostcolorinds = data(data(:,col-1)==lostcolor,col);
            allvalues_inlostcolorinds(allvalues_inlostcolorinds==0) = [];
            [color_mergedinto, lostcolor_vertcount] = mode(allvalues_inlostcolorinds);
            disp([num2str(lostcolor) ' merged into ' num2str(color_mergedinto) ': ' num2str(100*lostcolor_vertcount/numel(allvalues_inlostcolorinds)) '%'])
            mergecounter = mergecounter+1;
            Zplus(mergecounter,:) = [lostcolor color_mergedinto thresholds(col) ];
        end
    end
end

interval = mean(diff(thresholds));
colors_thiscol(colors_thiscol==1) = [];
for remainingcolornum = 1:length(colors_thiscol)
    apriori_inds(remainingcolornum) = find(apriori_merge_structure(:,1)==colors_thiscol(remainingcolornum));
end
addtoZplus = [apriori_merge_structure(sort(apriori_inds),:) (([1:length(colors_thiscol)]' * interval) + thresholds(end))];
for row = 1:size(addtoZplus,1)
    ind = find(Zplus(:,1)==addtoZplus(row,2));
    if ~isempty(ind)
        addtoZplus(row,2) = Zplus(ind,2);
    end
end

addtoZplus(end,3) = addtoZplus(end,3) + interval;

Zplus(mergecounter+1:end,:) = addtoZplus;


winner = Zplus(:,2);
loser = Zplus(:,1);


Ztemp_orig = Zplus(:,1:2);
Ztemp = Ztemp_orig;
for i = 1:length(colors)
    Ztemp(Ztemp_orig==colors(i)) = i;
end
Zplus(:,1:2) = Ztemp;


for color = colors(:)'
    [rows,columns] = find(Zplus(:,1:2)==color);
    uniquerows = unique(rows);
    for k = 1:length(rows)
        uniquerowind = find(uniquerows==rows(k));
        if uniquerowind > 1
            Zplus(rows(k),columns(k)) = (uniquerows(uniquerowind-1) + length(colors));
        end
    end
end

figure;
lines = dendrogram(Zplus);
Xlabels = get(gca,'XTickLabel');
Xlabels = colors(str2num(Xlabels));
tempfig = gcf;

figure;
set(gcf,'Position',[813 30 1102 805])
set(gcf,'Color',[.9 .9 .9])
set(gca,'Color',[.9 .9 .9])
for linenum = 1:length(lines)
    
   ordered_mergingcolors = Xlabels(logical((Xlabels==winner(linenum)) + (Xlabels==loser(linenum))));
    
   Xcoords =  get(lines(linenum),'XData');
   Ycoords =  get(lines(linenum),'YData');
   
   swapdirs = (Xcoords(1) > Xcoords(4));
   
   left_xcoords = [Xcoords(1) Xcoords(2) mean([Xcoords(1) Xcoords(4)])];
   left_ycoords = [Ycoords(1) Ycoords(2) Ycoords(2)];
   
   right_xcoords = [mean([Xcoords(1) Xcoords(4)]) Xcoords(3) Xcoords(4)];
   right_ycoords = [Ycoords(3) Ycoords(3) Ycoords(4)];
   
   if linenum==length(lines)
       left_xcoords = left_xcoords(1:2);
       left_ycoords = left_ycoords(1:2);
       right_xcoords = right_xcoords(2:3);
       right_ycoords = right_ycoords(2:3);
   end
   
   line(left_xcoords,left_ycoords,'Color',power_surf_colormap(ordered_mergingcolors(1+swapdirs),:),'Linewidth',4)
   line(right_xcoords,right_ycoords,'Color',power_surf_colormap(ordered_mergingcolors(2-swapdirs),:),'Linewidth',4)
   
end

set(gca,'XTickLabel',[])
%hline(thresholds(end)+(interval/2),'r--')
%ylim([0 Zplus(end,3)+(interval / 2)])
ylim([0 thresholds(end) + interval/2])
xlim([0 length(Xlabels)+.5])

close(tempfig)

dotsloc = strfind(recolored_networks_bythresh_file,'.');
outname = [recolored_networks_bythresh_file(1:dotsloc(1)-1) '_hierarchical_dendrogram.pdf'];
export_fig(gca,outname)