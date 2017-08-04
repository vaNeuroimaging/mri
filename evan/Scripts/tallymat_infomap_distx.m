function rawclrs = tallymat_infomap_distx(rawassn,distancemat,xdist,roifile)
%rawclrs = tallymat_infomap(rawassn,roifile)

if ~isnumeric(rawassn)
    mat = load(rawassn);
else
    mat = rawassn;
end




% %try to skip redundant thresholds with 99% overlap of the prev threshold
% col = 2;
% colcounter = 2;
% while col <= size(mat,2)
%     
%     col1mat = zeros(size(mat,1));
%     col2mat = zeros(size(mat,1));
%     for node = 1:size(mat,1)
%         if mat(node,col-1) > 0
%             col1mat(:,node) = single(mat(:,col-1)==mat(node,col-1));
%         end
%         if mat(node,col) > 0
%             col2mat(:,node) = single(mat(:,col)==mat(node,col));
%         end
%     end
%     
%     nodestouse = logical(col1mat + col2mat);
%     if (nnz(col1mat(nodestouse)==col2mat(nodestouse)) / nnz(nodestouse)) == 1; %>.99;
%         disp(['thresh ' num2str(colcounter) 'skipped'])
%         mat(:,col) = [];
%     else
%         col = col+1;
%     end
%     colcounter = colcounter+1;
% end
        


%Build tally matrix
tallymat = zeros(size(mat,1));
for col = 1:size(mat,2)
    for node = 1:size(mat,1)
        if mat(node,col) > 0
            tallymat(:,node) = tallymat(:,node) + single(mat(:,col)==mat(node,col));
        end
    end
end

tallymat = tallymat .* (distancemat > xdist);

clear distancemat


% %threshold tallymatrix
% thresh = .5;
% for node = 1:length(tallymat)
%     numconnected = nnz(mat(node,:)>0);
%      thisthresh = ceil(numconnected .* thresh);
%      belowthreshinds = logical(tallymat(node,:) < thisthresh);
%      tallymat(node,belowthreshinds) = 0; tallymat(belowthreshinds,node) = 0;
% 
%     %tallymat(node,:) = tallymat(node,:) ./ numconnected; %tallymat(:,node) = tallymat(:,node) ./ numconnected;
%         
% end

%tallymat(tallymat < thresh) = 0;
    

%run infomap
pajekfile = 'temp.net';
rawclrs = infomap_wrapper_faster(roifile,tallymat,pajekfile,1,1);
dlmwrite(['tallymat_assn.txt'],rawclrs,'\t');


%make figs of the tally matrix and resulting networks (if matrix is small
%enough)
if size(tallymat,1) < 2000
    
[ign sortorder] = sort(rawclrs,'ascend');

if length(sortorder)>1000
    sortorder = sortorder(sort(randperm([1:length(sortorder)],1000)));
end


figure

power_surf_colormap = [1 0 0;0 0 .6;1 1 0;1 .7 .4;0 .8 0;1 .6 1;0 .6 .6;0 0 0;.3 0 .6;.2 1 1;1 .5 0;.6 .2 1;0 .2 .4;.2 1 .2;0 0 1;1 1 .8;0 .4 0;.25 .25 .25];
if max(unique(rawclrs)) < (size(power_surf_colormap,1)-2)
    power_surf_colormap_touse = power_surf_colormap(1:max(unique(rawclrs))+2,:);
else
    power_surf_colormap_touse = [power_surf_colormap ; repmat(power_surf_colormap(end,:),(max(unique(rawclrs)) - (size(power_surf_colormap,1)-2)),1)];
end

h = subplot(1,2,1); imagesc(rawclrs(sortorder)); colormap(h,power_surf_colormap_touse)
freezeColors

h = subplot(1,2,2); imagesc(tallymat(sortorder,sortorder)); colormap(h,'default')
colorbar
end
pause(.1)