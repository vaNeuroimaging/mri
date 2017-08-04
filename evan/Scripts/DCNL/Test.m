for x = 1:size(seeddata,1)
    disp(num2str(x))
    for y = 1:size(seeddata,2)
        for z= 1:size(seeddata,3)
        
        
        distances = gsingle([gindexmatrix(:,1)-x gindexmatrix(:,2)-y gindexmatrix(:,3)-z]);
        euclidean = sqrt(sum((distances.^2),2));
        pointswithinsphere = double(find(euclidean<=spherevoxelradius));
        %coordswithinsphere = indexmatrix(pointswithinsphere,:);
        
        pointsbetween2D = pointsbetweenmatrix(pointswithinsphere,:);
        
        if ~isempty(find(pointsbetween2D,1))
        
        allpointsbetween = reshape(pointsbetween2D,(size(pointsbetween2D,1)*size(pointsbetween2D,2)),1);
        allpointsbetween = unique(allpointsbetween);
        allpointsbetween = allpointsbetween(find(allpointsbetween));
        
        
        %             allpointsbetween = zeros(10000,3);
        %             pointsbetweencounter = 0;
        %
        %             for i = 1:size(coordswithinsphere,1)
        %                 if ~isempty(pointsbetweenimage{coordswithinsphere(i,1),coordswithinsphere(i,2),coordswithinsphere(i,3)})
        %                     pointstoadd = pointsbetweenimage{coordswithinsphere(i,1),coordswithinsphere(i,2),coordswithinsphere(i,3)};
        %                     allpointsbetween(pointsbetweencounter+1:pointsbetweencounter+size(pointstoadd,1),:) = pointstoadd;
        %                     pointsbetweencounter = length(find(allpointsbetween(:,1)));
        %                     %allpointsbetween = [allpointsbetween; [pointsbetweenimage{coordswithinsphere(i,1),coordswithinsphere(i,2),coordswithinsphere(i,3)}]];
        %                 end
        %             end
        %             allpointsbetween = allpointsbetween(1:pointsbetweencounter,:);
        
        
        
        FAvals = FAimage(allpointsbetween);
        
        %             FAvals = zeros(length(allpointsbetween),1);
        %             for point = 1:size(allpointsbetween,1)
        %                 FAvals(point) = FAimage(allpointsbetween(point));
        %             end
        
        FAoutputimagesphere(x,y,z) = mean(FAvals);
        end
        end
    end
    
end