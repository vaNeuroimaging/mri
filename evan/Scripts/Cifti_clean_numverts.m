function cleaned = Cifti_clean_numverts(data,minsize)

if isstr(data)
    data = cifti_read(data);
end

neighbors = smartload('/data/cn4/evan/fsaverage_LR32k/Cifti_surf_neighbors_LR_normalwall.mat');

data(size(neighbors,1)+1:end,:) = 0;


cleaned = zeros(size(data));

for col = 1:size(data,2)
    disp(['Column ' num2str(col) ' of ' num2str(size(data,2))])
    coldata = data(:,col);
    
    blankinds = [];
    outputdata = coldata;
    
    IDs = unique(coldata); IDs(IDs==0) = [];
    
    for ID = IDs(:)'
        
        %find which verticies meet the threshold criteria
        data_inthresh = find((coldata >= ID-.001) .* (coldata <= ID+.001));
        
        %initialize the metric keeping track of unique cluster identifiers
        clustereddata = zeros(size(coldata));
        
        for vertex = data_inthresh'
            
            %find the neighbors of this vertex
            vertexneighbors = neighbors(vertex,:);
            
            %find which of those neighbors also pass the thresholds
            vertexneighbors_inthresh = intersect(data_inthresh,vertexneighbors);
            
            %find if those neighbors have already been assigned different cluster values
            uniqueneighborvals = unique(clustereddata(vertexneighbors_inthresh));
            uniqueneighborvals(uniqueneighborvals==0) = [];
            
            %if no neighbors have cluster identifiers, assign them the number of this vertex as a unique cluster identifier
            if isempty(uniqueneighborvals)
                clustereddata(vertexneighbors_inthresh) = vertex;
                %if there is only one previous cluster identifier present, make all the neighbors that value
            elseif length(uniqueneighborvals)==1
                clustereddata(vertexneighbors_inthresh) = uniqueneighborvals;
                %if there are multiple cluster identifier values in the neighborhood, merge them into one
            else
                for valuenum = 2:length(uniqueneighborvals)
                    clustereddata(clustereddata==uniqueneighborvals(valuenum)) = uniqueneighborvals(1);
                end
            end
            
        end
        
        %find out what the unique cluster identifier values are
        uniqueclustervals = unique(clustereddata);
        uniqueclustervals(uniqueclustervals==0) = [];
        
        
        for clusternum = 1:length(uniqueclustervals)
            
            if nnz(clustereddata==uniqueclustervals(clusternum)) < minsize
                outputdata(clustereddata==uniqueclustervals(clusternum)) = 0;
                blankinds = [blankinds; find(clustereddata==uniqueclustervals(clusternum))];
            end
            
        end
    end
    
    outputdata(coldata==0) = -1;
    
    while ~isempty(blankinds)
        temp_blankinds = blankinds;
        temp = outputdata;
        for ind = blankinds(:)'
            indneighs = neighbors(ind,2:7); indneighs(isnan(indneighs)) = [];
            neighvals = outputdata(indneighs); neighvals(neighvals==0) = [];
            if ~isempty(neighvals)
                temp(ind) = mode(neighvals);
                temp_blankinds(temp_blankinds==ind) = [];
            end
        end
        blankinds = temp_blankinds;
        outputdata = temp;
    end
    outputdata(outputdata==-1) = 0;
    cleaned(:,col) = outputdata;
end



