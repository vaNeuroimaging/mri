function compare_cifti_borders(borderfile1,borderfile2,distancethreshes,distancemat)
% compare_cifti_borders(borderfile1,borderfile2,[distancethreshes],[distancemat])


if ~exist('distancethreshes') || isempty(distancethreshes)
    distancethreshes = [6 8 10 12 14 16 18 20];
end


[file1path,temp,~] = fileparts(borderfile1);
[~,file1base,~] = fileparts(temp);

[file2path,temp,~] = fileparts(borderfile2);
[~,file2base,~] = fileparts(temp);


datastruct = ft_read_cifti_mod(borderfile1);
borders1 = datastruct.data;

temp = ft_read_cifti_mod(borderfile2);
borders2 = temp.data;

ncortverts = nnz(datastruct.brainstructure==1) + nnz(datastruct.brainstructure==2);
borders1(ncortverts+1:end) = 0;
borders2(ncortverts+1:end) = 0;



datastruct.dimord = 'scalar_pos';
columnnames{1} = 'Minimum distance, in mm';
columnnames{2} = 'Largest minimum distance by segment, in mm';
for i = 1:length(distancethreshes)
    columnnames{i+2} = ['Segments with at least one point farther than ' num2str(distancethreshes(i)) 'mm'];
end
datastruct.mapname = columnnames;


disp('Calculating distances')
if ~exist('distancemat','var') || isempty(distancemat)
    distances = smartload('/data/cn/data1/scripts/CIFTI_RELATED/Resources/Conte69_atlas-v2.LR.32k_fs_LR.wb/Cifti_geo_distances.mat');
else
    distances = smartload(distancemat);
end

comparisonout1 = zeros(size(borders1,1),length(distancethreshes)+2);
comparisonout1(:,1:2) = -1;
comparisonout2 = zeros(size(borders2,1),length(distancethreshes)+2);
comparisonout2(:,1:2) = -1;


mindists_1_to_2 = min(distances(find(borders1),find(borders2)),[],2);
comparisonout1(logical(borders1),1) = mindists_1_to_2;

mindists_2_to_1 = min(distances(find(borders1),find(borders2)),[],1);
comparisonout2(logical(borders2),1) = mindists_2_to_1;


disp('Identifying discrete segments for first border set')

neighbors = cifti_neighbors(borderfile1);

segments1 = identify_segments_cifti(borderfile1,neighbors,distances);
segIDs = unique(segments1); segIDs(segIDs==0) = [];

for ID = segIDs(:)'
    seginds = (segments1==ID);
    segdists = comparisonout1(seginds,1);
    segmax = max(segdists);
    comparisonout1(seginds,2) = segmax;
    
    comparisonout1(seginds,3:end) = 1;    
    comparisonout1(seginds,3:end) = comparisonout1(seginds,3:end) + repmat((segmax > distancethreshes),nnz(seginds),1);
    
end

datastruct.data = comparisonout1;
if isempty(file1path)
    outname = [file1base '_distanceto_' file2base];
else
    outname = [file1path '/' file1base '_distanceto_' file2base];
end
ft_write_cifti_mod(outname,datastruct)


disp('Identifying discrete segments for second border set')

segments2 = identify_segments_cifti(borderfile2,neighbors,distances);
segIDs = unique(segments2); segIDs(segIDs==0) = [];

for ID = segIDs(:)'
    seginds = (segments2==ID);
    segdists = comparisonout2(seginds,1);
    segmax = max(segdists);
    comparisonout2(seginds,2) = segmax;
    
    comparisonout2(seginds,3:end) = 1;    
    comparisonout2(seginds,3:end) = comparisonout2(seginds,3:end) + repmat((segmax > distancethreshes),nnz(seginds),1);
end


datastruct.data = comparisonout2;
if isempty(file2path)
    outname = [file2base '_distanceto_' file1base];
else
    outname = [file2path '/' file2base '_distanceto_' file1base];
end
ft_write_cifti_mod(outname,datastruct)

