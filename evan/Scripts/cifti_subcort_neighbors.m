function neighbors = cifti_subcort_neighbors(ciftifile)

slashloc = strfind(ciftifile,'/');
if isempty(slashloc)
    slashloc = 0;
end

thisdir = pwd;
if ~isempty(ciftifile(1:slashloc(end)))
cd(ciftifile(1:slashloc(end)))
end

if strcmp(ciftifile(end-13:end),'.dtseries.nii')
    evalc(['!wb_command -cifti-convert -to-gifti-ext ' ciftifile ' ' ciftifile((slashloc(end)+1):(end-13)) '.func.gii']);
    ciftifile = [ciftifile(slashloc(end)+1:end-13) '.func.gii'];
    deletestuff = 1;
else
    ciftifile = [ciftifile(slashloc(end)+1:end)];
    deletestuff = 0;
end


bufsize = 1048576;
ciftiheadertext = textread(ciftifile,'%s','delimiter','\r','bufsize',bufsize);

cd(thisdir)

xyzcoords = [];
findtext = '<VoxelIndicesIJK>';
endtext = '</VoxelIndicesIJK>';
readnext = 0;


for row = 1:length(ciftiheadertext)
    rowtext = ciftiheadertext{row};
    if length(rowtext)>length(findtext) && (strcmp(rowtext(1:length(findtext)),findtext))
        xyzcoords(end+1,:) = str2num(rowtext(length(findtext)+1:end));
        readnext = 1;
    elseif strcmp(rowtext,endtext)
        readnext = 0;
    elseif readnext==1
        xyzcoords(end+1,:) = str2num(rowtext);
    end
end

neighbors = zeros(size(xyzcoords,1),27);
neighbors(:) = NaN;

for ind = 1:size(neighbors,1);
    distmat = single(pdist2(xyzcoords,xyzcoords(ind,:)));
    these_neighbors = find(distmat < 2)';
    these_neighbors(these_neighbors==ind) = [];
    neighbors(ind,1:(length(these_neighbors)+1)) = [ind these_neighbors];
end

% distmat=single(squareform(pdist(xyzcoords)));
% [neighborsi,neighborsj] = find(distmat<2);
% 
% for ind = 1:size(neighbors,1);
%     these_neighbors = [ind setdiff(neighborsj(neighborsi==ind),ind)'];
%     neighbors(ind,1:length(these_neighbors)) = these_neighbors;
% end

if deletestuff
    delete(ciftifile);
    delete([ciftifile '.data']);
end
