function make_network_borders(filename)
%function make_network_borders(filename)

bufsize=16384;
% Read in node neighbor file generated from caret -surface-topology-neighbors
[neighbors(:,1) neighbors(:,2) neighbors(:,3) neighbors(:,4)...
    neighbors(:,5) neighbors(:,6) neighbors(:,7)] = ...
    textread('/data/cn4/laumannt/standard_mesh_atlases/Conte69_atlas.LR.32k_fs_LR_glasser/fsaverage_LR32k/node_neighbors.txt','%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);
neighbors = neighbors+1;

if strcmp(filename(end-12:end),'.dtseries.nii')
    
    cifti_to_gifti(filename,filename(1:end-13));
    
    filenames = {[filename(1:end-13) '_L.func.gii'], [filename(1:end-13) '_R.func.gii']};
    
elseif strcmp(filename(end-8:end),'.func.gii')
    
    filenames = {filename};
   
else
    
    disp('Input must be a .func.gii or a .dtseries.nii file')
    return
end

for filenum = 1:length(filenames)
    thisfilename = filenames{filenum};
    
    networks = gifti(thisfilename);
    networks = networks.cdata;
    
    output = zeros(size(networks));
    
    for i = 1:length(networks)
        nodeneighs = neighbors(i,:);
        nodeneighs(isnan(nodeneighs)) = [];
        if numel(unique(networks(nodeneighs)))>1
            output(i) = networks(i);
        end
    end
    
    save(gifti(single(output)),[thisfilename(1:end-9) '_borders.func.gii']);
    
end

if length(filenames)>1
    gifti_to_cifti([filenames{1}(1:end-9) '_borders.func.gii'],[filenames{2}(1:end-9) '_borders.func.gii'],[filename(1:end-13) '_borders']);
    delete([filename(1:end-13) '_*.func.gii'])
end
    
    
    
    
    
    