function cifti_to_gifti(ciftiname,outputnamestem)
%function cifti_to_gifti(ciftiname,[outputnamestem])

if ~exist('outputnamestem')
    outputnamestem = ciftiname(1:end-13);
end

data = cifti_read(ciftiname);

hems = {'L','R'};

datavertices = 0;

for hem = 1:length(hems)
    medialmask = ['/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.' hems{hem} '.32k_fs_LR.func.gii'];
    medialmaskdata = gifti(medialmask);
    nonmedialind = find(medialmaskdata.cdata==0);
    
    dataout = zeros(size(medialmaskdata.cdata,1),size(data,2));
    
    datavertices(hem+1) = length(nonmedialind) + datavertices(hem);
    
    dataindices = [(datavertices(hem)+1) : datavertices(hem+1)];
    
    dataout(nonmedialind,:) = data(dataindices,:);
    
    save(gifti(single(dataout)),[outputnamestem '_' hems{hem} '.func.gii']);
    
end