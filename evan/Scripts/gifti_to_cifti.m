function gifti_to_cifti(giftiL,giftiR,outname)
%function gifti_to_cifti(giftiL,giftiR,outnamestem)

templatefile = '/home/data/evan/Scripts/cifti_template.dtseries.nii';
ciftistruct = ft_read_cifti_mod(templatefile);
numverts = size(ciftistruct.data,1);

% medialmaskL = '/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.L.32k_fs_LR.func.gii';
% medialmaskdataL = gifti(medialmaskL);
% nonmedialindL = logical(medialmaskdataL.cdata==0);
% 
% 
% 
% medialmaskR = '/data/cn4/laumannt/32k_ConteAtlas_v2/medial_wall.R.32k_fs_LR.func.gii';
% medialmaskdataR = gifti(medialmaskR);
% nonmedialindR = logical(medialmaskdataR.cdata==0);

dataL = gifti(giftiL);
dataL = dataL.cdata;
surfLverts = size(dataL,1);
nonmedialindL = logical(ciftistruct.brainstructure(1:surfLverts) == 1);
dataL = dataL(nonmedialindL,:);

dataR = gifti(giftiR);
dataR = dataR.cdata;
surfRverts = size(dataR,1);
nonmedialindR = logical(ciftistruct.brainstructure((surfLverts+1) : (surfLverts+surfRverts)) == 2);
dataR = dataR(nonmedialindR,:);

combineddata = zeros(numverts,max(size(dataL,2),size(dataR,2)));

if size(dataL,2) > size(dataR,2)
    dataR = padarray(dataR,[0 (size(dataL,2) - size(dataR,2))]);
elseif size(dataL,2) < size(dataR,2)
    dataL = padarray(dataL,[0 (size(dataR,2) - size(dataL,2))]);
end

combineddata(1:(size(dataL,1) + size(dataR,1)),:) = [dataL; dataR];
ciftistruct.data = combineddata;

ft_write_cifti_mod(outname,ciftistruct)
%cifti_write_wHDR(combineddata,templatefile,outname);