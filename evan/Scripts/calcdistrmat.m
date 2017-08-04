function [plotarray] = calcdistrmat(mat,distbins,distmat)

numrois=size(mat,1);
numbins=length(distbins)-1;

binmask=zeros(numrois,numrois,numbins);
plotarray=zeros(1,numbins);
for i=1:numbins
    binmask(:,:,i)=(distmat >= distbins(i)) & (distmat < distbins(i+1));
    binmask=logical(binmask);
    plotarray(1,i)=mean(mat(binmask(:,:,i)));
end