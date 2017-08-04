function [img dims o]=nii_give_RAS(im)

% NIFTI qform and sform codes are relative to RAS
%
% does not matter which code used for showing the arrays
% but we do want to flip the arrays to be RAS or close to it
%
% also we need to scale the image, sometimes uint8 and a scaling factor are used to save space instead of single or double (float)


tmp=load_untouch_nii(im);
if tmp.hdr.dime.scl_slope~=0
    fprintf('Rescaling the image\n');
    tmp.img=tmp.img*tmp.hdr.dime.scl_slope+tmp.hdr.dime.scl_inter;
end

dims=tmp.hdr.dime.pixdim(2:4);

[jk o.sx]=system(['fslhd ' im ' | grep sform_xorient | awk ''{print $2}''' ]); o.sx(end)=[];
[jk o.sy]=system(['fslhd ' im ' | grep sform_yorient | awk ''{print $2}''' ]); o.sy(end)=[];
[jk o.sz]=system(['fslhd ' im ' | grep sform_zorient | awk ''{print $2}''' ]); o.sz(end)=[];
[jk o.scode]=system(['fslhd ' im ' | grep sform_code | awk ''{print $2}''' ]); o.scode(end)=[]; o.scode=str2num(o.scode);
[jk o.sRx]=system(['fslhd ' im ' | grep sto_xyz:1 | awk ''{print $2 " " $3 " " $4}''' ]); o.sRx(end)=[]; o.sR(1,1:3)=str2num(o.sRx);
[jk o.sRy]=system(['fslhd ' im ' | grep sto_xyz:2 | awk ''{print $2 " " $3 " " $4}''' ]); o.sRy(end)=[]; o.sR(2,1:3)=str2num(o.sRy);
[jk o.sRz]=system(['fslhd ' im ' | grep sto_xyz:3 | awk ''{print $2 " " $3 " " $4}''' ]); o.sRz(end)=[]; o.sR(3,1:3)=str2num(o.sRz);


[jk o.qx]=system(['fslhd ' im ' | grep qform_xorient | awk ''{print $2}''' ]); o.qx(end)=[];
[jk o.qy]=system(['fslhd ' im ' | grep qform_yorient | awk ''{print $2}''' ]); o.qy(end)=[];
[jk o.qz]=system(['fslhd ' im ' | grep qform_zorient | awk ''{print $2}''' ]); o.qz(end)=[];
[jk o.qcode]=system(['fslhd ' im ' | grep qform_code | awk ''{print $2}''' ]); o.qcode(end)=[]; o.qcode=str2num(o.qcode);
[jk o.qRx]=system(['fslhd ' im ' | grep qto_xyz:1 | awk ''{print $2 " " $3 " " $4}''' ]); o.qRx(end)=[]; o.qR(1,1:3)=str2num(o.qRx);
[jk o.qRy]=system(['fslhd ' im ' | grep qto_xyz:2 | awk ''{print $2 " " $3 " " $4}''' ]); o.qRy(end)=[]; o.qR(2,1:3)=str2num(o.qRy);
[jk o.qRz]=system(['fslhd ' im ' | grep qto_xyz:3 | awk ''{print $2 " " $3 " " $4}''' ]); o.qRz(end)=[]; o.qR(3,1:3)=str2num(o.qRz);

if o.scode>0
     [img perm]=perm_and_flip(tmp,o.sR);
elseif o.qcode>0
     [img perm]=perm_and_flip(tmp,o.qR);
else
    fprintf('ERROR WARNING in nii_give_RAS.m: qformcode and sformcode of %s are both 0 so conversion to RAS is unknown',niifile);
end

if numel(perm)==4
    perm(4)=[];
end
dims=dims(perm);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [img perm] = perm_and_flip(tmp,R)

% define the principal directions
signR=sign(R);
absR=abs(R);
maxR=max(absR,[],2);
maxR1=max(absR,[],1);

% if (sum(ismember(maxR,maxR1))==3)
%     for i=1:3
%         topR(i,:)=absR(i,:)==maxR(i);
%     end
% else

[s ss]=sort(absR(:));
s3(ss)=[1:9];
s3=reshape(s3,[3 3]);
topR=zeros(3,3);
[a b]=find(s3==max(s3(:)));
s3(a,:)=[0 0 0];
s3(:,b)=[0 0 0]';
topR(a,b)=1;
[a b]=find(s3==max(s3(:)));
s3(a,:)=[0 0 0];
s3(:,b)=[0 0 0]';
topR(a,b)=1;
[a b]=find(s3==max(s3(:)));
s3(a,:)=[0 0 0];
s3(:,b)=[0 0 0]';
topR(a,b)=1;



flipR=topR.*signR;
for i=1:3
    [jk perm(i) flipr(i) ]=find(flipR(i,:));
end

% have to preserve the time dimensions explicitly if present
d=size(tmp.img);
if numel(d)>3
    perm(4)=4;
end


img=permute(tmp.img,perm);

% now flip where needed to get RAS
for i=1:3
    if flipr(i)<0
        img=flipdim(img,i);
    end
end



    