function [allbrik] = create_264image()


% do it my way, make all ROIs equal size by centering on the voxels
ballDir = '/misc/data48/powerjd/bin/genbin/reffiles/plitt_264/balls';
for i=1:264
    ballfile=[ballDir '/ball_' num2str(i) '.nii.gz'];
    tmp=load_nii(ballfile);
    ballimg=tmp.img;
    if i==1
        allimg=ballimg;
        totimg=ballimg;
    else
        allimg=allimg+(ballimg*i);
        totimg=cat(4,totimg,ballimg); 
    end
    fprintf('%d\t%d\n',i,nnz(ballimg));
end

scriptdir=['/misc/data48/powerjd/bin/genbin/scripts'];
ofile=[scriptdir '/Power264_stack.nii.gz'];
tmp.img=allimg;
save_nii(tmp,ofile);

ofile=[scriptdir '/Power264_indiv.nii.gz'];
tmp.img=totimg;
tmp.hdr.dime.dim(1)=4;
tmp.hdr.dime.dim(5)=size(totimg,4);
save_nii(tmp,ofile);



% this is how Mark did it, centering exactly on the coordinates which have
% different numbers of voxels in some radius
% ballDir ='/misc/data16/barneska2/ABIDE/PowerEtAl_masks/3mm/';
% clear totimg allimg;
% for i=1:264
%     ballfile=[ballDir 'ball_' num2str(i) '+tlrc.HEAD'];
%     [jk ballimg hdr]=BrikLoad(ballfile);
%     if i==1
%         allimg=ballimg;
%         totimg=ballimg;
%     else
%         allimg=allimg+(ballimg*i);
%         totimg=cat(4,totimg,ballimg);
%     end
%     fprintf('%d\t%d\n',i,nnz(ballimg));
% end
% 
% cd('/misc/data48/powerjd/bin/genbin/scripts');
% Opt.Prefix='Power264_MarkROIs_stack';
% Opt.OverWrite='y';
% WriteBrik(allimg,hdr,Opt);
% system('3dAFNItoNIFTI -prefix Power264_MarkROIs_stack.nii.gz Power264_MarkROIs_stack+tlrc.BRIK');
% 
% Opt.Prefix='Power264_MarkROIs_indiv';
% Opt.OverWrite='y';
% hdr.DATASET_DIMENSIONS(4)=size(totimg,4);
% hdr.DATASET_RANK(2)=size(totimg,4);
% % hdr.BRIK_LABS=[];
% WriteBrik(totimg,hdr,Opt);
% system('3dAFNItoNIFTI -overwrite -prefix Power264_MarkROIs_indiv.nii.gz Power264_MarkROIs_indiv+tlrc.BRIK');
% 
% disp hi



