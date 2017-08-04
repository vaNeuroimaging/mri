



i=1;
k=1;
tmask=zeros(1,numel(QC(i).FD.AFNI{1,1}));
tmask(1:4)=1;
% tmask=tmask+QC(i).FD.AFNI{1,1}'>.2;
gmlimz=[-2 2];
motlimz=[-.5 .5];
vcol=6;
ord=3;

% subplot(2,1,1);
% tc1=removepoly(QC(i).EPI.GM{k,1},1,tmask);
% imagesc(tc1,gmlimz);
% colormap(gray);
% 
% subplot(2,1,2);
% hold on;
% plot(mean(tc1,1),'r');
% plot(QC(i).EPI.MEAN{k,1}-mean(QC(i).EPI.MEAN{k,1}),'k');
% hold off;
% ylim([-5 5]);

subplot(vcol,1,1);
tc1=removepoly(QC(i).EPI.Power264{k,1},1,tmask);
imagesc(tc1,gmlimz);
colormap(gray); 
colorbar;

tc3=removepoly(QC(i).EPI.Power264{k,1},ord,tmask);
subplot(vcol,1,2);
imagesc(tc3,gmlimz); 
colormap(gray); 
colorbar;

mot1=removepoly(QC(i).MOT.AFNI{k,1}',1,tmask);
subplot(vcol,1,3);
imagesc(mot1,motlimz); 
colormap(gray); 
colorbar;

mot3=removepoly(QC(i).MOT.AFNI{k,1}',ord,tmask);
subplot(vcol,1,4);
imagesc(mot3,motlimz); 
colormap(gray); 
colorbar;


tmask=QC(i).FD.AFNI{1,1}>.2;

b1=mot1(:,tmask)'\tc1(:,tmask)';
p1=mot1'*b1;
res1=tc1-p1';
subplot(vcol,1,5);
imagesc(res1,gmlimz); 
colormap(gray); 
colorbar;

b3=mot3(:,tmask)'\tc3(:,tmask)';
p3=mot3'*b3;
res3=tc3-p3';
subplot(vcol,1,6);
imagesc(res3,gmlimz); 
colormap(gray); 
colorbar;


% b1=mot1'\tc1';
% p1=mot1'*b1;
% res1=tc1-p1';
% subplot(vcol,1,5);
% imagesc(res1,gmlimz); 
% colormap(gray); 
% colorbar;
% 
% b3=mot3'\tc3';
% p3=mot3'*b3;
% res3=tc3-p3';
% subplot(vcol,1,6);
% imagesc(res3,gmlimz); 
% colormap(gray); 
% colorbar;
% 
% r4=mean(tc3,1);
% b3=r4'\tc3';
% p3=r4'*b3;
% res3=tc3-p3';
% imagesc(res3,gmlimz); 
% colormap(gray); 
% colorbar;
