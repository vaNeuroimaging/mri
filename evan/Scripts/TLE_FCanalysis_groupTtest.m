grp1_datalist = ['/data/cn4/scratch/tunde/LUIGI/FCPROCESS_V4/CTL/NEW_TMASKLIST_mod.txt'];

grp2_datalist = ['/data/cn4/scratch/tunde/LUIGI/FCPROCESS_V4/TLE/NEW_TMASKLIST_mod.txt'];

outfile = '/data/cn4/evan/TLE/CTLvTLE_Ttest';

datafolder = '/data/cn4/evan/TLE/';

etype = 'littleendian';


[g1subjects tmasks] = textread(grp1_datalist,'%s%s');
[g2subjects tmasks] = textread(grp2_datalist,'%s%s');

hems = {'l','r'};
%%
% for hemnum = 1:length(hems)
%     hem = hems{hemnum};
% 
% for s = 1:length(g1subjects)
%     
%     dataname = [datafolder g1subjects{s} '_' hem 'Hip.4dfp.img'];
%     g1data(:,s) = read_4dfpimg_HCP(dataname);
%     
% end
% 
% for s = 1:length(g2subjects)
%     
%     dataname = [datafolder g2subjects{s} '_' hem 'Hip.4dfp.img'];
%     g2data(:,s) = read_4dfpimg_HCP(dataname);
%     
% end
% 
% 
% [H,P,CI,STATS] = ttest2(g1data',g2data');
% 
% Pimg = P';
% Timg = STATS.tstat';
%     
% write_4dfpimg(Pimg,[outfile '_' hem 'Hip_P.4dfp.img'],etype);
% write_4dfpifh_333_MNI([outfile '_' hem 'Hip_P.4dfp.ifh'],1,etype);
% 
% write_4dfpimg(Timg,[outfile '_' hem 'Hip_T.4dfp.img'],etype);
% write_4dfpifh_333_MNI([outfile '_' hem 'Hip_T.4dfp.ifh'],1,etype);
% 
% end
%%
% for s = 1:length(g1subjects)
%     load([datafolder g1subjects{s} '_264.mat'])
%     g1mat(:,:,s) = matrix;
% end
% for s = 1:length(g2subjects)
%     load([datafolder g2subjects{s} '_264.mat'])
%     matrix(logical(diag(length(matrix),0))) = 0;
%     g2mat(:,:,s) = matrix;
% end
% 
% g1mat = reshape(g1mat,[size(g1mat,1) * size(g1mat,2) , size(g1mat,3)]);
% g2mat = reshape(g2mat,[size(g2mat,1) * size(g2mat,2) , size(g2mat,3)]);
% 
% [H,P,CI,STATS] = ttest2(g1mat',g2mat');
% 
% Pmat = reshape(P',[size(matrix,1) size(matrix,2)]);
% Pmat(logical(diag(ones(length(Pmat),1),0))) = 1;
% Tmat = reshape(STATS.tstat',[size(matrix,1) size(matrix,2)]);
% Tmat(logical(diag(ones(length(Tmat),1),0))) = 0;
% 
% 
% 
%  
% matrix = Pmat;
% save(['CTLvTLE_Ttest_P_264.mat'],'matrix');
% 
% matrix = Tmat;
% save(['CTLvTLE_Ttest_T_264.mat'],'matrix');

%% Resort nodes with new ordering
power_coords = load('/data/hcp-zfs/home/laumannt/Poldrome/shared_for_washu/BB264_coords.txt');
laum_coords = load('/data/cn4/laumannt/longRestingState/group_modules/264_consensus_hemsort_7112B.txt');

load(['CTLvTLE_Ttest_P_264.mat']); Pmat = matrix;
load(['CTLvTLE_Ttest_T_264.mat']); Tmat = matrix;

for i = 1:264
    for t = 1:264
        if isequal(laum_coords(i,:),power_coords(t,:))
            ind(i) = t;
        end
    end
end

Pmat = Pmat(ind,ind);
Tmat = Tmat(ind,ind);

figure_corrmat_powernetwork(Tmat,[-7 7],[-5.8 5.8])
