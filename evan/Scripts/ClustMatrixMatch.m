function [MBout,Cmap2,MBoutKey]=ClustMatrixMatch(MA,MB,cmap)

% This function matches the cluster assignments between 2 already sorted
% matrices. OrgClustMat takes too long for very large arrays (Modified
% Voxel-Wise, etc), so this must do better. MA is treated as the master
% set.
% Algorithm: (1) define each networks fingerprint as its membership across 
%           all edge densities chosen. For this to work, the matrices 
%           being compared must represent IMappings at equivalent kdens. 
%           (2) Then test to see if remaining Nets look anything like the
%           master set only at their most representative kden (i.e., the
%           kden with the most ROIs set to that network).
%           (3) If that also fails, the remainder are given new network
%           values and colored either as an alternate color to a known
%           network (if Dice is above threshold but network match was not
%           exclusively strongest) or set to a random grayscale value. This
%           distinction is also made in the MBoutKey 3rd column: if val is
%           non-integer, it represents the Dice (neg if not whole kden
%           fingerprint, if integer, it represents closest network, if
%           zero, this is an unknown network).
% The largest value for each matrix should represent unclassified ROIs and
% is therefore ignored.


%% Parameters
valsA=unique(MA(:));
NnetsA=length(valsA)-1;
valsB=unique(MB(:));
NnetsB=length(valsB)-1;
[Nroi,Nkden]=size(MA);
MBout=zeros(size(MB),'single');
metric=zeros(NnetsA,NnetsB,'single');
metric2=zeros(NnetsA,NnetsB,'single');
MBoutKey=zeros(NnetsB+1,3,'single');
MBoutKey(end,:)=[max(valsA(:)),max(valsA(:)),0];
TryAgain=[];
masterA=zeros(Nroi,NnetsA,'single');
NmaxKDenA=zeros(NnetsA,1);
testB=zeros(Nroi,NnetsB,'single');
NmaxKDenB=zeros(NnetsB,1);
Dth=0.1;
NetAmax=valsA(NnetsA);
n=0;
GrandMax=size(cmap,1);
Cmap2=cmap;

%% (0) Initial conditions
% Double check master set
if (NnetsA~=valsA(NnetsA)) % Ignore last value 
for j=1:NnetsA
    MA(MA==valsA(j))=j;
end
MA(MA==valsA(end))=GrandMax;
valsA=unique(MA(:));
MBoutKey(end,:)=[max(valsA(:)),max(valsA(:)),0];
NetAmax=valsA(NnetsA);
end
% Check test set
if (NnetsB~=valsB(NnetsB)) % Ignore last value 
for j=1:NnetsB
    MB(MB==valsB(j))=j;
end
MB(MB==valsB(end))=GrandMax;
valsB=unique(MB(:));
end


figure('Position',[115,50,1560,1000]);
subplot(2,6,[1,7]);
imagesc(sortrows(MA,Nkden));colormap(cmap);hold on;axis off
plot([Nkden-.5,Nkden-0.5],[0,Nroi],'k');
plot([Nkden+.5,Nkden+0.5],[0,Nroi],'k')
title('Master key')
freezeColors
subplot(2,6,[2,8]);
imagesc(sortrows(MB,Nkden));colormap(cmap);hold on;axis off
plot([Nkden-.5,Nkden-0.5],[0,Nroi],'k');
plot([Nkden+.5,Nkden+0.5],[0,Nroi],'k')
title('Un-Sorted Test set')
freezeColors

%% (1) Test global match of mutual maxima
% Calculate Dice Coeff across full matrix for each Network
disp('<< Testing Global Fingerprint')
for j=1:NnetsA
    foo=MA==valsA(j);
    for k=1:NnetsB
        foob=MB==valsB(k);
        metric(j,k)=2*sum(foo(:).*foob(:))./(sum(foo(:))+sum(foob(:)));
    end
end

subplot(2,6,[3:4]);
imagesc(metric,[0,1]);xlabel('test');ylabel('master')
title('Dice of Global Fingerprint');colormap(jet);colorbar
freezeColors

% Find mutual maxima, loop through MB first
for j=1:NnetsB
    [M,idx]=max(metric(:,j)); % max for test
    [~,idx2]=max(metric(idx,:)); % max for master
    if idx2==j 
        MBoutKey(j,:)=[j,idx,M]; % [test,master]
    else
        TryAgain=cat(1,TryAgain,[j,idx,idx2]);
    end
end
MAused=MBoutKey(MBoutKey(:,2)~=0,2);

subplot(2,6,[3:4]);hold on
for j=1:length(MAused)
    idx=find(MBoutKey(:,2)==MAused(j));
    plot(MBoutKey(idx,1),MBoutKey(idx,2),'w*')
end

%% (2) Find largest instance of each network
if length(TryAgain)>1
disp('<< Testing Individual kden Networks')
for j=1:NnetsA
    [~,Cidx]=max(sum(MA==valsA(j)));
    NmaxKDenA(j)=Cidx;
    masterA(:,j)=MA(:,NmaxKDenA(j))==valsA(j);
end
for j=1:NnetsB
    [~,Cidx]=max(sum(MB==valsB(j)));
    NmaxKDenB(j)=Cidx;
    testB(:,j)=MB(:,NmaxKDenB(j))==valsB(j);
    for k=1:NnetsA % DICE overlap
        metric2(k,j)=2*sum(masterA(:,k).*testB(:,j))./...
            (sum(masterA(:,k))+sum(testB(:,j)));
    end
end


subplot(2,6,[9:10]);
imagesc(metric2,[0,1]);xlabel('test');ylabel('master')
title('Dice of Kden-specific match');colormap(jet);colorbar
freezeColors

for j=1:size(TryAgain,1)
    [M,idx]=max(metric2(:,TryAgain(j,1))); % max for test
    [~,idx2]=max(metric2(idx,:)); % max for master
    if (idx2==TryAgain(j,1) && ~ismember(idx,MAused))
        MBoutKey(TryAgain(j,1),:)=[TryAgain(j,1),idx,M]; % [test,master]
        TryAgain(j,1)=-TryAgain(j,1);
    end
end
MAused=MBoutKey(MBoutKey(:,2)~=0,2);
end

subplot(2,6,[9:10]);hold on
new=find(TryAgain(:,1)<0);
for j=1:length(new)
    plot(MBoutKey(abs(TryAgain(new(j),1)),1),...
        MBoutKey(abs(TryAgain(new(j),1)),2),'w*')
end

temp=metric;
temp(:,MBoutKey(:,1)~=0)=0;
temp2=metric2;
temp2(:,MBoutKey(:,1)~=0)=0;

% figure;
% subplot(2,1,1);
% imagesc(temp);xlabel('test');ylabel('master')
% title('Global Fingerprint')
% subplot(2,1,2);
% imagesc(temp2);xlabel('test');ylabel('master')
% title('Kden-specific match')
% 

%% (3) Try 2nd and 3rd most similar, keep if unused. Else, Set as new value.
TryAgain(TryAgain(:,1)<0,:)=[];
if length(TryAgain)>1
disp('<< Testing Remaining Networks for match or new case')
for j=1:size(TryAgain,1)
    
    % Try 3 times to match with unique matrix to peak value
    M1=max(temp2(:));
    [idxM,idxT]=find(temp2==M1);                                   % try 1
    if ~ismember(idxM,MAused)
        if temp2(idxM,idxT)>=Dth
            MBoutKey(idxT,:)=[idxT,idxM,temp2(idxM,idxT)];
            temp2(:,idxT)=0;
        end
    else temp2(idxM,idxT)=0;
        idxM2=find(temp2(:,idxT)==max(temp2(:,idxT)));             % try 2
        if ~ismember(idxM2,MAused)
            if temp2(idxM2,idxT)>=Dth
                MBoutKey(idxT,:)=[idxT,idxM2,temp2(idxM2,idxT)];
                temp2(:,idxT)=0;
            end
        else temp2(idxM2,idxT)=0;
            idxM3=find(temp2(:,idxT)==max(temp2(:,idxT)));         % try 3
            if ~ismember(idxM3,MAused)
                if temp2(idxM3,idxT)>=Dth
                    MBoutKey(idxT,:)=[idxT,idxM3,temp2(idxM3,idxT)];
                    temp2(:,idxT)=0;
                end
            else temp2(idxM3,idxT)=0; % Set new network
                
    % If still no unique match, check dice with best match and
    % set as alternative to that network. Else, set to new
    % unknown network. Either case, set new colormap value
                n=n+1;
                if (NetAmax+n)>=valsA(end)
error('Max network number reached. Increase Unclassified Index.')
                else
                    MBoutKey(idxT,:)=[idxT,NetAmax+n,0];
                    temp2(:,idxT)=0;
                    if M1>=Dth % Alt to existing net
                        MBoutKey(idxT,:)=[idxT,NetAmax+n,idxM];
                        Cmap2(NetAmax+n,:)=cmap(idxM,:).*(0.8^n);
                    else        % New Unclassified network
                        MBoutKey(idxT,:)=[idxT,NetAmax+n,0];
                        Cmap2(NetAmax+n,:)=[1,1,1].*rand;
                    end
                end
            end
        end
    end
end
MAused=MBoutKey(MBoutKey(:,2)~=0,2);
end

% Last check that something didn't slip through (may have zero overlap)
uhoh=find(MBoutKey(:,1)==0);
if any(uhoh)
    for j=1:size(uhoh,1)
        n=n+1;
        if (NetAmax+n)>=valsA(end)
error('Max network number reached. Increase Unclassified Index.')
        else% New Unclassified network
            MBoutKey(uhoh(j),:)=[uhoh(j),NetAmax+n,0];
            Cmap2(NetAmax+n,:)=[1,1,1].*rand;
        end
    end
end

%% Relabel networks
disp('<< Relabeling')
for j=1:size(MBoutKey,1)
    MBout(MB==MBoutKey(j,1))=MBoutKey(j,2);
end

subplot(2,6,[5,11]);
imagesc(sortrows(MA,Nkden));colormap(Cmap2);hold on;axis off
plot([Nkden-.5,Nkden-0.5],[0,Nroi],'k');
plot([Nkden+.5,Nkden+0.5],[0,Nroi],'k')
title('Master key')
freezeColors
subplot(2,6,[6,12]);
imagesc(sortrows(MBout,Nkden));colormap(Cmap2);hold on;axis off
plot([Nkden-.5,Nkden-0.5],[0,Nroi],'k');
plot([Nkden+.5,Nkden+0.5],[0,Nroi],'k')
title('Re-Sorted Test set')
freezeColors

subplot(2,6,[3:4]);colormap(jet);colorbar;freezeColors
subplot(2,6,[9:10]);colormap(jet);colorbar;freezeColors
set(gcf,'Color','w')