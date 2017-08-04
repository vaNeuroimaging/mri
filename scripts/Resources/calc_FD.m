function calc_FD(motionparamsfile)

brainradius = 50;
contiguoustimepoints = 5;
FDthresh = .2;

if iscell(motionparamsfile)
    for i = 1:length(motionparamsfile)
        temp_params = load(motionparamsfile{i});
        if i==1
            thisrun_params = temp_params;
        else
            thisrun_params = [thisrun_params;temp_params];
        end
    end
else
    thisrun_params = load(motionparamsfile);
end

thisrun_rot = thisrun_params(:,1:3);
thisrun_rot_mm = thisrun_rot * brainradius;
thisrun_params(:,1:3) = thisrun_rot_mm;
thisrun_params_delta = [zeros(1,6) ; diff(thisrun_params)];

%demean and detrend motion regressors
thisrun_params_demeandetrend = demean_detrend(thisrun_params');
thisrun_params_delta_demeandetrend = demean_detrend(thisrun_params_delta');


%Calculate FD and tmask
thisrun_FD = [sum(abs(thisrun_params_delta),2)]';
thisrun_tmask = (thisrun_FD') < FDthresh;

%Apply contiguous frame criterion
chunkedges = [1; diff([thisrun_tmask;0])];
chunkstarts = find(chunkedges==1); chunkstops = (find(chunkedges==-1) -1);
chunksizes = zeros(length(chunkstarts),1);
chunkID = zeros(size(thisrun_tmask));
for chunknum = 1:length(chunkstarts)
    chunkID(chunkstarts(chunknum) : chunkstops(chunknum)) = chunknum;
    chunksizes(chunknum) = nnz(chunkID==chunknum);
    if chunksizes(chunknum) < contiguoustimepoints
        thisrun_tmask(chunkID==chunknum) = 0;
    end
end

string = ['Mean FD: ' num2str(mean(thisrun_FD)) '; ' num2str(nnz(thisrun_tmask) / numel(thisrun_tmask) * 100) '% (' num2str(nnz(thisrun_tmask)) ' of ' num2str(numel(thisrun_tmask)) ') frames retained'];
disp(string)
figure;
plot(thisrun_FD,'k')
hold on
plot((ones(length(thisrun_FD),1) .* FDthresh),'--r')
nantmask = single(~thisrun_tmask); nantmask(logical(nantmask)) = NaN;
plot(nantmask+.001,'b','Linewidth',5)
title(string)