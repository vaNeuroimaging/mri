%This script originally sent from Xiong jiang in Max Riesenhuber's lab.
%Used to check for SNR in EPI images in out-of-brain voxels
%  determined by all voxels with -86>x<86.
%Designed for his group of 11 subjects with 6 sessions each.
%Output is an 11x6x3 matrix, with the third dimension indicating
%  mean voxel intensity, and the std of all mean voxel intensities,
%  and the mean/std (which is the proposed approximation of SNR).

%subjSNR = zeros(11,6,3);
subs = {'161','189','110' '166' '102' '113' '172' '187' '207' '211' '214' '229' '232' '233' '242' '250' '254' '255' '272' '274' '101' '118' '120' '122' '125' '127' '132' '138' '147' '150' '151' '154' '156' '159' '160' '162' '181' '182' '199' '202' '215' '221' '225' '269'};
runs = {'Nback','Rest'};
subjSNR = zeros(length(subs),length(runs),3);

for i=1:length(subs)
    subs{i}
    for j=1:length(runs)
        cwd=['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subs{i} '/' runs{j} '/'];
        cd(cwd);
        clear files;
            files=dir(['D*.img']);
        len=length(files);
        clear xyz;
        clear y;
        for k=1:len
            clear V;
            V=spm_vol(files(k).name);
            [y(:,:,:,k) xyz] = spm_read_vols(V);
        end
        
        clear yy;
        yy = reshape(y, size(y,1)*size(y,2)*size(y,3), size(y,4));
        clear gs;
        gs = find(abs(xyz(1,:) > 86));
        subjSNR(i,j,1) = mean(mean(yy(gs,:)));
        subjSNR(i,j,2) = std(mean(yy(gs,:)));
        subjSNR(i,j,3) = subjSNR(i,j,1) / subjSNR(i,j,2);
    end
end

subjSNR(:,:,3)

save('/fmri/data3/Evan/Gene-Rest-Nback/Analysis/subjSNR.mat', 'subjSNR');
        
    

