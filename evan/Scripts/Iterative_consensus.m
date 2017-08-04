%pool = parpool(20);
xdistance = 30;
thresholdarray = .01:.005:.1;
outdir = '/home/data/evan/Temp/testing_infomap/';
numiters = 20;
noise_mag = .3;

for iter = 1:numiters
    disp(iter)
    this_corrmat = corrmat + rand(size(corrmat)) * noise_mag;
    
    Run_Infomap_nopar(this_corrmat, distances, xdistance, thresholdarray, 0, outdir)
    assignments = load('rawassn.txt');
    agreementmat = ones(size(assignments,1)) - squareform(pdist(assignments,'hamming'));
    
    delete('rawassn.txt');
    Run_Infomap_nopar(agreementmat, ones(size(agreementmat)), 0, 1, 0, outdir)
    
    assignments2(:,iter) = load('rawassn.txt');
    
    
end

agreementmat = ones(size(assignments2,1)) - squareform(pdist(assignments2,'hamming'));
uncertainmat = (agreementmat > 0) & (agreementmat < 1);

% uncertainmat = (agreementmat > 0) & (agreementmat < 1);
% 
% while any(uncertainmat(:))
%     
%     delete('rawassn.txt');
%     Run_Infomap_nopar(agreementmat, zeros(size(agreementmat)), 1000, ones(1,numiters), 0, outdir)
%     assignments = load('rawassn.txt');
%     agreementmat = ones(size(assignments,1)) - squareform(pdist(assignments,'hamming'));
%     uncertainmat = (agreementmat > 0) & (agreementmat < 1);
%     
% end