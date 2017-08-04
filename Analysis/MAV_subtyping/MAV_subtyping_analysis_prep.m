subjects = textread('/home/data/Analysis/MAV_subtyping/MAV_subjectlist.txt','%s');
communities = load('/home/data/atlases/Group_parcellation/Parcel_Communities.txt');
communityIDs = unique(communities); communityIDs(communityIDs<1) = [];

%[outcomesubs,outcomes_init] = textread('/home/data/Analysis/MAV_subtyping/MAV_BAI.txt','%s%f');
outcomes = load('/home/data/Analysis/MAV_subtyping/MAV_PCL.txt');
%outcomes = ones(length(subjects),1) .* NaN;

tracts = load_untouch_nii_2D('/home/data/Analysis/MAV_analysis/JHU-ICBM-labels-2mm_selected.nii.gz');
%tracts = load_untouch_nii_2D('JHU-ICBM-labels-2mm-ero.nii.gz');
tractIDs = unique(tracts.img); tractIDs(tractIDs<1) = []; tractIDs(tractIDs==6) = []; tractIDs(tractIDs==47) = []; tractIDs(tractIDs==48) = [];
tractFAvals = zeros(length(subjects),length(tractIDs));

for s = 1:length(subjects)
    
    %matchind = strmatch(subjects{s}, outcomesubs);
    %outcomes(s,1) = outcomes_init(matchind);
    
    corrmat = smartload(['/home/data/subjects/' subjects{s} '/connectome/RSFC_Parcels_LR_corr.mat']);
    
    %allfeatures(s,:) = corrmat(triu(true(size(corrmat)),1));
    
    
    
    
    
    communitymat = zeros(length(communityIDs));
    
    for c1 = 1:length(communityIDs)
        for c2 = c1 : length(communityIDs)
            indsmat = false(size(corrmat));
            indsmat(communities==communityIDs(c1),communities==communityIDs(c2)) = true;
            indsmat = indsmat & triu(true(size(corrmat)),1);
            
            communitymat(c1,c2) = mean(corrmat(indsmat));
            
        end
    end
    
    
    FAdata = load_untouch_nii_2D(['/home/data/subjects/' subjects{s} '/DTI/DTI_avg_ec_FA_MNI.nii.gz']);
    
    WMmask = load_untouch_nii_2D(['/home/data/subjects/' subjects{s} '/DTI/aparc+aseg_cerebralwm_mask_222.nii.gz']);
    
    for IDnum = 1:length(tractIDs)
        tractFAvals(s,IDnum) = mean(FAdata.img((tracts.img==tractIDs(IDnum)) & logical(WMmask.img)));
    end
    
    allfeatures(s,:) = [communitymat(triu(true(size(communitymat)),1))' tractFAvals(s,:)];
    
end

allfeatures = [outcomes allfeatures];

%dlmwrite('/home/data/Analysis/MAV_subtyping/RSFC_features.txt',allfeatures,'delimiter',',');
save('/home/data/Analysis/MAV_subtyping/RSFC_features.mat','allfeatures')





