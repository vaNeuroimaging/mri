cd /fmri/data3/Evan/Gene-Rest-Nback/Scripts

warning off

subs = {'101','102','113','118','120','122','125','127','132','138','147','150','151','154','156','159','160','161','162','166','172','181','182','187','202','207','211','214','215','221','225','229','232','233','242','250','254','255','272','274','112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','247','295','297','300','305','309','334','339','340','343','359','362','383','395','396','374'};

genematrix = textread('/fmri/data3/Evan/Gene-Rest-Nback/Data/Genes.txt');

for sub = 1:length(subs);
    allCOMTgenotype(sub) = genematrix(find(genematrix==str2num(subs{sub})),3); 
    allDATgenotype(sub) = genematrix(find(genematrix==str2num(subs{sub})),2);
    allDRD4genotype(sub) = genematrix(find(genematrix==str2num(subs{sub})),4);
    allBDNFgenotype(sub) = genematrix(find(genematrix==str2num(subs{sub})),5);
    allSRTTgenotype(sub) = genematrix(find(genematrix==str2num(subs{sub})),6);
end

DATsubs = subs(find(allDATgenotype>0));
DATgenotype = allDATgenotype(find(allDATgenotype>0));

DRD4subs = subs(find(allDRD4genotype>0));
DRD4genotype = allDRD4genotype(find(allDRD4genotype>0));

BDNFsubs = subs(find(allBDNFgenotype>0 & allBDNFgenotype<3));
BDNFgenotype = allBDNFgenotype(find(allBDNFgenotype>0 & allBDNFgenotype<3));

SRTTsubs = subs(find(allSRTTgenotype>0));
SRTTgenotype = allSRTTgenotype(find(allSRTTgenotype>0));

BDNFmetxsubs = subs(find(allBDNFgenotype>0));
BDNFmetxgenotype = allBDNFgenotype(find(allBDNFgenotype>0));
BDNFmetxgenotype(find(BDNFmetxgenotype==3)) = 2;

COMTwDATgenotype = allCOMTgenotype(find(allDATgenotype>0));

BDNFwDATsubs = intersect(DATsubs,BDNFsubs);
BDNFwDATgenotype = allBDNFgenotype(find(allBDNFgenotype>0 & allBDNFgenotype<3 & allDATgenotype>0));
DATwBDNFgenotype = allDATgenotype(find(allBDNFgenotype>0 & allBDNFgenotype<3 & allDATgenotype>0));

COMTwBDNFgenotype = allCOMTgenotype(find(allBDNFgenotype>0 & allBDNFgenotype<3));




behavname = {'RT','RT3B','Acc','Acc100%','Acc3B','ICV','ICV_3Back','Inattentive','H-I','BIS'};

behavmatrix = textread('/fmri/data3/Evan/Gene-Rest-Nback/Analysis/Behavior/AllMeasures.txt');

for behavnum = 1:length(behavname)
    for sub = 1:length(subs)
        behavior{behavnum}(sub) = behavmatrix(find(behavmatrix(:,1)==str2num(subs{sub})),behavnum+1);
    end
end




DATcounter = [1 1];
load SPM8_2Sample_Ttest_template;
for subnum = 1:length(DATsubs)
    subj = DATsubs{subnum};
    
    if DATgenotype(subnum) == 1
        matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1{DATcounter(DATgenotype(subnum))} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/Nonstationarity/dACC.nii,1'];
    else
        matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2{DATcounter(DATgenotype(subnum))} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/Nonstationarity/dACC.nii,1'];
    end
    DATcounter(DATgenotype(subnum)) = DATcounter(DATgenotype(subnum))+1;
    
end
matlabbatch{1}.spm.stats.factorial_design.dir{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Nonstationarity/dACC_DAT/';

try rmdir(matlabbatch{1}.spm.stats.factorial_design.dir{1},'s'); catch; end
mkdir(matlabbatch{1}.spm.stats.factorial_design.dir{1});

save([matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat'], 'matlabbatch');
spm_jobman('run',[matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat']);

clear matlabbatch




% load SPM8_1Sample_Ttest_withCov_template;
% subcounter = 0;
% for subnum = 1:length(subs)
%     subj = subs{subnum};
%        subcounter = subcounter+1;
%         
%        matlabbatch{1}.spm.stats.factorial_design.des.t1.scans{subcounter} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/Nonstationarity/dACC.nii,1'];
%        matlabbatch{1}.spm.stats.factorial_design.cov.c(subcounter) = behavior{10}(subnum);
% end
% matlabbatch{1}.spm.stats.factorial_design.dir{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Nonstationarity/dACC_vs_BIS/';
% 
% try rmdir(matlabbatch{1}.spm.stats.factorial_design.dir{1},'s'); catch; end
% mkdir(matlabbatch{1}.spm.stats.factorial_design.dir{1});
% 
% save([matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat'], 'matlabbatch');
% spm_jobman('run',[matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat']);
% 
% clear matlabbatch




% BDNFCounter = [1 1];
% load SPM8_LevelXDAT_template
% 
% matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).name = 'Side';
% matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).name = 'BDNF';
% 
% 
% for subnum = 1:length(BDNFmetxsubs)
%     subj = BDNFmetxsubs{subnum};
%     
%     matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(BDNFmetxgenotype(subnum)).levels = [1;BDNFmetxgenotype(subnum)];
%     matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(BDNFmetxgenotype(subnum)).scans{BDNFCounter(BDNFmetxgenotype(subnum)),1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/SPM8/R_Hipp_Rest/con_0001.img,1'];
%     matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(BDNFmetxgenotype(subnum)+2).levels = [2;BDNFmetxgenotype(subnum)];
%     matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(BDNFmetxgenotype(subnum)+2).scans{BDNFCounter(BDNFmetxgenotype(subnum)),1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/SPM8/L_Hipp_Rest/con_0001.img,1'];    
%     
%     BDNFCounter(BDNFmetxgenotype(subnum)) = BDNFCounter(BDNFmetxgenotype(subnum))+1;
%     
% end
% matlabbatch{1}.spm.stats.factorial_design.dir{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Seed-based/BDNFmetxXSide_Hipp_Rest_78/';
% 
% try;rmdir(matlabbatch{1}.spm.stats.factorial_design.dir{1},'s');catch;end
% mkdir(matlabbatch{1}.spm.stats.factorial_design.dir{1});
% 
% save([matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat'], 'matlabbatch');
% spm_jobman('run',[matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat']);
% 
% clear matlabbatch





%%COMTXBDNF
% genecounter = [1 1;1 1;1 1];
% load('SPM8_COMTxDAT_template');
% for subnum = 1:length(BDNFsubs)
%     subj = BDNFsubs{subnum};
%     matlabbatch{1}.spm.stats.factorial_design.des.fd.icell((COMTwBDNFgenotype(subnum)-1)*2 + BDNFgenotype(subnum)).levels = [COMTwBDNFgenotype(subnum) BDNFgenotype(subnum)];
%     matlabbatch{1}.spm.stats.factorial_design.des.fd.icell((COMTwBDNFgenotype(subnum)-1)*2 + BDNFgenotype(subnum)).scans{genecounter(COMTwBDNFgenotype(subnum),BDNFgenotype(subnum)),1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/SPM8/L_Hipp_Rest/con_0001.img,1'];
%     genecounter(COMTwBDNFgenotype(subnum),BDNFgenotype(subnum)) = genecounter(COMTwBDNFgenotype(subnum),BDNFgenotype(subnum))+1;
% end
% matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).name = 'BDNF';
% matlabbatch{1}.spm.stats.factorial_design.dir{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Seed-based/L_Hipp_Rest_COMTXBDNF/';
% try;rmdir(matlabbatch{1}.spm.stats.factorial_design.dir{1},'s');catch;end
% mkdir(matlabbatch{1}.spm.stats.factorial_design.dir{1});
% save([matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat'], 'matlabbatch');
% spm_jobman('run',[matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat']);
% clear matlabbatch





% load SPM8_1Sample_Ttest_template;
% for subnum = 1:length(subs)
%     subj = subs{subnum};
%     
%     matlabbatch{1}.spm.stats.factorial_design.des.t1.scans{subnum,1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/SPM8/Negative_DL6_Rest/con_0001.img,1'];
% 
% end
% matlabbatch{1}.spm.stats.factorial_design.dir{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Seed-based/Negative_DL6_Rest';
% 
% try rmdir(matlabbatch{1}.spm.stats.factorial_design.dir{1},'s'); catch; end
% mkdir(matlabbatch{1}.spm.stats.factorial_design.dir{1});
% 
% save([matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat'], 'matlabbatch');
% spm_jobman('run',[matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat']);
% 
% clear matlabbatch




%%DTI
% load SPM8_1Sample_Ttest_withCov_template;
% allTPNvsTNN = textread('/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ICA/ICA_61DAT_unconstrained.gica/TPNvsTNN.txt');
% subcounter = 0;
% for subnum = 1:length(DATsubs)
%     subj = DATsubs{subnum};
%     
%     TPNvsTNN = allTPNvsTNN((find(allTPNvsTNN(:,1)==str2num(subj))),2);
%     
%     if ~strcmp(subj,'166')
%        subcounter = subcounter+1;
%         
%        matlabbatch{1}.spm.stats.factorial_design.des.t1.scans{subcounter} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/TBSS/BDNF_77/stats/' subj '_skeletonized.nii,1'];
%        matlabbatch{1}.spm.stats.factorial_design.cov.c(subcounter) = TPNvsTNN;
%         
%     end
%     
% end
% matlabbatch{1}.spm.stats.factorial_design.dir{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/TBSS/BDNF_77/Group_Analysis/FCcorrelation_60/';
% 
% try rmdir(matlabbatch{1}.spm.stats.factorial_design.dir{1},'s'); catch; end
% mkdir(matlabbatch{1}.spm.stats.factorial_design.dir{1});
% 
% save([matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat'], 'matlabbatch');
% spm_jobman('run',[matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat']);
% 
% clear matlabbatch







% %BDNFxDAT
% genecounter = [1 1;1 1];
% load('SPM8_COMTxDAT_template');
% for subnum = 1:length(BDNFwDATsubs)
%     subj = BDNFwDATsubs{subnum};
%     
%     matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).levels = 2;
%     matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).name = 'BDNF';
%     matlabbatch{1}.spm.stats.factorial_design.des.fd.icell = matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(1:4);
%     
%     matlabbatch{1}.spm.stats.factorial_design.des.fd.icell((BDNFwDATgenotype(subnum)-1)*2 + DATwBDNFgenotype(subnum)).levels = [BDNFwDATgenotype(subnum) DATwBDNFgenotype(subnum)];
%     matlabbatch{1}.spm.stats.factorial_design.des.fd.icell((BDNFwDATgenotype(subnum)-1)*2 + DATwBDNFgenotype(subnum)).scans{genecounter(BDNFwDATgenotype(subnum),DATwBDNFgenotype(subnum)),1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/SPM8/L_Hipp_Rest/con_0001.img,1'];
%     genecounter(BDNFwDATgenotype(subnum),DATwBDNFgenotype(subnum)) = genecounter(BDNFwDATgenotype(subnum),DATwBDNFgenotype(subnum))+1;
% end
% matlabbatch{1}.spm.stats.factorial_design.dir{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Seed-based/L_Hipp_Rest_BDNFxDAT/';
% try;rmdir(matlabbatch{1}.spm.stats.factorial_design.dir{1},'s');catch;end
% mkdir(matlabbatch{1}.spm.stats.factorial_design.dir{1});
% save([matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat'], 'matlabbatch');
% spm_jobman('run',[matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat']);
% clear matlabbatch




% %%DTI
% BDNFcounter = [1 1];
% load SPM8_2Sample_Ttest_template;
% for subnum = 1:length(BDNFmetxsubs)
%     subj = BDNFmetxsubs{subnum};
%     
%     if ~strcmp(subj,'166')
%         
%         if BDNFmetxgenotype(subnum) == 1
%             matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1{BDNFcounter(BDNFmetxgenotype(subnum))} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/TBSS/BDNF_77/stats/' subj '_skeletonized.nii,1'];
%             
%         else
%             matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2{BDNFcounter(BDNFmetxgenotype(subnum))} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/TBSS/BDNF_77/stats/' subj '_skeletonized.nii,1'];
%             
%         end
%         BDNFcounter(BDNFmetxgenotype(subnum)) = BDNFcounter(BDNFmetxgenotype(subnum))+1;
%     end
%     
% end
% matlabbatch{1}.spm.stats.factorial_design.dir{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/TBSS/BDNF_77/Group_Analysis/';
% 
% try rmdir(matlabbatch{1}.spm.stats.factorial_design.dir{1},'s'); catch; end
% mkdir(matlabbatch{1}.spm.stats.factorial_design.dir{1});
% 
% save([matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat'], 'matlabbatch');
% spm_jobman('run',[matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat']);
% 
% clear matlabbatch
% 
% 
% 
% %%VBM group analysis
% BDNFcounter = [1 1];
% load SPM8_2Sample_Ttest_withCov_template;
% for subnum = 1:length(BDNFmetxsubs)
%     subj = BDNFmetxsubs{subnum};
%     
%     if BDNFmetxgenotype(subnum) == 1
%         matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1{BDNFcounter(BDNFmetxgenotype(subnum))} = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subj '/Struct/smwrc1Registered_MPRAGE_thr3.nii,1'];
%         
%         t = VBM_get_TCV(['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subj '/Struct/smwrc1Registered_MPRAGE_thr3.nii,1'; '/fmri/data3/Evan/Gene-Rest-Nback/Data/' subj '/Struct/smwrc1Registered_MPRAGE_thr3.nii,1']);
%         
%         TCV_gr1(BDNFcounter(BDNFmetxgenotype(subnum)),1) = sum(t);
%         
%     else
%         matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2{BDNFcounter(BDNFmetxgenotype(subnum))} = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subj '/Struct/smwrc1Registered_MPRAGE_thr3.nii,1'];
%         
%         t = VBM_get_TCV(['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subj '/Struct/smwrc1Registered_MPRAGE_thr3.nii,1'; '/fmri/data3/Evan/Gene-Rest-Nback/Data/' subj '/Struct/smwrc1Registered_MPRAGE_thr3.nii,1']);
%         
%         TCV_gr2(BDNFcounter(BDNFmetxgenotype(subnum)),1) = sum(t);
%         
%     end
%     BDNFcounter(BDNFmetxgenotype(subnum)) = BDNFcounter(BDNFmetxgenotype(subnum))+1;
%     
% end
% matlabbatch{1}.spm.stats.factorial_design.dir{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/VBM/Group_Analysis/BDNF_77/';
% 
% matlabbatch{1}.spm.stats.factorial_design.cov.c = [TCV_gr1;TCV_gr2];
% matlabbatch{1}.spm.stats.factorial_design.cov.cname = 'TCV';
% 
% try rmdir(matlabbatch{1}.spm.stats.factorial_design.dir{1},'s'); catch; end
% mkdir(matlabbatch{1}.spm.stats.factorial_design.dir{1});
% 
% save([matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat'], 'matlabbatch');
% spm_jobman('run',[matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat']);
% 
% clear matlabbatch



% BDNFcounter = [1 1];
% load SPM8_2Sample_Ttest_template;
% for subnum = 1:length(BDNFmetxsubs)
%     subj = BDNFmetxsubs{subnum};
%     
%     if ~strcmp(subj,'166')
%        
%         if BDNFmetxgenotype(subnum) == 1
%             matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1{BDNFcounter(BDNFmetxgenotype(subnum))} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/SPM8/R_Hipp_Rest/con_0001.img,1'];
%             
%         else
%             matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2{BDNFcounter(BDNFmetxgenotype(subnum))} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/SPM8/R_Hipp_Rest/con_0001.img,1'];
%             
%         end
%         BDNFcounter(BDNFmetxgenotype(subnum)) = BDNFcounter(BDNFmetxgenotype(subnum))+1;
%     
%     end
%     
% end
% matlabbatch{1}.spm.stats.factorial_design.dir{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Seed-based/R_Hipp_Rest_BDNF_76_metx/';
% 
% try rmdir(matlabbatch{1}.spm.stats.factorial_design.dir{1},'s'); catch; end
% mkdir(matlabbatch{1}.spm.stats.factorial_design.dir{1});
% 
% save([matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat'], 'matlabbatch');
% spm_jobman('run',[matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat']);
% 
% clear matlabbatch


% genecounter = [1 1;1 1;1 1];
% load('SPM8_COMTxDAT_template');
% for subnum = 1:length(DATsubs)
%     subj = DATsubs{subnum};
%     matlabbatch{1}.spm.stats.factorial_design.des.fd.icell((COMTwDATgenotype(subnum)-1)*2 + DATgenotype(subnum)).levels = [COMTwDATgenotype(subnum) DATgenotype(subnum)];
%     matlabbatch{1}.spm.stats.factorial_design.des.fd.icell((COMTwDATgenotype(subnum)-1)*2 + DATgenotype(subnum)).scans{genecounter(COMTwDATgenotype(subnum),DATgenotype(subnum)),1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/SPM8/DualRegression_Nback_ICA76/beta_0009.img,1'];
%     genecounter(COMTwDATgenotype(subnum),DATgenotype(subnum)) = genecounter(COMTwDATgenotype(subnum),DATgenotype(subnum))+1;
% end
% matlabbatch{1}.spm.stats.factorial_design.dir{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/DualRegression/COMTxDAT_DR_Nback_DMN_76/';
% try;rmdir(matlabbatch{1}.spm.stats.factorial_design.dir{1},'s');catch;end
% mkdir(matlabbatch{1}.spm.stats.factorial_design.dir{1});
% save([matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat'], 'matlabbatch');
% spm_jobman('run',[matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat']);
% clear matlabbatch





% load SPM8_Paired_Ttest_template;
% for subnum = 1:length(DATsubs)
%     subj = DATsubs{subnum};
%     
%     matlabbatch{1}.spm.stats.factorial_design.des.pt.pair(subnum) = matlabbatch{1}.spm.stats.factorial_design.des.pt.pair(1);
%     
%     matlabbatch{1}.spm.stats.factorial_design.des.pt.pair(subnum).scans{1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/SPM8/R_IPL_Nback/con_0001.img,1'];
%     matlabbatch{1}.spm.stats.factorial_design.des.pt.pair(subnum).scans{2} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/SPM8/R_IPL_Rest/con_0001.img,1'];
% 
% end
% matlabbatch{1}.spm.stats.factorial_design.dir{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Seed-based/R_IPL_Nback_vs_Rest_61';
% 
% try rmdir(matlabbatch{1}.spm.stats.factorial_design.dir{1},'s'); catch; end
% mkdir(matlabbatch{1}.spm.stats.factorial_design.dir{1});
% 
% save([matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat'], 'matlabbatch');
% spm_jobman('run',[matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat']);
% 
% clear matlabbatch




% DATCounter = [1 1];
% load SPM8_LevelXDAT_template
% for subnum = 1:length(DATsubs)
%     subj = DATsubs{subnum};
%     
%     matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(DATgenotype(subnum)).levels = [1;DATgenotype(subnum)];
%     matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(DATgenotype(subnum)).scans{DATCounter(DATgenotype(subnum)),1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/SPM8/ActivatedStriatumgating_wTPN//Highhalf/con_0001.img,1'];
%     matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(DATgenotype(subnum)+2).levels = [2;DATgenotype(subnum)];
%     matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(DATgenotype(subnum)+2).scans{DATCounter(DATgenotype(subnum)),1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/SPM8/ActivatedStriatumgating_wTPN//Lowhalf/con_0001.img,1'];    
%     
%     DATCounter(DATgenotype(subnum)) = DATCounter(DATgenotype(subnum))+1;
%     
% end
% matlabbatch{1}.spm.stats.factorial_design.dir{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Seed-based/ActivatedStriatumgating_wTPNxDAT_59/';
% 
% try;rmdir(matlabbatch{1}.spm.stats.factorial_design.dir{1},'s');catch;end
% mkdir(matlabbatch{1}.spm.stats.factorial_design.dir{1});
% 
% save([matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat'], 'matlabbatch');
% spm_jobman('run',[matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat']);
% 
% clear matlabbatch










% AccCounter = [1 1];
% load SPM8_LevelXAcc_template
% for subnum = 1:length(DATsubs)
%     subj = DATsubs{subnum};
%     
%     %cell 1 = inaccurate; 2 = accurate
%     matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(Accuracy(subnum)+1).levels = [1;Accuracy(subnum)+1];
%     matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(Accuracy(subnum)+1).scans{AccCounter(Accuracy(subnum)+1),1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/SPM8/Caudategating/Highhalf/con_0001.img,1'];
%     matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(Accuracy(subnum)+3).levels = [2;Accuracy(subnum)+1];
%     matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(Accuracy(subnum)+3).scans{AccCounter(Accuracy(subnum)+1),1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/SPM8/Caudategating/Lowhalf/con_0001.img,1'];    
%     
%     AccCounter(Accuracy(subnum)+1) = AccCounter(Accuracy(subnum)+1)+1;
%     
% end
% matlabbatch{1}.spm.stats.factorial_design.dir{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Seed-based/CaudateLevelXAcc_61/';
% 
% try;rmdir(matlabbatch{1}.spm.stats.factorial_design.dir{1},'s');catch;end
% mkdir(matlabbatch{1}.spm.stats.factorial_design.dir{1});
% 
% save([matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat'], 'matlabbatch');
% spm_jobman('run',[matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat']);
% 
% clear matlabbatch

% DATcounter = [1 1];
% load /fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Seed-based/CaudateLevel_halves.mat;
% for subnum = 1:length(DATsubs)
%     subj = DATsubs{subnum};
%     
%     matlabbatch{1}.spm.stats.factorial_design.des.pt.pair(subnum).scans{1,1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/SPM8/Caudategating_Rest/Highhalf/con_0001.img,1'];
%     matlabbatch{1}.spm.stats.factorial_design.des.pt.pair(subnum).scans{2,1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/SPM8/Caudategating_Rest/Lowhalf/con_0001.img,1 '];
%     
% end
% matlabbatch{1}.spm.stats.factorial_design.dir{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Seed-based/CaudateLevel_halves_Rest/';
% 
% try rmdir(matlabbatch{1}.spm.stats.factorial_design.dir{1},'s'); catch; end
% mkdir(matlabbatch{1}.spm.stats.factorial_design.dir{1});
% 
% save([matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat'], 'matlabbatch');
% spm_jobman('run',[matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat']);
% 
% clear matlabbatch















































% analysestorun = {'Nback_VStr_DATpeak_sphere_LVWMregress'};
% betas = {'1'};
% outputnames = {'Nback_VStr_DATpeak_sphere_LVWMregress'};
% 
% 
% for analysisnum = 1:length(analysestorun);
%     analysisname = analysestorun{analysisnum};
%     beta = betas{analysisnum};
    
    
%     clear matlabbatch
% 

   
% load Group_analysis_Load;
% for subnum = 1:length(subs)
%     subj = subs{subnum};
%     
%     jobs{1}.stats{1}.factorial_design.des.fd.icell(1).scans{subnum,1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/Cond/con_0001.img,1'];
%     jobs{1}.stats{1}.factorial_design.des.fd.icell(1).levels = [1];
% 
%     jobs{1}.stats{1}.factorial_design.des.fd.icell(2).scans{subnum,1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/Cond/con_0002.img,1'];
%     jobs{1}.stats{1}.factorial_design.des.fd.icell(2).levels = [2];
%     
%     jobs{1}.stats{1}.factorial_design.des.fd.icell(3).scans{subnum,1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/Cond/con_0003.img,1'];
%     jobs{1}.stats{1}.factorial_design.des.fd.icell(3).levels = [3];
%     
%     
% end
% jobs{1}.stats{1}.factorial_design.dir{1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Activation/Load40'];
% jobs{1}.stats{2}.fmri_est.spmmat{1} = [jobs{1}.stats{1}.factorial_design.dir{1} '/SPM.mat'];
% jobs{1}.stats{3}.con.spmmat{1} = jobs{1}.stats{2}.fmri_est.spmmat{1};
% 
% try rmdir(jobs{1}.stats{1}.factorial_design.dir{1},'s'); catch; end
% mkdir(jobs{1}.stats{1}.factorial_design.dir{1});
% 
% save([jobs{1}.stats{1}.factorial_design.dir{1} '.mat'], 'jobs', 'jobhelps');
% spm_jobman('run',[jobs{1}.stats{1}.factorial_design.dir{1} '.mat']) 


% genecounter = [1 1;1 1;1 1];
% load('Group_analysis_COMTXDAT');
% for subnum = 1:length(DATsubs)
%     subj = DATsubs{subnum};
%     jobs{1}.stats{1}.factorial_design.des.fd.icell((COMTwDATgenotype(subnum)-1)*2 + DATgenotype(subnum)).scans{genecounter(COMTwDATgenotype(subnum),DATgenotype(subnum)),1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/' analysisname '/beta_000' beta '.img,1'];
%     genecounter(COMTwDATgenotype(subnum),DATgenotype(subnum)) = genecounter(COMTwDATgenotype(subnum),DATgenotype(subnum))+1;
%     jobs{1}.stats{1}.factorial_design.des.fd.icell((COMTwDATgenotype(subnum)-1)*2 + DATgenotype(subnum)).levels = [COMTwDATgenotype(subnum) DATgenotype(subnum)];
% end
% jobs{1}.stats{1}.factorial_design.dir{1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Seed-based/' analysisname '_COMTxDAT'];
% jobs{1}.stats{2}.fmri_est.spmmat{1} = [jobs{1}.stats{1}.factorial_design.dir{1} '/SPM.mat'];
% jobs{1}.stats{3}.con.spmmat{1} = jobs{1}.stats{2}.fmri_est.spmmat{1};
% 
% try rmdir(jobs{1}.stats{1}.factorial_design.dir{1},'s');catch;end
% mkdir(jobs{1}.stats{1}.factorial_design.dir{1});
% 
% save([jobs{1}.stats{1}.factorial_design.dir{1} '.mat'], 'jobs', 'jobhelps');
% spm_jobman('run',[jobs{1}.stats{1}.factorial_design.dir{1} '.mat'])
% 
% clear jobs
% 
% 
% clear jobs
% 
% load('Group_analysis_LinearAllele');
% for subnum = 1:length(DATsubs)
%     subj = DATsubs{subnum};
%     
%     jobs{1}.stats{1}.factorial_design.des.mreg.scans{subnum} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/' analysisname '/beta_000' beta '.img,1'];
%     jobs{1}.stats{1}.factorial_design.des.mreg.mcov.c(subnum) = -(DATgenotype(subnum)-3) + COMTwDATgenotype(subnum);
%     %-(DATgenotype(subnum)-3)*3 + COMTwDATgenotype(subnum);
%     
%    
% end
% jobs{1}.stats{1}.factorial_design.dir{1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Seed-based/' analysisname 'LinearAllele'];
% jobs{1}.stats{2}.fmri_est.spmmat{1} = [jobs{1}.stats{1}.factorial_design.dir{1} '/SPM.mat'];
% jobs{1}.stats{3}.con.spmmat{1} = jobs{1}.stats{2}.fmri_est.spmmat{1};
% 
% try rmdir(jobs{1}.stats{1}.factorial_design.dir{1},'s');catch;end
% mkdir(jobs{1}.stats{1}.factorial_design.dir{1});
% 
% save([jobs{1}.stats{1}.factorial_design.dir{1} '.mat'], 'jobs', 'jobhelps');
% spm_jobman('run',[jobs{1}.stats{1}.factorial_design.dir{1} '.mat'])
% 
% clear jobs



% clear jobs
% 
% genecounter = [1 1];
% load('Group_analysis_DAT_mask');
% for subnum = 1:length(DATsubs)
%     subj = DATsubs{subnum};
%     
%     if DATgenotype(subnum)==1
%         jobs{1}.stats{1}.factorial_design.des.t2.scans1{genecounter(1)} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/' analysisname '/con_000' beta '.img,1'];
%         genecounter(1) = genecounter(1)+1;
%     else
%         jobs{1}.stats{1}.factorial_design.des.t2.scans2{genecounter(2)} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/' analysisname '/con_000' beta '.img,1'];
%         genecounter(2) = genecounter(2)+1;
%     end
%     
% end
% jobs{1}.stats{1}.factorial_design.masking.em{1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/Striatum_dialated.nii,1'];
% jobs{1}.stats{1}.factorial_design.dir{1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Activation/' outputnames{analysisnum} '_DAT'];
% jobs{1}.stats{2}.fmri_est.spmmat{1} = [jobs{1}.stats{1}.factorial_design.dir{1} '/SPM.mat'];
% jobs{1}.stats{3}.con.spmmat{1} = jobs{1}.stats{2}.fmri_est.spmmat{1};
% 
% try rmdir(jobs{1}.stats{1}.factorial_design.dir{1},'s');catch;end
% mkdir(jobs{1}.stats{1}.factorial_design.dir{1});
% 
% save([jobs{1}.stats{1}.factorial_design.dir{1} '.mat'], 'jobs', 'jobhelps');
% spm_jobman('run',[jobs{1}.stats{1}.factorial_design.dir{1} '.mat'])
% 
% clear jobs

% 
% analysisname = 'Slidingwindow_connectivitygating_ACCvsPCC';
% 
% 
% clear jobs
% 
% genecounter = [1 1];
% load('Group_analysis_DAT');
% for subnum = 1:length(DATsubs)
%     subj = DATsubs{subnum};
%     
%     if DATgenotype(subnum)==1
%         jobs{1}.stats{1}.factorial_design.des.t2.scans1{genecounter(1)} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/SPM8/' analysisname '/con_0001.img,1'];
%         genecounter(1) = genecounter(1)+1;
%     else
%         jobs{1}.stats{1}.factorial_design.des.t2.scans2{genecounter(2)} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/SPM8/' analysisname '/con_0001.img,1'];
%         genecounter(2) = genecounter(2)+1;
%     end
%     
% end
% jobs{1}.stats{1}.factorial_design.dir{1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Seed-based/' analysisname  'xDAT'];
% jobs{1}.stats{2}.fmri_est.spmmat{1} = [jobs{1}.stats{1}.factorial_design.dir{1} '/SPM.mat'];
% jobs{1}.stats{3}.con.spmmat{1} = jobs{1}.stats{2}.fmri_est.spmmat{1};
% 
% try rmdir(jobs{1}.stats{1}.factorial_design.dir{1},'s');catch;end
% mkdir(jobs{1}.stats{1}.factorial_design.dir{1});
% 
% save([jobs{1}.stats{1}.factorial_design.dir{1} '.mat'], 'jobs', 'jobhelps');
% spm_jobman('run',[jobs{1}.stats{1}.factorial_design.dir{1} '.mat'])
% 
% clear jobs









% genecounter = [1 1;1 1;1 1];
% load Group_analysis_COMTXDATXLoad;
% for subnum = 1:length(DATsubs)
%     subj = DATsubs{subnum};
%     
%     jobs{1}.stats{1}.factorial_design.des.fd.icell((COMTwDATgenotype(subnum)-1)*6 + (DATgenotype(subnum)-1)*3 + 1).scans{genecounter(COMTwDATgenotype(subnum),DATgenotype(subnum)),1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/' analysisname '/beta_0003.img,1'];
%     jobs{1}.stats{1}.factorial_design.des.fd.icell((COMTwDATgenotype(subnum)-1)*6 + (DATgenotype(subnum)-1)*3 + 1).levels = [COMTwDATgenotype(subnum) DATgenotype(subnum) 1];
% 
%     jobs{1}.stats{1}.factorial_design.des.fd.icell((COMTwDATgenotype(subnum)-1)*6 + (DATgenotype(subnum)-1)*3 + 2).scans{genecounter(COMTwDATgenotype(subnum),DATgenotype(subnum)),1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/' analysisname '/beta_0004.img,1'];
%     jobs{1}.stats{1}.factorial_design.des.fd.icell((COMTwDATgenotype(subnum)-1)*6 + (DATgenotype(subnum)-1)*3 + 2).levels = [COMTwDATgenotype(subnum) DATgenotype(subnum) 2];
% 
%     jobs{1}.stats{1}.factorial_design.des.fd.icell((COMTwDATgenotype(subnum)-1)*6 + (DATgenotype(subnum)-1)*3 + 3).scans{genecounter(COMTwDATgenotype(subnum),DATgenotype(subnum)),1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/' analysisname '/beta_0005.img,1'];
%     jobs{1}.stats{1}.factorial_design.des.fd.icell((COMTwDATgenotype(subnum)-1)*6 + (DATgenotype(subnum)-1)*3 + 3).levels = [COMTwDATgenotype(subnum) DATgenotype(subnum) 3];
% 
%     
%     genecounter(COMTwDATgenotype(subnum),DATgenotype(subnum)) = genecounter(COMTwDATgenotype(subnum),DATgenotype(subnum))+1;
% 
% 
% end
% jobs{1}.stats{1}.factorial_design.dir{1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Seed-based/' analysisname 'COMTxDATxLoad'];
% jobs{1}.stats{2}.fmri_est.spmmat{1} = [jobs{1}.stats{1}.factorial_design.dir{1} '/SPM.mat'];
% jobs{1}.stats{3}.con.spmmat{1} = jobs{1}.stats{2}.fmri_est.spmmat{1};
% 
% try rmdir(jobs{1}.stats{1}.factorial_design.dir{1},'s'); catch; end
% mkdir(jobs{1}.stats{1}.factorial_design.dir{1});
% 
% save([jobs{1}.stats{1}.factorial_design.dir{1} '.mat'], 'jobs', 'jobhelps');
% spm_jobman('run',[jobs{1}.stats{1}.factorial_design.dir{1} '.mat'])


% load Group_analysis;
% for subnum = 1:length(subs)
%     subj = subs{subnum};
%     jobs{1}.stats{1}.factorial_design.des.t1.scans{subnum} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/LinearLoad/con_0001.img,1'];
% end
% jobs{1}.stats{1}.factorial_design.dir{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Activation/LinearLoad';
%   jobs{1}.stats{2}.fmri_est.spmmat{1} = [jobs{1}.stats{1}.factorial_design.dir{1} '/SPM.mat'];
%   jobs{1}.stats{3}.con.spmmat{1} = [jobs{1}.stats{1}.factorial_design.dir{1} '/SPM.mat'];
% 
% try rmdir(jobs{1}.stats{1}.factorial_design.dir{1},'s'); catch; end
% mkdir(jobs{1}.stats{1}.factorial_design.dir{1});
% 
% save(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/LinearLoad.mat'], 'jobs', 'jobhelps');
% spm_jobman('run',['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/LinearLoad.mat'])
% 
% clear jobs

% DATcounter = [1 1];
% load Group_analysis_DAT;
% for subnum = 1:length(DATsubs)
%     subj = DATsubs{subnum};
%     jobs{1}.stats{1}.factorial_design.des.fd.icell(DATgenotype(subnum)).scans{DATcounter(DATgenotype(subnum)),1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/Cond/con_0004.img,1'];
%     DATcounter(DATgenotype(subnum)) = DATcounter(DATgenotype(subnum))+1;
%     jobs{1}.stats{1}.factorial_design.des.fd.icell(DATgenotype(subnum)).levels = DATgenotype(subnum);
% 
% end
% jobs{1}.stats{1}.factorial_design.dir{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Activation/DAT_2v1';
% jobs{1}.stats{2}.fmri_est.spmmat{1} = [jobs{1}.stats{1}.factorial_design.dir{1} '/SPM.mat'];
% jobs{1}.stats{3}.con.spmmat{1} = jobs{1}.stats{2}.fmri_est.spmmat{1};
% 
% try rmdir(jobs{1}.stats{1}.factorial_design.dir{1},'s'); catch; end
% mkdir(jobs{1}.stats{1}.factorial_design.dir{1});
% 
% save(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/DAT_2v1.mat'], 'jobs', 'jobhelps');
% spm_jobman('run',['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/DAT_2v1.mat'])
% 
% clear jobs
% 
% DATcounter = [1 1];
% load Group_analysis_DAT;
% for subnum = 1:length(DATsubs)
%     subj = DATsubs{subnum};
%     jobs{1}.stats{1}.factorial_design.des.fd.icell(DATgenotype(subnum)).scans{DATcounter(DATgenotype(subnum)),1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/Cond/con_0005.img,1'];
%     DATcounter(DATgenotype(subnum)) = DATcounter(DATgenotype(subnum))+1;
%     jobs{1}.stats{1}.factorial_design.des.fd.icell(DATgenotype(subnum)).levels = DATgenotype(subnum);
% 
% end
% jobs{1}.stats{1}.factorial_design.dir{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Activation/DAT_3v1';
% jobs{1}.stats{2}.fmri_est.spmmat{1} = [jobs{1}.stats{1}.factorial_design.dir{1} '/SPM.mat'];
% jobs{1}.stats{3}.con.spmmat{1} = jobs{1}.stats{2}.fmri_est.spmmat{1};
% 
% try rmdir(jobs{1}.stats{1}.factorial_design.dir{1},'s'); catch; end
% mkdir(jobs{1}.stats{1}.factorial_design.dir{1});
% 
% save(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/DAT_3v1.mat'], 'jobs', 'jobhelps');
% spm_jobman('run',['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/DAT_3v1.mat'])
% 
% clear jobs
% 
% DATcounter = [1 1];
% load Group_analysis_DAT;
% for subnum = 1:length(DATsubs)
%     subj = DATsubs{subnum};
%     jobs{1}.stats{1}.factorial_design.des.fd.icell(DATgenotype(subnum)).scans{DATcounter(DATgenotype(subnum)),1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/Cond/con_0006.img,1'];
%     DATcounter(DATgenotype(subnum)) = DATcounter(DATgenotype(subnum))+1;
%     jobs{1}.stats{1}.factorial_design.des.fd.icell(DATgenotype(subnum)).levels = DATgenotype(subnum);
% 
% end
% jobs{1}.stats{1}.factorial_design.dir{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Activation/DAT_3v2';
% jobs{1}.stats{2}.fmri_est.spmmat{1} = [jobs{1}.stats{1}.factorial_design.dir{1} '/SPM.mat'];
% jobs{1}.stats{3}.con.spmmat{1} = jobs{1}.stats{2}.fmri_est.spmmat{1};
% 
% try rmdir(jobs{1}.stats{1}.factorial_design.dir{1},'s'); catch; end
% mkdir(jobs{1}.stats{1}.factorial_design.dir{1});
% 
% save(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/DAT_3v2.mat'], 'jobs', 'jobhelps');
% spm_jobman('run',['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/DAT_3v2.mat'])
% 
% clear jobs

% load LinearCOMT;
% for subnum = 1:length(subs)
%     subj = subs{subnum};
%     jobs{1}.stats{1}.factorial_design.des.mreg.scans{subnum} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/Cond/con_0005.img,1'];
%     jobs{1}.stats{1}.factorial_design.des.mreg.mcov.c(subnum) = COMTgenotype(subnum);
% end
% jobs{1}.stats{1}.factorial_design.dir{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Activation/COMT_IndCondTask_linear';
% jobs{1}.stats{2}.fmri_est.spmmat{1} = [jobs{1}.stats{1}.factorial_design.dir{1} '/SPM.mat'];
% jobs{1}.stats{3}.con.spmmat{1} = jobs{1}.stats{2}.fmri_est.spmmat{1};
% 
% try rmdir(jobs{1}.stats{1}.factorial_design.dir{1},'s');catch;end
% mkdir(jobs{1}.stats{1}.factorial_design.dir{1});
% 
% save([jobs{1}.stats{1}.factorial_design.dir{1} '.mat'], 'jobs', 'jobhelps');
% spm_jobman('run',[jobs{1}.stats{1}.factorial_design.dir{1} '.mat'])
% 
% clear jobs
% 
% COMTcounter = [1 1 1];
% load Group_analysis_COMT;
% for subnum = 1:length(subs)
%     subj = subs{subnum};
%     jobs{1}.stats{1}.factorial_design.des.fd.icell(COMTgenotype(subnum)).scans{COMTcounter(COMTgenotype(subnum)),1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/Cond/con_0005.img,1'];
%     COMTcounter(COMTgenotype(subnum)) = COMTcounter(COMTgenotype(subnum))+1;
%     jobs{1}.stats{1}.factorial_design.des.fd.icell(COMTgenotype(subnum)).levels = COMTgenotype(subnum);
% 
% end
% jobs{1}.stats{1}.factorial_design.dir{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Activation/COMT_IndCondTask';
% jobs{1}.stats{2}.fmri_est.spmmat{1} = [jobs{1}.stats{1}.factorial_design.dir{1} '/SPM.mat'];
% jobs{1}.stats{3}.con.spmmat{1} = jobs{1}.stats{2}.fmri_est.spmmat{1};
% 
% try rmdir(jobs{1}.stats{1}.factorial_design.dir{1},'s'); catch; end
% mkdir(jobs{1}.stats{1}.factorial_design.dir{1});
% 
% save([jobs{1}.stats{1}.factorial_design.dir{1} '.mat'], 'jobs', 'jobhelps');
% spm_jobman('run',[jobs{1}.stats{1}.factorial_design.dir{1} '.mat'])
% 




% 
% load Group_analysis;
% for subnum = 1:length(subs)
%     subj = subs{subnum};
%     jobs{1}.stats{1}.factorial_design.des.t1.scans{subnum} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/con_0003.img,1'];
% end
% jobs{1}.stats{1}.factorial_design.dir{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Task';
%   jobs{1}.stats{2}.fmri_est.spmmat{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Task/SPM.mat';
%   jobs{1}.stats{3}.con.spmmat{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Task/SPM.mat';
% 
% rmdir(jobs{1}.stats{1}.factorial_design.dir{1},'s');
% mkdir(jobs{1}.stats{1}.factorial_design.dir{1});
%
% save(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Task.mat'], 'jobs', 'jobhelps');
% spm_jobman('run',['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Task.mat'])
% 
% clear jobs
% 
% 
% load Group_analysis;
% for subnum = 1:length(subs)
%     subj = subs{subnum};
%     jobs{1}.stats{1}.factorial_design.des.t1.scans{subnum} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/con_0004.img,1'];
% end
% jobs{1}.stats{1}.factorial_design.dir{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Reverse_Task';
%   jobs{1}.stats{2}.fmri_est.spmmat{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Reverse_Task/SPM.mat';
%   jobs{1}.stats{3}.con.spmmat{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Reverse_Task/SPM.mat';
% 
% rmdir(jobs{1}.stats{1}.factorial_design.dir{1},'s');
% mkdir(jobs{1}.stats{1}.factorial_design.dir{1});
%
% save(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Reverse_Task.mat'], 'jobs', 'jobhelps');
% spm_jobman('run',['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Reverse_Task.mat'])
% 
% clear jobs
% 
% 
% 
% load Group_analysis;
% for subnum = 1:length(subs)
%     subj = subs{subnum};
%     jobs{1}.stats{1}.factorial_design.des.t1.scans{subnum} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/con_0001.img,1'];
% end
% jobs{1}.stats{1}.factorial_design.dir{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Load';
%   jobs{1}.stats{2}.fmri_est.spmmat{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Load/SPM.mat';
%   jobs{1}.stats{3}.con.spmmat{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Load/SPM.mat';
% 
% rmdir(jobs{1}.stats{1}.factorial_design.dir{1},'s');
% mkdir(jobs{1}.stats{1}.factorial_design.dir{1});
%
% save(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Load.mat'], 'jobs', 'jobhelps');
% spm_jobman('run',['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Load.mat'])
% 
% clear jobs
% 
% load Group_analysis;
% for subnum = 1:length(subs)
%     subj = subs{subnum};
%     jobs{1}.stats{1}.factorial_design.des.t1.scans{subnum} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/con_0002.img,1'];
% end
% jobs{1}.stats{1}.factorial_design.dir{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Reverse_Load';
%   jobs{1}.stats{2}.fmri_est.spmmat{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Reverse_Load/SPM.mat';
%   jobs{1}.stats{3}.con.spmmat{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Reverse_Load/SPM.mat';
% 
% rmdir(jobs{1}.stats{1}.factorial_design.dir{1},'s');
% mkdir(jobs{1}.stats{1}.factorial_design.dir{1});
%
% save(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Reverse_Load.mat'], 'jobs', 'jobhelps');
% spm_jobman('run',['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Reverse_Load.mat'])
% 
% clear jobs
% 
% COMTcounter = [1 1 1];
% load Group_analysis_COMT;
% for subnum = 1:length(subs)
%     subj = subs{subnum};
%     jobs{1}.stats{1}.factorial_design.des.fd.icell(genotypes(subnum)).scans{COMTcounter(genotypes(subnum)),1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/con_0003.img,1'];
%     COMTcounter(genotypes(subnum)) = COMTcounter(genotypes(subnum))+1;
%     jobs{1}.stats{1}.factorial_design.des.fd.icell(genotypes(subnum)).levels = genotypes(subnum);
% 
% end
% jobs{1}.stats{1}.factorial_design.dir{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/COMT';
% jobs{1}.stats{2}.fmri_est.spmmat{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/COMT/SPM.mat';
% jobs{1}.stats{3}.con.spmmat{1} ='/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/COMT/SPM.mat';
% 
% rmdir(jobs{1}.stats{1}.factorial_design.dir{1},'s');
% mkdir(jobs{1}.stats{1}.factorial_design.dir{1});
%
% save(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/COMT.mat'], 'jobs', 'jobhelps');
% spm_jobman('run',['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/COMT.mat'])
% 
% clear jobs
% 
% COMTcounter = [1 1 1];
% load Group_analysis_COMT;
% for subnum = 1:length(subs)
%     subj = subs{subnum};
%     jobs{1}.stats{1}.factorial_design.des.fd.icell(genotypes(subnum)).scans{COMTcounter(genotypes(subnum)),1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/con_0001.img,1'];
%     COMTcounter(genotypes(subnum)) = COMTcounter(genotypes(subnum))+1;
%     jobs{1}.stats{1}.factorial_design.des.fd.icell(genotypes(subnum)).levels = genotypes(subnum);
% 
% end
% jobs{1}.stats{1}.factorial_design.dir{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/COMTXLoad';
% jobs{1}.stats{2}.fmri_est.spmmat{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/COMTXLoad/SPM.mat';
% jobs{1}.stats{3}.con.spmmat{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/COMT/SPM.mat';
% 
% rmdir(jobs{1}.stats{1}.factorial_design.dir{1},'s');
% mkdir(jobs{1}.stats{1}.factorial_design.dir{1});
%
% save(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/COMTXLoad.mat'], 'jobs', 'jobhelps');
% spm_jobman('run',['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/COMTXLoad.mat'])
% 
% clear jobs
% 
% DATcounter = [1 1];
% load Group_analysis_DAT;
% for subnum = 1:length(DATsubs)
%     subj = DATsubs{subnum};
%     jobs{1}.stats{1}.factorial_design.des.fd.icell(DATgenotype(subnum)).scans{DATcounter(DATgenotype(subnum)),1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/con_0003.img,1'];
%     DATcounter(DATgenotype(subnum)) = DATcounter(DATgenotype(subnum))+1;
%     jobs{1}.stats{1}.factorial_design.des.fd.icell(DATgenotype(subnum)).levels = DATgenotype(subnum);
% 
% end
% jobs{1}.stats{1}.factorial_design.dir{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/DAT';
% jobs{1}.stats{2}.fmri_est.spmmat{1} = [jobs{1}.stats{1}.factorial_design.dir{1} '/SPM.mat'];
% jobs{1}.stats{3}.con.spmmat{1} = jobs{1}.stats{2}.fmri_est.spmmat{1};
% 
% rmdir(jobs{1}.stats{1}.factorial_design.dir{1},'s');
% mkdir(jobs{1}.stats{1}.factorial_design.dir{1});
%
% save(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/DAT.mat'], 'jobs', 'jobhelps');
% spm_jobman('run',['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/DAT.mat'])
% 

% 
% clear jobs
% 
% genecounter = [1 1;1 1;1 1];
% load Group_analysis_COMTXDAT;
% for subnum = 1:length(DATsubs)
%     subj = DATsubs{subnum};
%     jobs{1}.stats{1}.factorial_design.des.fd.icell((COMTwDATgenotype(subnum)-1)*2 + DATgenotype(subnum)).scans{genecounter(COMTwDATgenotype(subnum),DATgenotype(subnum)),1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/con_0003.img,1'];
%     genecounter(COMTwDATgenotype(subnum),DATgenotype(subnum)) = genecounter(COMTwDATgenotype(subnum),DATgenotype(subnum))+1;
%     jobs{1}.stats{1}.factorial_design.des.fd.icell((COMTwDATgenotype(subnum)-1)*2 + DATgenotype(subnum)).levels = [COMTwDATgenotype(subnum) DATgenotype(subnum)];
% 
% end
% jobs{1}.stats{1}.factorial_design.dir{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/COMTXDAT';
% jobs{1}.stats{2}.fmri_est.spmmat{1} = [jobs{1}.stats{1}.factorial_design.dir{1} '/SPM.mat'];
% jobs{1}.stats{3}.con.spmmat{1} = jobs{1}.stats{2}.fmri_est.spmmat{1};
% 
% rmdir(jobs{1}.stats{1}.factorial_design.dir{1},'s');
% mkdir(jobs{1}.stats{1}.factorial_design.dir{1});
%
% save(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/COMTXDAT.mat'], 'jobs', 'jobhelps');
% spm_jobman('run',['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/COMTXDAT.mat'])
% 
% clear jobs
% 
% 
% 
% genecounter = [1 1;1 1;1 1];
% load Group_analysis_COMTXDAT;
% for subnum = 1:length(DATsubs)
%     subj = DATsubs{subnum};
%     jobs{1}.stats{1}.factorial_design.des.fd.icell((COMTwDATgenotype(subnum)-1)*2 + DATgenotype(subnum)).scans{genecounter(COMTwDATgenotype(subnum),DATgenotype(subnum)),1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/con_0001.img,1'];
%     genecounter(COMTwDATgenotype(subnum),DATgenotype(subnum)) = genecounter(COMTwDATgenotype(subnum),DATgenotype(subnum))+1;
%     jobs{1}.stats{1}.factorial_design.des.fd.icell((COMTwDATgenotype(subnum)-1)*2 + DATgenotype(subnum)).levels = [COMTwDATgenotype(subnum) DATgenotype(subnum)];
% 
% end
% jobs{1}.stats{1}.factorial_design.dir{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/COMTXDATXLoad';
% jobs{1}.stats{2}.fmri_est.spmmat{1} = [jobs{1}.stats{1}.factorial_design.dir{1} '/SPM.mat'];
% jobs{1}.stats{3}.con.spmmat{1} = jobs{1}.stats{2}.fmri_est.spmmat{1};
% 
% rmdir(jobs{1}.stats{1}.factorial_design.dir{1},'s');
% mkdir(jobs{1}.stats{1}.factorial_design.dir{1});
%
% save(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/COMTXDATXLoad.mat'], 'jobs', 'jobhelps');
% spm_jobman('run',['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/COMTXDATXLoad.mat'])
% 
% clear jobs


% DRD4counter = [1 1];
% load Group_analysis_DAT;
% for subnum = 1:length(DRD4subs)
%     subj = subs{subnum};
%     jobs{1}.stats{1}.factorial_design.des.fd.icell(DRD4genotype(subnum)).scans{DRD4counter(DRD4genotype(subnum)),1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/con_0003.img,1'];
%     DRD4counter(DRD4genotype(subnum)) = DRD4counter(DRD4genotype(subnum))+1;
%     jobs{1}.stats{1}.factorial_design.des.fd.icell(DRD4genotype(subnum)).levels = DRD4genotype(subnum);
% 
% end
% jobs{1}.stats{1}.factorial_design.des.fd.fact.name = 'DRD4';
% jobs{1}.stats{1}.factorial_design.dir{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/DRD4';
% jobs{1}.stats{2}.fmri_est.spmmat{1} = [jobs{1}.stats{1}.factorial_design.dir{1} '/SPM.mat'];
% jobs{1}.stats{3}.con.spmmat{1} = jobs{1}.stats{2}.fmri_est.spmmat{1};
% 
% rmdir(jobs{1}.stats{1}.factorial_design.dir{1},'s');
% mkdir(jobs{1}.stats{1}.factorial_design.dir{1});
% 
% save(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/DRD4.mat'], 'jobs', 'jobhelps');
% spm_jobman('run',['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/DRD4.mat'])
% 
% clear jobs
% 
% 
% DRD4counter = [1 1];
% load Group_analysis_DAT;
% for subnum = 1:length(DRD4subs)
%     subj = subs{subnum};
%     jobs{1}.stats{1}.factorial_design.des.fd.icell(DRD4genotype(subnum)).scans{DRD4counter(DRD4genotype(subnum)),1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/con_0001.img,1'];
%     DRD4counter(DRD4genotype(subnum)) = DRD4counter(DRD4genotype(subnum))+1;
%     jobs{1}.stats{1}.factorial_design.des.fd.icell(DRD4genotype(subnum)).levels = DRD4genotype(subnum);
% 
% end
% jobs{1}.stats{1}.factorial_design.des.fd.fact.name = 'DRD4';
% jobs{1}.stats{1}.factorial_design.dir{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/DRD4XLoad';
% jobs{1}.stats{2}.fmri_est.spmmat{1} = [jobs{1}.stats{1}.factorial_design.dir{1} '/SPM.mat'];
% jobs{1}.stats{3}.con.spmmat{1} = jobs{1}.stats{2}.fmri_est.spmmat{1};
% 
% rmdir(jobs{1}.stats{1}.factorial_design.dir{1},'s');
% mkdir(jobs{1}.stats{1}.factorial_design.dir{1});
% 
% save(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/DRD4XLoad.mat'], 'jobs', 'jobhelps');
% spm_jobman('run',['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/DRD4XLoad.mat'])
% 
% clear jobs
% 
% 
% BDNFcounter = [1 1 1];
% load Group_analysis_DAT;
% for subnum = 1:length(BDNFsubs)
%     subj = subs{subnum};
%     jobs{1}.stats{1}.factorial_design.des.fd.icell(BDNFgenotype(subnum)).scans{BDNFcounter(BDNFgenotype(subnum)),1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/con_0003.img,1'];
%     BDNFcounter(BDNFgenotype(subnum)) = BDNFcounter(BDNFgenotype(subnum))+1;
%     jobs{1}.stats{1}.factorial_design.des.fd.icell(BDNFgenotype(subnum)).levels = BDNFgenotype(subnum);
% 
% end
% jobs{1}.stats{1}.factorial_design.des.fd.fact.name = 'BDNF';
% jobs{1}.stats{1}.factorial_design.dir{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/BDNF';
% jobs{1}.stats{2}.fmri_est.spmmat{1} = [jobs{1}.stats{1}.factorial_design.dir{1} '/SPM.mat'];
% jobs{1}.stats{3}.con.spmmat{1} = jobs{1}.stats{2}.fmri_est.spmmat{1};
% 
% rmdir(jobs{1}.stats{1}.factorial_design.dir{1},'s');
% mkdir(jobs{1}.stats{1}.factorial_design.dir{1});
% 
% save(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/BDNF.mat'], 'jobs', 'jobhelps');
% spm_jobman('run',['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/BDNF.mat'])
% 
% clear jobs
% 
% 
% BDNFcounter = [1 1 1];
% load Group_analysis_DAT;
% for subnum = 1:length(BDNFsubs)
%     subj = subs{subnum};
%     jobs{1}.stats{1}.factorial_design.des.fd.icell(BDNFgenotype(subnum)).scans{BDNFcounter(BDNFgenotype(subnum)),1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/con_0001.img,1'];
%     BDNFcounter(BDNFgenotype(subnum)) = BDNFcounter(BDNFgenotype(subnum))+1;
%     jobs{1}.stats{1}.factorial_design.des.fd.icell(BDNFgenotype(subnum)).levels = BDNFgenotype(subnum);
% 
% end
% jobs{1}.stats{1}.factorial_design.des.fd.fact.name = 'BDNF';
% jobs{1}.stats{1}.factorial_design.dir{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/BDNFXLoad';
% jobs{1}.stats{2}.fmri_est.spmmat{1} = [jobs{1}.stats{1}.factorial_design.dir{1} '/SPM.mat'];
% jobs{1}.stats{3}.con.spmmat{1} = jobs{1}.stats{2}.fmri_est.spmmat{1};
% 
% rmdir(jobs{1}.stats{1}.factorial_design.dir{1},'s');
% mkdir(jobs{1}.stats{1}.factorial_design.dir{1});
% 
% save(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/BDNFXLoad.mat'], 'jobs', 'jobhelps');
% spm_jobman('run',['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/BDNFXLoad.mat'])
% 
% clear jobs
% 
% 
% SRTTcounter = [1 1 1];
% load Group_analysis_COMT;
% for subnum = 1:length(SRTTsubs)
%     subj = subs{subnum};
%     jobs{1}.stats{1}.factorial_design.des.fd.icell(SRTTgenotype(subnum)).scans{SRTTcounter(SRTTgenotype(subnum)),1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/con_0003.img,1'];
%     SRTTcounter(SRTTgenotype(subnum)) = SRTTcounter(SRTTgenotype(subnum))+1;
%     jobs{1}.stats{1}.factorial_design.des.fd.icell(SRTTgenotype(subnum)).levels = SRTTgenotype(subnum);
% 
% end
% jobs{1}.stats{1}.factorial_design.des.fd.fact.name = 'SRTT';
% jobs{1}.stats{1}.factorial_design.dir{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/SRTT';
% jobs{1}.stats{2}.fmri_est.spmmat{1} = [jobs{1}.stats{1}.factorial_design.dir{1} '/SPM.mat'];
% jobs{1}.stats{3}.con.spmmat{1} = jobs{1}.stats{2}.fmri_est.spmmat{1};
% 
% rmdir(jobs{1}.stats{1}.factorial_design.dir{1},'s');
% mkdir(jobs{1}.stats{1}.factorial_design.dir{1});
% 
% save(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/SRTT.mat'], 'jobs', 'jobhelps');
% spm_jobman('run',['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/SRTT.mat'])
% 
% clear jobs
% 
% 
% SRTTcounter = [1 1 1];
% load Group_analysis_COMT;
% for subnum = 1:length(SRTTsubs)
%     subj = subs{subnum};
%     jobs{1}.stats{1}.factorial_design.des.fd.icell(SRTTgenotype(subnum)).scans{SRTTcounter(SRTTgenotype(subnum)),1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/con_0001.img,1'];
%     SRTTcounter(SRTTgenotype(subnum)) = SRTTcounter(SRTTgenotype(subnum))+1;
%     jobs{1}.stats{1}.factorial_design.des.fd.icell(SRTTgenotype(subnum)).levels = SRTTgenotype(subnum);
% 
% end
% jobs{1}.stats{1}.factorial_design.des.fd.fact.name = 'SRTT';
% jobs{1}.stats{1}.factorial_design.dir{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/SRTTXLoad';
% jobs{1}.stats{2}.fmri_est.spmmat{1} = [jobs{1}.stats{1}.factorial_design.dir{1} '/SPM.mat'];
% jobs{1}.stats{3}.con.spmmat{1} = jobs{1}.stats{2}.fmri_est.spmmat{1};
% 
% rmdir(jobs{1}.stats{1}.factorial_design.dir{1},'s');
% mkdir(jobs{1}.stats{1}.factorial_design.dir{1});
% 
% save(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/SRTTXLoad.mat'], 'jobs', 'jobhelps');
% spm_jobman('run',['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/SRTTXLoad.mat'])
% 
% clear jobs

