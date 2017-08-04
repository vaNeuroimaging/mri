cd /fmri/data3/Evan/Gene-Rest-Nback/Scripts

warning off

subs = {'112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','247','295','297','300','305','309','334','339','340','343','359','362','383','395','396','374','400','401','402','406','407','410','412','415','416','417','420'};
%'208','251','343',
%{'101','102','110','113','118','120','122','125','127','132','138','147','150','151','154','156','159','160','161','162','166','172','181','182','187','189','199','202','207','211','214','215','221','225','229','232','233','242','250','254','255','269','272','274','112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','247','295','297','300','305','309','334','339','340','343','359','362','383','395','396','374','400','401','402','406','407','410','412','415','416','417','420'};
%{'112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','247','295','297','300','305','309','334','339','340','343','359','362','383','395','396','374','400','401','402','406','407','410','412','415','416','417','420'};
%{'101','102','113','118','120','122','125','127','132','138','147','150','151','154','156','159','160','161','162','166','172','181','182','187','202','207','211','214','215','221','225','229','232','233','242','250','254','255','272','274','112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','247','295','297','300','305','309','334','339','340','343','359','362','383','395','396','374','400','401','402','406','407','410','412','415','416','417','420'};
%
% {'110','189','199','269'}
%
%{'208','222','227','253','256','270','292','301','305','309','343','362','374','383','395','415','416','417','420'}; %Chelsea's peeps '258', 

%
%
%'166','182'
%

%

% 
%'181','182',,
%'110','189','199','269',

genematrix = textread('/fmri/data3/Evan/Gene-Rest-Nback/Data/Genes.txt');

for sub = 1:length(subs);
    
    genotypes(sub,:) = genematrix(find(genematrix(:,1)==str2num(subs{sub})),:);
    
    %Subnum=1, DAT=2, COMT=3, DRD4=4, BDNF=5, SRTT=6, SRTTadj=7, DRD2=8, Gender=9, Menstrual=10 
    

end





behavname = {'RT','RT3B','Acc','Acc100%','Acc3B','ICV','ICV_3Back','Inattentive','H-I','BIS','STAI_Trait','Prescan_2B_Acc','Prescan_3B_Acc','Prescan_Overall_Acc'};

behavmatrix = textread('/fmri/data3/Evan/Gene-Rest-Nback/Analysis/Behavior/AllMeasures.txt');

for behavnum = 1:length(behavname)
    for sub = 1:length(subs)
        behavior{behavnum}(sub) = behavmatrix(find(behavmatrix(:,1)==str2num(subs{sub})),behavnum+1);
    end
end


Hipvolmatrix = textread('/fmri/data3/Evan/Gene-Rest-Nback/Analysis/Hip_Tracing/HipPctTCV.txt');

for sub = 1:length(subs)
    HipVol(sub) = Hipvolmatrix(find(Hipvolmatrix(:,1)==str2num(subs{sub})),2);
end


% 
% Howardbehavname = {'SCCT','ASRT'};
% 
% behavmatrix = textread('/fmri/data3/Evan/Gene-Rest-Nback/Analysis/Behavior/HowardMeasures.txt');
% 
% for behavnum = 1:length(Howardbehavname)
%     for sub = 1:length(subs)
%         Howardbehavior{behavnum}(sub) = behavmatrix(find(behavmatrix(:,1)==str2num(subs{sub})),behavnum+1);
%     end
% end

% 
% Connectname = {'DC-dmPFC','DC-lPar','DC-laIns','DC-ldlPFC','DC-raIns','DC-rdlPFC'};
% 
% matrix = textread('/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Seed-based/DAT_DC_FirstRest/Connectivity.txt');
% 
% for behavnum = 1:length(Connectname)
%     for sub = 1:length(subs)
%         Connectvals{behavnum}(sub) = matrix(find(matrix(:,1)==str2num(subs{sub})),behavnum+1);
%     end
% end







%-----------------------------------------------------------------------------------------------------
%%TWO-SAMPLE T-TEST

   seeds = {'VSi_bilat'};
   for seed = 1:length(seeds)
   seedname = seeds{seed};


load SPM8_2Sample_Ttest_template;
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1 = [];
matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2 = [];
for subnum = 1:length(subs)
    subj = subs{subnum};
    
    if genotypes(subnum,2) == 1
        matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1{length(matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1)+1,1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/VoxelwiseConnectivity/' seedname '_FirstRest_scrub.nii'];
        %['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/PPI/PPI_' seedname '_3B_vs_1B/con_0001.img'];
        %['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/VoxelwiseConnectivity/' seedname '_Nback.nii,1'];;
        %['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/TBSS/BDNF_92/stats/' subj '_skeletonized.nii,1'];
        
    elseif genotypes(subnum,2) == 2 %|| genotypes(subnum,5) == 3
        matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2{length(matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2)+1,1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/VoxelwiseConnectivity/' seedname '_FirstRest_scrub.nii'];
        %['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/PPI/PPI_' seedname '_3B_vs_1B/con_0001.img'];
        %['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/VoxelwiseConnectivity/' seedname '_Nback.nii,1'];;
        %['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/TBSS/BDNF_92/stats/' subj '_skeletonized.nii,1'];
        
        
    end
    
end

%matlabbatch{1}.spm.stats.factorial_design.masking.em{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/Wholebrain_EPI_p3_funcspace.nii';

matlabbatch{1}.spm.stats.factorial_design.dir{1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Seed-based/DAT_' seedname '_FirstRest_scrub/'];
%['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/PPI/DAT_' seedname '_3B_vs_1B_92/'];
%['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Seed-based/DAT_' seedname '_Nback_92/'];
%['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/TBSS/BDNF_92/'];
% 
%


try rmdir(matlabbatch{1}.spm.stats.factorial_design.dir{1},'s'); catch; end
mkdir(matlabbatch{1}.spm.stats.factorial_design.dir{1});

save([matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat'], 'matlabbatch');
spm_jobman('run',[matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat']);

clear matlabbatch
 end





% %-----------------------------------------------------------------------------------------------------
% %ONE-SAMPLE T-TEST
% 
% % seeds = {'DC_back','DC_down','DC_forward','DC_out','DC_up'};
% % %'DC','dACC','PCC','vmPFC',
% % for seed = 1:length(seeds)
% % seedname = seeds{seed};
% 
% load SPM8_1Sample_Ttest_template;
% subcounter = 0;
% for subnum = 1:length(subs)
%     subj = subs{subnum};
%     
%     %if (genotypes(subnum,2) == 1 || genotypes(subnum,2) == 2) %&& genotypes(subnum,7) ~= 0;
%         
%         subcounter = subcounter+1;
%         
%         matlabbatch{1}.spm.stats.factorial_design.des.t1.scans{subcounter,1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/SPM8/Cond/con_0001.img,1'];
%         
%         
%     %end
%     
% end
% 
% matlabbatch{1}.spm.stats.factorial_design.masking.em{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ICA/Old/ICA_lowthresh20dim_HBMpaper.gica/ROIs/DMN-TPN.nii';
% 
% matlabbatch{1}.spm.stats.factorial_design.dir{1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Activation/TaskvRest_50_DMN-TPN_alphasimmasktest/'];
% 
% try rmdir(matlabbatch{1}.spm.stats.factorial_design.dir{1},'s'); catch; end
% mkdir(matlabbatch{1}.spm.stats.factorial_design.dir{1});
% 
% save([matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat'], 'matlabbatch');
% spm_jobman('run',[matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat']);
% 
% clear matlabbatch
% 
% %end






% %-----------------------------------------------------------------------------------------------------
% %SIMPLE CORRELATION
% 
% 
% 
% seeds = {'laIFG','lpIFG'};
% %{'dmPFC','raIns','laIns','lPFC','SMA'};
% %'DC','dACC','PCC','vmPFC',
% for seed = 1:length(seeds)
% seedname = seeds{seed};
% 
% %for behavnum = [8 9 10 14];
% 
% load SPM8_1Sample_Ttest_withCov_template;
% 
% subcounter = 0;
% for subnum = 1:length(subs)
%     subj = subs{subnum};
%     roi = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Seed-based/DAT_VSi_bilat_FirstRest_scrub/', seedname , '_roi.mat'];
%     rois = maroi('load_cell', roi);
%     Y = getdata(rois{1},['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/VoxelwiseConnectivity/VSi_bilat_FirstRest_scrub.nii'],'l');
%     mean_within_roi = mean(Y(find(~isnan(Y))));
%     %if behavior{behavnum}(subnum) ~= 9999
%        subcounter = subcounter+1;
%         
%        matlabbatch{1}.spm.stats.factorial_design.des.t1.scans{subcounter} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/VoxelwiseConnectivity/VSiconnect/' seedname '_FirstRest_scrub_VSiregress.nii'];
%        %['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/SPM8/Cond/con_0001.img,1'];
%        %['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/VoxelwiseConnectivity/PCC_aal_FirstRest.nii'];
%        
%        %matlabbatch{1}.spm.stats.factorial_design.cov.c(subcounter) = behavior{behavnum}(subnum);
%        matlabbatch{1}.spm.stats.factorial_design.cov.c(subcounter) = mean_within_roi;
%     %end
% end
% matlabbatch{1}.spm.stats.factorial_design.dir{1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Seed-based/VSi_connect_seeds/VSiC_connectval_to_' seedname '_FirstRest_scrub_partialVSisignal/'];
% matlabbatch{1}.spm.stats.factorial_design.masking.em{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/Wholebrain_limited.img';
% 
% try rmdir(matlabbatch{1}.spm.stats.factorial_design.dir{1},'s'); catch; end
% mkdir(matlabbatch{1}.spm.stats.factorial_design.dir{1});
% 
% save([matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat'], 'matlabbatch');
% spm_jobman('run',[matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat']);
% 
% clear matlabbatch
% 
% %end
% end






% %-----------------------------------------------------------------------------------------------------
% %%3X2 ANOVA
% 
% %COMTXDAT
% load('SPM8_COMTxDAT_template');
% matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).name = 'BDNF';
% 
% for i = 1:6; matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(i).scans = []; end
% 
% 
% for sub = 1:length(subs)
%     subj = subs{sub};
%     
%     if (genotypes(sub,3) == 1 || genotypes(sub,3) == 2 || genotypes(sub,3) == 3 ) && (genotypes(sub,5) == 1 || genotypes(sub,5) == 2 || genotypes(sub,5) == 3)
%         
%         cellnum = (genotypes(sub,3)-1)*2 + (genotypes(sub,5)>1) +1;
%     
%         matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(cellnum).levels = [genotypes(sub,3) ((genotypes(sub,5)>1) +1)];
%         matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(cellnum).scans{length(matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(cellnum).scans)+1,1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/VoxelwiseConnectivity/DC_bilat_FirstRest.nii,1'];
%     
%     end
%     
% end
% matlabbatch{1}.spm.stats.factorial_design.dir{1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Seed-based/COMTxBDNF_DC_bilat_FirstRest/'];
% try rmdir(matlabbatch{1}.spm.stats.factorial_design.dir{1},'s');catch;end
% mkdir(matlabbatch{1}.spm.stats.factorial_design.dir{1});
% save([matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat'], 'matlabbatch');
% spm_jobman('run',[matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat']);
% clear matlabbatch
% 




% %-----------------------------------------------------------------------------------------------------
% %3(btwn) X 2(btwn) X 3(within) ANOVA with covariate
% 
% %COMTXDATXLoad covary Gender
% 
% 
% load('SPM8_COMTxDATxGender_covariate_template');
% 
% matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(3).name = 'Load';
% matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(3).dept = 1;
% matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(3).levels = 3;
% 
% 
% matlabbatch{1}.spm.stats.factorial_design.cov.cname = 'Gender';
% 
%  for i = 1:18; matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(i).scans = []; covariate{i} = []; end
% 
% for sub = 1:length(subs)
%     subj = subs{sub};
%     
%     if (genotypes(sub,3) == 1 || genotypes(sub,3) == 2 || genotypes(sub,3) == 3 ) && (genotypes(sub,2) == 1 || genotypes(sub,2) == 2)
%         
%         cellnum = (genotypes(sub,3)-1)*6 + (genotypes(sub,2)-1)*3 + 1;
%         
%         matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(cellnum).levels = [genotypes(sub,3) genotypes(sub,2) 1];
%         matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(cellnum).scans{length(matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(cellnum).scans)+1,1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/SPM8/Cond/con_0002.img,1'];
%         covariate{cellnum}(size(covariate{cellnum})+1,1) = genotypes(sub,9);
%         
%         
%         cellnum = (genotypes(sub,3)-1)*6 + (genotypes(sub,2)-1)*3 + 2;
%         
%         matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(cellnum).levels = [genotypes(sub,3) genotypes(sub,2) 2];
%         matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(cellnum).scans{length(matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(cellnum).scans)+1,1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/SPM8/Cond/con_0003.img,1'];
%         covariate{cellnum}(size(covariate{cellnum})+1,1) = genotypes(sub,9);
%         
%         
%         cellnum = (genotypes(sub,3)-1)*6 + (genotypes(sub,2)-1)*3 + 3;
%         
%         matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(cellnum).levels = [genotypes(sub,3) genotypes(sub,2) 3];
%         matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(cellnum).scans{length(matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(cellnum).scans)+1,1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/SPM8/Cond/con_0004.img,1'];
%         covariate{cellnum}(size(covariate{cellnum})+1,1) = genotypes(sub,9);
%         
%     end
%     
% end
% 
% finalcovariate = []; for i = 1:length(covariate); finalcovariate = [finalcovariate; covariate{i}]; end
% matlabbatch{1}.spm.stats.factorial_design.cov.c = finalcovariate;
% matlabbatch{1}.spm.stats.factorial_design.dir{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Activation/COMTXDATXLoad_covaryGender_92/';
% %matlabbatch{1}.spm.stats.factorial_design.masking.em{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/aal_allregions.nii';
% try rmdir(matlabbatch{1}.spm.stats.factorial_design.dir{1},'s');catch;end
% mkdir(matlabbatch{1}.spm.stats.factorial_design.dir{1});
% save([matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat'], 'matlabbatch');
% spm_jobman('run',[matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat']);
% clear matlabbatch






% %-----------------------------------------------------------------------------------------------------
% %%2X2 ANOVA
% 
% %COMTXDAT
% load('SPM8_COMTxDAT_template');
% matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).name = 'BDNF';
% matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).levels = 2;
% matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).name = 'Gender';
% 
% matlabbatch{1}.spm.stats.factorial_design.des.fd.icell = matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(1:4);
% 
% for i = 1:4; matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(i).scans = []; end
% 
% 
% for sub = 1:length(subs)
%     subj = subs{sub};
%     
%     if (genotypes(sub,5) == 1 || genotypes(sub,5) == 2 || genotypes(sub,5) == 3 )
%         
%         cellnum = (genotypes(sub,5)>1)*2 + genotypes(sub,9);
%     
%         matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(cellnum).levels = [(genotypes(sub,5)>1)+1 genotypes(sub,9)];
%         matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(cellnum).scans{length(matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(cellnum).scans)+1,1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/VoxelwiseConnectivity/aal_Hipp_both_FirstRest.nii,1'];
%     
%     end
%     
% end
% matlabbatch{1}.spm.stats.factorial_design.dir{1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Seed-based/BDNFxGender_Hipp_FirstRest/'];
% try rmdir(matlabbatch{1}.spm.stats.factorial_design.dir{1},'s');catch;end
% mkdir(matlabbatch{1}.spm.stats.factorial_design.dir{1});
% save([matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat'], 'matlabbatch');
% spm_jobman('run',[matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat']);
% clear matlabbatch



% %-----------------------------------------------------------------------------------------------------
% %COMT(linear)XDAT(linear)
% 
% 
% load('SPM8_COMTxDAT_linear_template');
% 
% subcount = 0;
% for sub = 1:length(subs)
%     subj = subs{sub};
%     
%     
%     if (genotypes(sub,3) == 1 || genotypes(sub,3) == 2 || genotypes(sub,3) == 3 ) && (genotypes(sub,2) == 1 || genotypes(sub,2) == 2) && ~strcmp(subj,'166')
%         
%         subcount = subcount+1;
%         
%         matlabbatch{1}.spm.stats.factorial_design.des.t1.scans{subcount,1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/TBSS/DAT_91/stats/' subj '_skeletonized.nii,1'];
%         
%         DATregressor(subcount,1) = 2*genotypes(sub,2)-3;
%         COMTregressor(subcount,1) = genotypes(sub,3);
%         DATxCOMTregressor(subcount,1) = (2*genotypes(sub,2)-3)*genotypes(sub,3);
%         
%     end
%     
%     matlabbatch{1}.spm.stats.factorial_design.cov(1).c = [DATregressor(:,1)];
%     matlabbatch{1}.spm.stats.factorial_design.cov(2).c = [COMTregressor(:,1)];
%     matlabbatch{1}.spm.stats.factorial_design.cov(3).c = [DATxCOMTregressor(:,1)];
%     
% end
% matlabbatch{1}.spm.stats.factorial_design.dir{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/TBSS/COMTXDAT_linear_91/';
% %matlabbatch{1}.spm.stats.factorial_design.masking.em{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/aal_allregions.nii';
% try rmdir(matlabbatch{1}.spm.stats.factorial_design.dir{1},'s');catch;end
% mkdir(matlabbatch{1}.spm.stats.factorial_design.dir{1});
% save([matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat'], 'matlabbatch');
% spm_jobman('run',[matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat']);
% clear matlabbatch




% %-----------------------------------------------------------------------------------------------------
% % 3(within) x linear x linear ANOVA
% 
% %COMT(linear)XDAT(linear)XLoad
% 
% 
% load('SPM8_COMTxDAT_linear_xLoad_template');
% 
% subcount = 0;
% for sub = 1:length(subs)
%     subj = subs{sub};
%     
%     
%     if (genotypes(sub,3) == 1 || genotypes(sub,3) == 2 || genotypes(sub,3) == 3 ) && (genotypes(sub,2) == 1 || genotypes(sub,2) == 2)
%         
%         subcount = subcount+1;
%         
%         matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(1).levels = 1;
%         matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(1).scans{subcount,1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/SPM8/Cond/con_0002.img,1'];
%         
%         matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(2).levels = 2;
%         matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(2).scans{subcount,1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/SPM8/Cond/con_0003.img,1'];
%         
%         matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(3).levels = 3;
%         matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(3).scans{subcount,1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/SPM8/Cond/con_0004.img,1'];
%         
%         DATregressor(subcount,1) = 2*genotypes(sub,2)-3;
%         COMTregressor(subcount,1) = genotypes(sub,3);
%         DATxCOMTregressor(subcount,1) = (2*genotypes(sub,2)-3)*genotypes(sub,3);
%         
%     end
%     
%     matlabbatch{1}.spm.stats.factorial_design.cov(1).c = [DATregressor(:,1);DATregressor(:,1);DATregressor(:,1)];
%     matlabbatch{1}.spm.stats.factorial_design.cov(2).c = [COMTregressor(:,1);COMTregressor(:,1);COMTregressor(:,1)];
%     matlabbatch{1}.spm.stats.factorial_design.cov(3).c = [DATxCOMTregressor(:,1);DATxCOMTregressor(:,1);DATxCOMTregressor(:,1)];
%     
% end
% matlabbatch{1}.spm.stats.factorial_design.dir{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Activation/COMTXDAT_linear_XLoad_92/';
% %matlabbatch{1}.spm.stats.factorial_design.masking.em{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/aal_allregions.nii';
% try rmdir(matlabbatch{1}.spm.stats.factorial_design.dir{1},'s');catch;end
% mkdir(matlabbatch{1}.spm.stats.factorial_design.dir{1});
% save([matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat'], 'matlabbatch');
% spm_jobman('run',[matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat']);
% clear matlabbatch
% 









% %-----------------------------------------------------------------------------------------------------
% %ONE-WAY ANOVA
% 
% %  seeds = {'DC','VSi','DCP','VRP'};
% %  for seed = 1:length(seeds)
% %  seedname = seeds{seed};
% 
% 
% load SPM8_Oneway_template;
% 
% matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(1).scans = [];
% matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(2).scans = [];
% matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(3).scans = [];
% 
% for subnum = 1:length(subs)
%     subj = subs{subnum};
%     
%     if (genotypes(subnum,7) == 1 || genotypes(subnum,7) == 2 || genotypes(subnum,7) == 3)
%         numprevscans = size(matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(genotypes(subnum,7)).scans,1);
%         matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(genotypes(subnum,7)).scans{numprevscans+1,1} = ...
%             ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/VoxelwiseConnectivity/aal_Amygdala_bilat_FirstRest.nii'];
%         
%     end
%     
%     
% end
% 
% matlabbatch{1}.spm.stats.factorial_design.masking.em{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/Wholebrain_limited.img';
% 
% matlabbatch{1}.spm.stats.factorial_design.dir{1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Seed-based/5HTTadj_aal_Amygdala_bilat_FirstRest/'];
% 
% try rmdir(matlabbatch{1}.spm.stats.factorial_design.dir{1},'s'); catch; end
% mkdir(matlabbatch{1}.spm.stats.factorial_design.dir{1});
% 
% save([matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat'], 'matlabbatch');
% spm_jobman('run',[matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat']);
% 
% clear matlabbatch
% 
% % end











% %-----------------------------------------------------------------------------------------------------
% %TWO SAMPLE T-TEST WITH COVARIATE
% 
% for behavnum = [2];
% 
% load SPM8_2Sample_Ttest_withCov_template;
% matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1 = [];
% matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2 = [];
% 
% for subnum = 1:length(subs)
%     subj = subs{subnum};
%     
%     
%     
%     if genotypes(subnum,2) == 1
%         matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1{length(ma
%         tlabbatch{1}.spm.stats.factorial_design.des.t2.scans1)+1,1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/VoxelwiseConnectivity/DC_bilat_FirstRest.nii'];
%         
%         cov_gr1(length(matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1)) = Howardbehavior{behavnum}(subnum);
%         
%     elseif genotypes(subnum,2) == 2
%         matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2{length(matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2)+1,1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/VoxelwiseConnectivity/DC_bilat_FirstRest.nii'];
%         
%         cov_gr2(length(matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2)) = Howardbehavior{behavnum}(subnum);
%     end
%         
%    
%     
% end
% matlabbatch{1}.spm.stats.factorial_design.dir{1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Seed-based/DC_FirstRest_vs_' Howardbehavname{behavnum} '_byDAT/'];
% 
% matlabbatch{1}.spm.stats.factorial_design.cov.c = [cov_gr1 cov_gr2];
% matlabbatch{1}.spm.stats.factorial_design.cov.cname = Howardbehavname{behavnum};
% 
% try rmdir(matlabbatch{1}.spm.stats.factorial_design.dir{1},'s'); catch; end
% mkdir(matlabbatch{1}.spm.stats.factorial_design.dir{1});
% 
% save([matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat'], 'matlabbatch');
% spm_jobman('run',[matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat']);
% 
% clear matlabbatch
% 
% end










% %-----------------------------------------------------------------------------------------------------
% %2X2 ANOVA
% 
% load('SPM8_COMTxDAT_template');
% 
%  matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).levels = 2;
% % matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).name = 'DRD2';
%  matlabbatch{1}.spm.stats.factorial_design.des.fd.icell = matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(1:4);
% for i = 1:4; matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(i).scans = []; end
% 
% for sub = 1:length(subs)
%     subj = subs{sub};
%     
%     if (genotypes(sub,2) == 1 || genotypes(sub,2) == 2 ) && (genotypes(sub,3) == 2 || genotypes(sub,3) == 3)
%         
%         cellnum = (genotypes(sub,3)-2)*2 + genotypes(sub,2);
%         
%         matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(cellnum).levels = [genotypes(sub,3)-1 genotypes(sub,2)];
%         matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(cellnum).s
%         cans{length(matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(cellnum).scans)+1,1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/SPM8/Cond/con_0005.img,1'];
%         
%         
%     end
% end
% matlabbatch{1}.spm.stats.factorial_design.dir{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Activation/COMTXDAT_2Bv1B_92_noVV/';
% try;rmdir(matlabbatch{1}.spm.stats.factorial_design.dir{1},'s');catch;end
% mkdir(matlabbatch{1}.spm.stats.factorial_design.dir{1});
% save([matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat'], 'matlabbatch');
% spm_jobman('run',[matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat']);
% clear matlabbatch






% %-----------------------------------------------------------------------------------------------------
% %3(btwn) X 2(btwn) X 3(within) ANOVA
% 
% %COMTXDATXLoad
% 
% 
% load('SPM8_COMTxDATxGender_template');
% 
% matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(3).name = 'Load';
% matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(3).dept = 1;
% matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(3).levels = 3;
% matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).levels = 2;
% 
%  for i = 1:12; matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(i).scans = []; end
% 
% for sub = 1:length(subs)
%     subj = subs{sub};
%     
%     if (genotypes(sub,3) == 2 || genotypes(sub,3) == 3 ) && (genotypes(sub,2) == 1 || genotypes(sub,2) == 2)
%         %genotypes(sub,3) == 1 || 
%         cellnum = (genotypes(sub,3)-2)*6 + (genotypes(sub,2)-1)*3 + 1;
%         
%         matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(cellnum).levels = [genotypes(sub,3)-1 genotypes(sub,2) 1];
%         matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(cellnum).scans{length(matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(cellnum).scans)+1,1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/SPM8/Cond/con_0002.img,1'];
%         
%         
%         cellnum = (genotypes(sub,3)-2)*6 + (genotypes(sub,2)-1)*3 + 2;
%         
%         matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(cellnum).levels = [genotypes(sub,3)-1 genotypes(sub,2) 2];
%         matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(cellnum).scans{length(matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(cellnum).scans)+1,1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/SPM8/Cond/con_0003.img,1'];
%         
%         
%         cellnum = (genotypes(sub,3)-2)*6 + (genotypes(sub,2)-1)*3 + 3;
%         
%         matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(cellnum).levels = [genotypes(sub,3)-1 genotypes(sub,2) 3];
%         matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(cellnum).scans{length(matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(cellnum).scans)+1,1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/SPM8/Cond/con_0004.img,1'];
%             
%     end
%     
% end
% matlabbatch{1}.spm.stats.factorial_design.dir{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Activation/COMTXDATXLoad_92_noVV/';
% %matlabbatch{1}.spm.stats.factorial_design.masking.em{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/aal_allregions.nii';
% try rmdir(matlabbatch{1}.spm.stats.factorial_design.dir{1},'s');catch;end
% mkdir(matlabbatch{1}.spm.stats.factorial_design.dir{1});
% save([matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat'], 'matlabbatch');
% spm_jobman('run',[matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat']);
% clear matlabbatch












% %-----------------------------------------------------------------------------------------------------
% %SIMPLE CORRELATION WITH OTHER REGION'S ACTIVATION
% 
% 
% 
% seeds = {'PCC'};
% for seed = 1:length(seeds)
% seedname = seeds{seed};
% 
% load SPM8_1Sample_Ttest_withCov_template;
% 
% roiname = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Seed-based/DAT_DC_FirstRest/laIns_roi.mat';
% rois = maroi('load_cell', roiname);
% 
% subcounter = 0;
% for subnum = 1:length(subs)
%     subj = subs{subnum};
%        subcounter = subcounter+1;
%         
%        matlabbatch{1}.spm.stats.factorial_design.des.t1.scans{subcounter} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/VoxelwiseConnectivity/' seedname '_DeLuca_Nback.nii,1'];
%        
%        Y = getdata(rois{1},['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/VoxelwiseConnectivity/DC_bilat_FirstRest.nii'],'l');
%        
%        matlabbatch{1}.spm.stats.factorial_design.cov.c(subcounter) =  mean(Y(find(~isnan(Y))));
%        
% end
% matlabbatch{1}.spm.stats.factorial_design.dir{1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Seed-based/' seedname '_Nback_vs_DC-laInsConnect/'];
% matlabbatch{1}.spm.stats.factorial_design.masking.em{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/Wholebrain_limited.img';
% 
% try rmdir(matlabbatch{1}.spm.stats.factorial_design.dir{1},'s'); catch; end
% mkdir(matlabbatch{1}.spm.stats.factorial_design.dir{1});
% 
% save([matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat'], 'matlabbatch');
% spm_jobman('run',[matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat']);
% 
% clear matlabbatch
% 
% end







% %-----------------------------------------------------------------------------------------------------
% %%3X2X2 ANOVA
% 
% load('SPM8_COMTxDATxGender_template');
% 
% % matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).name = 'DRD2';
% % matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).levels = 2;
% % matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(3).name = 'DRD4';
% % matlabbatch{1}.spm.stats.factorial_design.des.fd.icell = matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(1:8);
% 
% for i = 1:12; matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(i).scans = []; end
% 
% 
% for sub = 1:length(subs)
%     subj = subs{sub};
%     
%     if (genotypes(sub,2) == 1 || genotypes(sub,2) == 2) %&& (genotypes(sub,4) == 1 || genotypes(sub,4) == 2)
%         
%         cellnum = (genotypes(sub,3)-1)*4 + (genotypes(sub,2)-1)*2 + genotypes(sub,8);
%     
%         matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(cellnum).levels = [genotypes(sub,3) genotypes(sub,2) genotypes(sub,8)];
%         matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(cellnum).scans{length(matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(cellnum).scans)+1,1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/SPM8/Cond/con_0009.img,1'];
%     
%     end
%     
% end
% matlabbatch{1}.spm.stats.factorial_design.dir{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Activation/COMTxDATxGender_2B3B>1B/';
% %matlabbatch{1}.spm.stats.factorial_design.masking.em{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/aal_allregions.nii';
% try rmdir(matlabbatch{1}.spm.stats.factorial_design.dir{1},'s');catch;end
% mkdir(matlabbatch{1}.spm.stats.factorial_design.dir{1});
% save([matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat'], 'matlabbatch');
% spm_jobman('run',[matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat']);
% clear matlabbatch















% %-----------------------------------------------------------------------------------------------------
% %PAIRED T-TEST
% 
% load SPM8_Paired_Ttest_template;
% subcounter = 0;
% for subnum = 1:length(subs)
%     subj = subs{subnum};
%     
%     %if (genotypes(sub,2) == 1 || genotypes(sub,2) == 2)
%     
%     subcounter = subcounter+1;
%     
%     matlabbatch{1}.spm.stats.factorial_design.des.pt.pair(subcounter) = matlabbatch{1}.spm.stats.factorial_design.des.pt.pair(1);
%     
%     matlabbatch{1}.spm.stats.factorial_design.des.pt.pair(subcounter).scans{1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/VoxelwiseConnectivity/DCP_bilat_FirstRest.nii,1'];
%     matlabbatch{1}.spm.stats.factorial_design.des.pt.pair(subcounter).scans{2} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/VoxelwiseConnectivity/DCP_bilat_Rest.nii,1'];
%     
%     %end
% 
% end
% matlabbatch{3}.spm.stats.con.consess{1}.tcon.name = 'First > Rest';
% matlabbatch{3}.spm.stats.con.consess{1}.tcon.convec = [zeros(1,subcounter) 1 -1];
% matlabbatch{3}.spm.stats.con.consess{2}.tcon.name = 'Rest > First';
% matlabbatch{3}.spm.stats.con.consess{2}.tcon.convec = [zeros(1,subcounter) -1 1];
% 
% matlabbatch{1}.spm.stats.factorial_design.masking.em{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/aal_allregions.nii';
% 
% matlabbatch{1}.spm.stats.factorial_design.dir{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Seed-based/DCP_FirstRestVsRest/';
% 
% try rmdir(matlabbatch{1}.spm.stats.factorial_design.dir{1},'s'); catch; end
% mkdir(matlabbatch{1}.spm.stats.factorial_design.dir{1});
% 
% save([matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat'], 'matlabbatch');
% spm_jobman('run',[matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat']);
% 
% clear matlabbatch















%-----------------------------------------------------------------------------------------------------
%%ONE-WAY REPEATED MEASURES ANOVA
% 
% load SPM8_Oneway_Load_template;
% for subnum = 1:length(subs)
%     subj = subs{subnum};
%     
%     matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(1).scans{subnum,1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/Nonstationarity/RestVTask/dACC_Nback_1Back.nii,1'];
%     matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(2).scans{subnum,1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/Nonstationarity/RestVTask/dACC_Nback_2Back.nii,1'];
%     matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(3).scans{subnum,1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/Nonstationarity/RestVTask/dACC_Nback_3Back.nii,1'];
% 
% end
% 
% matlabbatch{1}.spm.stats.factorial_design.masking.em{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/Wholebrain_limited.img';
% 
% matlabbatch{1}.spm.stats.factorial_design.dir{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Nonstationarity/dACC_Nback_Load/';
% 
% try rmdir(matlabbatch{1}.spm.stats.factorial_design.dir{1},'s'); catch; end
% mkdir(matlabbatch{1}.spm.stats.factorial_design.dir{1});
% 
% save([matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat'], 'matlabbatch');
% spm_jobman('run',[matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat']);
% 
% clear matlabbatch







%-----------------------------------------------------------------------------------------------------
%2X2 WITHIN-BETWEEN ANOVA
% 
% 
% seedname = 'aDMN-vmPFC';
% 
% load SPM8_LevelXDAT_template
% for i = 1:4;
% matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(i).scans = []; end
% 
% matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(1).name = 'Run';
% matlabbatch{1}.spm.stats.factorial_design.des.fd.fact(2).name = 'DAT';
% 
% 
% for subnum = 1:length(subs)
%     subj = subs{subnum};
%     
%     if genotypes(subnum,2) == 1;
%     
%     matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(1).levels = [1;1];
%     matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(1).scans{length(matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(1).scans)+1,1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/VoxelwiseConnectivity/' seedname '_Rest.nii,1'];
%     matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(3).levels = [2;1];
%     matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(3).scans{length(matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(3).scans)+1,1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/VoxelwiseConnectivity/' seedname '_Nback.nii,1'];  
%     
%     elseif genotypes(subnum,2) == 2
%         
%     matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(2).levels = [1;2];
%     matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(2).scans{length(matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(2).scans)+1,1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/VoxelwiseConnectivity/' seedname '_Rest.nii,1'];
%     matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(4).levels = [2;2];
%     matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(4).scans{length(matlabbatch{1}.spm.stats.factorial_design.des.fd.icell(4).scans)+1,1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/VoxelwiseConnectivity/' seedname '_Nback.nii,1'];    
%     
%     end
%     
% end
% matlabbatch{1}.spm.stats.factorial_design.dir{1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Seed-based/' seedname '_StateXDAT/'];
% 
% matlabbatch{1}.spm.stats.factorial_design.masking.em{1} = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/aal_allregions.nii';
% 
% try;rmdir(matlabbatch{1}.spm.stats.factorial_design.dir{1},'s');catch;end
% mkdir(matlabbatch{1}.spm.stats.factorial_design.dir{1});
% 
% save([matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat'], 'matlabbatch');
% spm_jobman('run',[matlabbatch{1}.spm.stats.factorial_design.dir{1}(1:end-1) '.mat']);
% 
% clear matlabbatch










%-----------------------------------------------------------------------------------------------------
%%TWO SAMPLE T-TEST (BDNF)

%%DTI
% load SPM8_2Sample_Ttest_template;
% matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1 = [];
% matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2 = [];
% for subnum = 1:length(subs)
%     subj = subs{subnum};
%     
%     if ~strcmp(subj,'166') && (genotypes(sub,5) == 1 || genotypes(sub,5) == 2 || genotypes(sub,5) == 3)
%         
%         if genotypes(sub,5) == 1
%             matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1{length(matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1)+1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/TBSS/BDNF_77/stats/' subj '_skeletonized.nii,1'];
%             
%         else
%             matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2{length(matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2)+1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/TBSS/BDNF_77/stats/' subj '_skeletonized.nii,1'];
%             
%         end
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


%-----------------------------------------------------------------------------------------------------
%%TWO SAMPLE T-TEST WITH COVARIATE (VBM)
% 
%%VBM group analysis
% BDNFcounter = [1 1];
% load SPM8_2Sample_Ttest_withCov_template;
% matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1 = [];
% matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2 = [];
% for subnum = 1:length(subs)
%     subj = subs{subnum};
%     
%     if (genotypes(sub,5) == 1 || genotypes(sub,5) == 2 || genotypes(sub,5) == 3)
%     
%     if genotypes(sub,5) == 1
%         matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1{length(matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1)+1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subj '/Struct/smwrc1Registered_MPRAGE_thr3.nii,1'];
%         
%         t = VBM_get_TCV(['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subj '/Struct/smwrc1Registered_MPRAGE_thr3.nii,1'; '/fmri/data3/Evan/Gene-Rest-Nback/Data/' subj '/Struct/smwrc1Registered_MPRAGE_thr3.nii,1']);
%         
%         TCV_gr1(length(matlabbatch{1}.spm.stats.factorial_design.des.t2.scans1),1) = sum(t);
%         
%     else
%         matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2{length(matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2)+1} = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subj '/Struct/smwrc1Registered_MPRAGE_thr3.nii,1'];
%         
%         t = VBM_get_TCV(['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subj '/Struct/smwrc1Registered_MPRAGE_thr3.nii,1'; '/fmri/data3/Evan/Gene-Rest-Nback/Data/' subj '/Struct/smwrc1Registered_MPRAGE_thr3.nii,1']);
%         
%         TCV_gr2(length(matlabbatch{1}.spm.stats.factorial_design.des.t2.scans2)) = sum(t);
%     end
%         
%     end
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












