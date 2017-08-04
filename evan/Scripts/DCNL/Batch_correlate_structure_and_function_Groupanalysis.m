warning off

subjects = {'112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','247','295','297','300','305','309','334','339','340','343','359','362','383','395','396','374','400','401','402','406','407','410','412','415','416','417','420'};




FCoutputimgname = '4mmspace_vmPFC_FC.nii';
SCoutputimgname_FA = '4mmspace_vmPFC_SC.nii';
%SCoutputimgname_Fiberdensity = 'Fiberdensity_to_Caudate_sphere.nii';

CorrelOutputfilename_FA = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/Func_Struct/correl_vmPFC_FC_FirstRestvsFA.nii';
%CorrelOutputfilename_Fiberdensity = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/Func_Struct/correl_Caudate_FC_FirstRestvsFiberdensity_sphere.nii';
TOutputfilename = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/Func_Struct/T_vmPFC_SCvsFC_FirstRest.nii';

disp('Loading Subjects')

for subject = 1:length(subjects)
    subjid = subjects{subject};
    
    %try
    tempdata = load_nii(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subjid '/Func_Struct/' FCoutputimgname]);
%     catch
%         gunzip(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subjid '/Func_Struct/' FCoutputimgname '.gz']);
%         tempdata = load_nii(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subjid '/Func_Struct/' FCoutputimgname]);
%     end
    
    fcdata(subject,:,:,:) = tempdata.img;
    
    tempdata = load_nii(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subjid '/Func_Struct/' SCoutputimgname_FA]);
    scdataFA(subject,:,:,:) = tempdata.img;
    
%     tempdata = load_nii(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subjid '/Func_Struct/' SCoutputimgname_Fiberdensity]);
%     scdataFiberdensity(subject,:,:,:) = tempdata.img;
    
    
end


reshapedfcdata = reshape(fcdata,[length(subjects) size(fcdata,2)*size(fcdata,3)*size(fcdata,4)]);



% disp(['Conducting correlation with ' SCoutputimgname_Fiberdensity])
% 
% CorrelOutput = tempdata;
% reshapedscdata = reshape(scdataFiberdensity,[length(subjects) size(scdataFiberdensity,2)*size(scdataFiberdensity,3)*size(scdataFiberdensity,4)]);
% 
% for vox = 1:size(reshapedfcdata,2)
%     if length(find(reshapedscdata(:,vox)>0)) > (length(subjects)/5);
%         substouse = find(reshapedscdata(:,vox)>0);
%         [rho, p] = corr(reshapedfcdata(substouse,vox),reshapedscdata(substouse,vox));
%         correlation(vox) = p;
%         direction(vox) = sign(rho);
%     else
%         correlation(vox) = 1;
%         direction(vox) = 0;
%     end
% end
% 
% correlationoutput = (1-reshape(correlation,[size(fcdata,2) size(fcdata,3) size(fcdata,4)]))  .* reshape(direction,[size(fcdata,2) size(fcdata,3) size(fcdata,4)]);
% CorrelOutput.img = correlationoutput;
% 
% save_nii(CorrelOutput,CorrelOutputfilename_Fiberdensity);
% 
% clear direction Correloutput correlationoutput correlation



disp(['Conducting correlation with ' SCoutputimgname_FA])


CorrelOutput = tempdata;


reshapedscdata = reshape(scdataFA,[length(subjects) size(scdataFA,2)*size(scdataFA,3)*size(scdataFA,4)]);


for vox = 1:size(reshapedfcdata,2)
    if length(find(reshapedscdata(:,vox)>0)) > (length(subjects)/5);
        substouse = find(reshapedscdata(:,vox)>0);
        [rho, p] = corr(reshapedfcdata(substouse,vox),reshapedscdata(substouse,vox));
        correlation(vox) = p;
        direction(vox) = sign(rho);
    else
        correlation(vox) = 1;
        direction(vox) = 0;
    end
end

correlationoutput = (1-reshape(correlation,[size(fcdata,2) size(fcdata,3) size(fcdata,4)]))  .* reshape(direction,[size(fcdata,2) size(fcdata,3) size(fcdata,4)]);
CorrelOutput.img = correlationoutput;

save_nii(CorrelOutput,CorrelOutputfilename_FA);

clear direction correlationoutput correlation





disp('Conducting T-tests')


TOutput = tempdata;

for vox = 1:size(reshapedfcdata,2)
    if length(find(reshapedscdata(:,vox)>0)) > 2  && length(find(reshapedscdata(:,vox)>0)) < (length(subjects)-2);
        
        [H, p, CI, stats] = ttest2(reshapedfcdata(find(reshapedscdata(:,vox)>0),vox), reshapedfcdata(find(reshapedscdata(:,vox)==0),vox));
        
        ttest(vox) = p;
        direction(vox) = sign(stats.tstat);
    else
        ttest(vox) = 1;
        direction(vox) = 0;
    end
end

ttestoutput = (1-reshape(ttest,[size(fcdata,2) size(fcdata,3) size(fcdata,4)])) .* reshape(direction,[size(fcdata,2) size(fcdata,3) size(fcdata,4)]);
TOutput.img = ttestoutput;

save_nii(TOutput,TOutputfilename);
