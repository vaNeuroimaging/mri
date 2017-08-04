warning off

subjects = {'101','102','113','118','120','122','125','127','132','138','147','150','151','154','156','159','160','161','162','166','172','181','182','187','202','207','211','214','215','221','225','229','232','233','242','250','254','255','272','274','112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','247','295','297','300','305','309','334','339','340','343','359','362','383','395','396','374'};


CorrelOutputfilename = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Nonstationarity/correl_with_connectivity_dACC_Nback_regress_vs_Rest.nii';


FCoutputimgname = 'dACC_Nback.nii';
NSoutputimgname = 'dACC_Nback_regress.nii';

for subject = 1:length(subjects)
    subjid = subjects{subject};
    
    tempdata = load_nii(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subjid '/Nonstationarity/RestVTask/' NSoutputimgname]);
    nsdatanback(subject,:,:,:) = tempdata.img;
    
    tempdata = load_nii(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subjid '/VoxelwiseConnectivity/RestVTask/' FCoutputimgname]);
    fcdatanback(subject,:,:,:) = tempdata.img;
    
end


FCoutputimgname = 'dACC_Rest.nii';
NSoutputimgname = 'dACC_Rest_regress.nii';

for subject = 1:length(subjects)
    subjid = subjects{subject};
    
    tempdata = load_nii(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subjid '/Nonstationarity/RestVTask/' NSoutputimgname]);
    nsdatarest(subject,:,:,:) = tempdata.img;
    
    tempdata = load_nii(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subjid '/VoxelwiseConnectivity/RestVTask/' FCoutputimgname]);
    fcdatarest(subject,:,:,:) = tempdata.img;
    
end

fcdata = fcdatanback - fcdatarest;
nsdata = nsdatanback - nsdatarest;


CorrelOutput = tempdata;


reshapedfcdata = reshape(fcdata,[length(subjects) size(fcdata,2)*size(fcdata,3)*size(fcdata,4)]);
reshapednsdata = reshape(nsdata,[length(subjects) size(nsdata,2)*size(nsdata,3)*size(nsdata,4)]);


for vox = 1:size(reshapedfcdata,2)
        [rho, p] = corr(reshapedfcdata(:,vox),reshapednsdata(:,vox));
        correlation(vox) = p;
end

correlationoutput = 1-reshape(correlation,[size(fcdata,2) size(fcdata,3) size(fcdata,4)]);
CorrelOutput.img = correlationoutput;

save_nii(CorrelOutput,CorrelOutputfilename);


