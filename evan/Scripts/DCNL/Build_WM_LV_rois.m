%Subjects
subjects = {'402'};
%'101','102','113','118','120','122','125','127','132','138','147','150','151','154','156','159','160','161','162','166','172','181','182','187','189','202','207','211','214','215','221','225','229','232','233','242','250','254','255','272','274','112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','283','327','281','292','301','307','322' '247','295','297','300','305','309','334','339','340','343','359','362','383','395','396','374','400','401','402','406','407','410','412','415','416','417','110','189','199','269'     
    
%Name of output ROIs.  ROI names will be SUBJ#_NAMEINTHISVARIABLE_roi.mat
filetypes = {'WM','CSF'};

%Header of normalized segmented WM and CSF images
filenames = {'wc2*','wc3*'};

%Location to output ROIs
outputpath = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ROI_analysis/Ind_ROIs/';



for subnum = 1:length(subjects)
    subj = subjects{subnum};
    disp(subj)
    
     %Location of structural data
    datapath = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subj '/Struct/'];
%     rawdata = dir([datapath '*MPRAGE*.nii']);
% 
%     
%     clear matlabbatch
%     load('/fmri/data3/Evan/Gene-Rest-Nback/Scripts/SPM8_Segment_Template.mat')
%     
%     matlabbatch{1}.spm.spatial.preproc.data{1} = [datapath rawdata(1).name ',1'];
%     
%     batchfilename = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subj '/Struct/Segment.mat'];
% 
%     save(batchfilename, 'matlabbatch');
% 
%     spm_jobman('run',batchfilename)

    for filetype = 1:length(filetypes)
        
        data = dir([datapath filenames{filetype} '.nii']);
        if isempty(data)
            data = dir([datapath filenames{filetype} '.img']);
        end

        outputname = [outputpath subj '_' filetypes{filetype} '_roi.mat'];

        o = [];
        d = [];

        imgname = [datapath data(1).name];

        [p f e] = fileparts(imgname);
        binf = 1;

        func = 'img >= .99';

        d = f; l = f;
        if ~isempty(func)
            d = [d ' func: ' func];
            l = [l '_f_' func];
        end
        if binf
            d = [d ' - binarized'];
            l = [l '_bin'];
        end
        o = maroi_image(struct('vol', spm_vol(imgname), 'binarize',binf,...
            'func', func));

        % convert to matrix format to avoid delicacies of image format
        o = maroi_matrix(o);


        o = descrip(o,d);
        o = label(o,l);

        varargout = {saveroi(o, outputname)};

    end

end