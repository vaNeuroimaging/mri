%Interrogate_ROI
%
%Extracts the activity within defined ROIs for the subjects.  Saves
%those activations in in one big file that can be easily put into Excel.
%
%ROIs and subjects are specified at the top of the script
%
%
%Created by E. Gordon 5/08/08


warning off

directory = pwd;



subjects = {'101','102','110','113','118','120','122','125','127','132','138','147','150','151','154','156','159','160','161','162','166','172','187','189','199','202','207','211','214','215','221','225','229','232','233','242','250','254','255','269','272','274','112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','247','295','297','300','305','309','334','339','340','343','359','362','383','395','396','374','400','401','402','406','407','410','412','415','416','417','420'};
%{'112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','247','295','297','300','305','309','334','339','340','343','359','362','383','395','396','374','400','401','402','406','407','410','412','415','416','417','420'};
%
%{'101','102','110','113','118','120','122','125','127','132','138','147','150','151','154','156','159','160','161','162','166','172','187','189','199','202','207','211','214','215','221','225','229','232','233','242','250','254','255','269','272','274','112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','247','295','297','300','305','309','334','339','340','343','359','362','383','395','396','374','400','401','402','406','407','410','412','415','416','417','420'};
%'101' '102'  '113' '118' '120' '122' '125' '127' '132' '138' '147' '150' '151' '154' '156' '159' '160' '161' '162' '166' '172' '181' '182' '187' '202' '207' '211' '214' '215' '221' '225' '229' '232' '233' '242' '250' '254' '255' '272' '274' 
% ,'199','269' '189' '110'
conditions =  {'1B','2B','3B'}; %should be the same number as the number of con images you're interrogating per subject.  These names will be put into the output file.
conimagenums = [2 3 4];
%

outputdir = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Activation/COMTXDATXLoad_92_INTERESTING/'; %where the output will be saved

outputfilename = 'ActivationOutput.txt';

delete([outputdir outputfilename]);
fid = fopen([outputdir outputfilename],'at'); %open the output file for writing
fprintf(fid,'%s\t\%s\t\%s\t\%s\t\%s\n\r\','Subject','Run','Condition','ROI','Value'); %write the output file header
fclose(fid);
dlmwrite([outputdir outputfilename],' ','-append');  %go to a new line in the output file


runs = {'Nback'};
    
 %roi_names = {'PCC_aal','vmPFC_aal'};

for run = 1:length(runs)




for subject = 1:length(subjects)

    subjsid = subjects{subject};
    
    data_dir = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subjsid '/SPM8/Cond/']; %where the con images are stored for each subject
    %
    %; 
        
    for condition = 1:length(conditions);
        
        %analysis_dir = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/Seed-based/DAT_DC_bilat_FirstRest_scrub/'];  %where the ROIs are
        analysis_dir = outputdir;
        
        clear roi_names
        allrois = dir([analysis_dir '*_roi.mat']);
        for i = 1:length(allrois)
            roi_names{i} = allrois(i).name(1:end-8);
        end
        
        
        
        for roinum = 1:length(roi_names)
            disp([subjsid ' : ' conditions{condition} ' : ' roi_names{roinum}])
            roiname = [analysis_dir, roi_names{roinum} , '_roi.mat'];
            rois = maroi('load_cell', roiname);

            getdataworked = 0;
            while getdataworked == 0
            try 
                Y = getdata(rois{1},[data_dir 'con_' sprintf('%04i',conimagenums(condition)) '.img'],'l'); %[data_dir conditions{condition} '_' runs{run} '_scrub_DCregress.nii'],'l');
                getdataworked = 1;
            catch
                disp('getdata failed')
            end
            end
            
            mean_within_roi = mean(Y(find(~isnan(Y))));
            texttowrite = [subjsid,'   ',runs{run},'   ',conditions{condition},'   ',roi_names{roinum},'   ',num2str(mean_within_roi)];  %save the data as a string to be written to the output
            dlmwrite([outputdir outputfilename],texttowrite,'-append','delimiter','');%write the data to the output file
        end
    end

end

end


cd(directory);
