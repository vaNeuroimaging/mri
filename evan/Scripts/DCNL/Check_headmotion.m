subs = {'101','102','110','113','118','120','122','125','127','132','138','147','150','151','154','156','159','160','161','162','166','172','181','182','187','189','199','202','207','211','214','215','221','225','229','232','233','242','250','254','255','269','272','274','112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','247','295','297','300','305','309','334','339','340','343','359','362','383','395','396','374','400','401','402','406','407','410','412','415','416','417','420'};

runs = {'Nback','Rest','FirstRest'};

outputname = '/fmri/data3/Evan/Gene-Rest-Nback/Data/Headmotion/Headmotion.txt';

    delete(outputname);
    fid = fopen(outputname,'at');
    fprintf(fid,'%s\t\%s\t\%s\t\%s\t\%s\t\%s\t\%s\t\%s\t\%s\t\%s\t\%s\t\%s\t\%s\t\%s\t\%s\t\%s\t\%s\n\r\','Run','Subject','TR','X','Y','X','Roll','Pitch','Yaw','Abs_X','Abs_Y','Abs_Z','Abs_Roll_degrees','Abs_Pitch_degrees','Abs_Yaw_degrees','Euclidean_translation','Euclidean_rotation');
    fclose(fid);
    dlmwrite(outputname,' ','-append');


for run = 1:length(runs)
    
    warning off
    for subnum = 1:length(subs)
        subj = subs{subnum};
        disp([subj ': ' runs{run}])
        hmfile = dir(['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subj '/SPM8/' runs{run} '/rp*.txt']);
        
        if ~isempty(hmfile)
            data = textread(['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subj '/SPM8/' runs{run} '/' hmfile.name]);
            for TR = 1:size(data,1)
                texttowrite = [runs{run} '  ' subj '  ' num2str(TR) '  ' num2str(data(TR,1)) '  ' num2str(data(TR,2)) '  ' num2str(data(TR,3)) '  ' num2str(data(TR,4)) '  ' num2str(data(TR,5)) ...
                     '  ' num2str(data(TR,6)) '  ' num2str(abs(data(TR,1))) '  ' num2str(abs(data(TR,2))) '  ' num2str(abs(data(TR,3))) '  ' num2str(abs(data(TR,4)*180/pi)) '  ' num2str(abs(data(TR,5)*180/pi)) '  ' num2str(abs(data(TR,6)*180/pi))...
                     '  ' num2str(sqrt((data(TR,1)^2 + data(TR,2)^2 + data(TR,3)^2))) '  ' num2str(sqrt((data(TR,4)^2 + data(TR,5)^2 + data(TR,6)^2)))];
                dlmwrite(outputname,texttowrite,'-append','delimiter','');
            end
        end
        
    end
    
end
