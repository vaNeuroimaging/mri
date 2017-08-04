subjects = {'101','102','110','113','118','120','122','125','127','132','138','147','150','151','154','156','159','160','161','162','166','172','181','182','187','189','199','202','207','211','214','215','221','225','229','232','233','242','250','254','255','269','272','274','112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','247','295','297','300','305','309','334','339','340','343','359','362','383','395','396','374','400','401','402','406','407','410','412','415','416','417','420'};
%{'101','102','110','113','118','120','122','125','127','132','138','147','150','151','154','156','159','160','161','162','166','172','181','182','187','189','199','202','207','211','214','215','221','225','229','232','233','242','250','254','255','269','272','274','112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','247','295','297','300','305','309','334','339','340','343','359','362','383','395','396','374','400','401','402','406','407','410','412','415','416','417','420'};
%{}
%
%

%Names of runs to interrogate
runs = {'FirstRest','Rest','Nback'};


outputdir = '/fmri/data3/Evan/Gene-Rest-Nback/Data/Headmotion/'; %where the output will be saved

outputfilename = 'FramewiseDisplacement.txt';

delete([outputdir outputfilename]);
fid = fopen([outputdir outputfilename],'at'); %open the output file for writing
fprintf(fid,'%s\t\%s\t\%s\n\r\','Subject','Run','TR','FD'); %write the output file header
fclose(fid);
dlmwrite([outputdir outputfilename],' ','-append');  %go to a new line in the output file



for subject = 1:length(subjects)
    
    subjsid = subjects{subject};
    
    for run = 1:length(runs)
        
        data_dir = ['/fmri/data3/Evan/Gene-Rest-Nback/Data/', subjsid, '/SPM8/' runs{run} '/'];
        
        if exist(data_dir,'dir')
            
            disp(['Subject ' subjsid ': ' runs{run}])
            
            FD = Calc_FD(data_dir);
            
            for time = 1:length(FD)
                
                texttowrite = [subjsid,'   ',runs{run},'   ',num2str(time),'   ',num2str(FD(time))];  %save the data as a string to be written to the output
                success = 0;
                while success == 0;
                try dlmwrite([outputdir outputfilename],texttowrite,'-append','delimiter','');%write the data to the output file
                    success = 1;
                catch
                end
                end
            end
        end
    end
end