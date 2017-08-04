%[~,~,output] = xlsread('ENIGMA_Subject_Info.xlsx');

[~,~,MAVoveralldata] = xlsread('/home/data/Analysis/ENIGMA/MAV_AssessmentData_cleaned_forENIGMA.xlsx','Data Summary');
[~,~,MAVPCLdata] = xlsread('/home/data/Analysis/ENIGMA/MAV_AssessmentData_cleaned_forENIGMA.xlsx','PCL-5');
[~,~,MAVhistdata] = xlsread('/home/data/Analysis/ENIGMA/MAV_AssessmentData_cleaned_forENIGMA.xlsx','Patient History');
[~,~,MAVCTQdata] = xlsread('/home/data/Analysis/ENIGMA/MAV_AssessmentData_cleaned_forENIGMA.xlsx','CTQ');
[~,~,MAVBDIdata] = xlsread('/home/data/Analysis/ENIGMA/MAV_AssessmentData_cleaned_forENIGMA.xlsx','BDI-II');
[~,~,MAVAUDITdata] = xlsread('/home/data/Analysis/ENIGMA/MAV_AssessmentData_cleaned_forENIGMA.xlsx','AUDIT');


ROBIdata = smartload('/home/data/EEG/processed/Robi/robiNeuropsych.mat');
ROBIdata = table2cell(ROBIdata);

subjects = dir('*_freesurfer');
outindex = 0;
for s = 1:length(subjects)
    subjectname = subjects(s).name(1:end-11);
    
    
    if ~isempty(find(strcmp(subjectname,MAVoveralldata(:,1)))) || ~isempty(find(strcmp(subjectname,ROBIdata(:,2))))
    outindex = outindex+1;
    
    output{outindex,1} = subjectname;
    output{outindex,4} = 'ND';
    output{outindex,5} = 'ND';
    output{outindex,6} = 'PCL-5';
    output{outindex,7} = 'ND';
    output{outindex,8} = 'ND';
    output{outindex,9} = 'ND';
    output{outindex,10} = 'ND';
    output{outindex,11} = 'ND';
    output{outindex,14} = 'ND';
    output{outindex,15} = 'ND';
    output{outindex,19} = 'ND';
    output{outindex,20} = 'ND';
    output{outindex,21} = 'ND';
    output{outindex,24} = 'ND';
    output{outindex,25} = 'ND';
    output{outindex,29} = 'ND';
    output{outindex,31} = 'ND';
    output{outindex,32} = 'Philips';
    output{outindex,33} = 'Achieva';
    output{outindex,34} = 16;
    output{outindex,35} = 'MP-RAGE';
    output{outindex,36} = '1x1x1';
    output{outindex,37} = '288x288';
    output{outindex,38} = 'sagittal';
    output{outindex,39} = 2400;
    output{outindex,40} = 3.08;
    output{outindex,41} = 90;
    output{outindex,43} = 1;
    output{outindex,44} = 'FS6.0';
    output{outindex,45} = 'CentOS6.7';
        
    
    if ~isempty(find(strcmp(subjectname,MAVoveralldata(:,1))))
        subindex = find(strcmp(subjectname,MAVPCLdata(:,1)));
        output{outindex,2} = MAVPCLdata{subindex,23};
        if strcmp(output{outindex,2},'No PTSD')
            output{outindex,2} = 'control';
        end
        output{outindex,3} = MAVPCLdata{subindex,22};
        
        subindex = find(strcmp(subjectname,MAVhistdata(:,1)) & (cellfun(@isstr,MAVhistdata(:,4))));
        if ~isempty(subindex)
            age = 2017 - str2num(MAVhistdata{subindex,4}(end-3:end));
            output{outindex,12} = age;
            
            genders = {'F','M'};
            output{outindex,13} = genders{str2num(MAVhistdata{subindex,6})};
            
            subindex = find(strcmp(subjectname,MAVhistdata(:,1)) & strcmp(MAVhistdata(:,14),'9'));
            if ~isempty(subindex)
                output{outindex,22} = MAVhistdata{subindex,15};
            else
                output{outindex,22} = 'no';
            end
            subindex = find(strcmp(subjectname,MAVhistdata(:,1)));
            temp = MAVhistdata(subindex,17)';
            medstring = [];
            for i = 1:length(temp)
                if ~isnan(temp{i})
                    medstring = [medstring ' ' temp{i}];
                end
            end
            output{outindex,23} = medstring;
            if isempty(medstring)
                output{outindex,23} = 'no';
            end
        end
        
        
        
        %race = {'white','white','white','white','white','white','white','white','white','white','white','Asian','Arabic','Native American','African American','Other'};
        %raceout = [];
        
        subindex = find(strcmp(subjectname,MAVCTQdata(:,1)));
        if isempty(MAVCTQdata{subindex,2})
            output{outindex,16} = 'ND';
            output{outindex,17} = 'ND';
            output{outindex,18} = 'ND';
        else
            numexperienced = nnz(strcmp('Moderate/Severe',MAVCTQdata(subindex,[31 33 35 37 39]))) + nnz(strcmp('Severe/Extreme',MAVCTQdata(subindex,[31 33 35 37 39])));
            numexperienced(numexperienced>2) = 2;
            output{outindex,16} = numexperienced;
            output{outindex,17} = max(cell2mat(MAVCTQdata(subindex,[30 32 34 36 38])));
            output{outindex,18} = 'CTQ';
        end
        
        
        subindex = find(strcmp(subjectname,MAVBDIdata(:,1)));
        if ~isempty(MAVBDIdata{subindex,2})
            output{outindex,26} = MAVBDIdata{subindex,26};
            output{outindex,27} = MAVBDIdata{subindex,25};
            output{outindex,28} = 'BDI-II';
        else
            output{outindex,26} = 'ND';
            output{outindex,27} = 'ND';
            output{outindex,28} = 'ND';
        end
        
        subindex = find(strcmp(subjectname,MAVAUDITdata(:,1)));
        output{outindex,30} = MAVBDIdata{subindex,12};
        
    elseif ~isempty(find(strcmp(subjectname,ROBIdata(:,2))));
        subindex = find(strcmp(subjectname,ROBIdata(:,2)));
        
        output{outindex,2} = 'PTSD';
        output{outindex,3} = sum(cell2mat(ROBIdata(subindex,285:304)));
        
        if ~isempty(ROBIdata{subindex,62})
            age = 2017 - (1900+str2num(ROBIdata{subindex,62}(end-1:end)));
            output{outindex,12} = age;
        end
        
        if ~isempty(ROBIdata{subindex,65}) && ~isnan(ROBIdata{subindex,65})
            genders = {'F','M'};
            output{outindex,13} = genders{ROBIdata{subindex,65}};
        end
        
        output{outindex,22} = ROBIdata{subindex,103};
        output{outindex,23} = ROBIdata{subindex,105};
        if isempty(output{outindex,23})
            output{outindex,23} = 'no';
        end
        
        output{outindex,16} = 'ND';
        output{outindex,17} = 'ND';
        output{outindex,18} = 'ND';
        output{outindex,26} = 'ND';
        output{outindex,27} = 'ND';
        output{outindex,28} = 'ND';
        output{outindex,30} = 'ND';
        
    end
    end
end

for i = 1:numel(output)
    if isempty(output{i}) || (isnumeric(output{i}) && isnan(output{i}))
        output{i} = 'ND';
    end
    if ischar(output{i}) && any(strfind(output{i},','))
        output{i}(strfind(output{i},',')) = '';
    end
end

copyfile('ENIGMA_Subject_Info_template.csv','ENIGMA_Subject_Info.csv')

for i = 1:size(output,1)
    outputtext = [];
    for j = 1:size(output,2)
        outputtext = [outputtext num2str(output{i,j}) ','];
    end

dlmwrite('ENIGMA_Subject_Info.csv',outputtext(1:end-1),'delimiter','','-append');
end
        
        
        
        
        
        
        
        