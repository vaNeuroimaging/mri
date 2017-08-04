subjects = {'101','102','113','118','120','122','125','127','132','138','147','150','151','154','156','159','160','161','162','166','172','181','182','187','202','207','211','214','215','221','225','229','232','233','242','250','254','255','272','274','112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','247','295','297','300','305','309','334','339','340','343','359','362','383','395','396','374'};
%


outputname = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/Hip_Tracing/volumes.txt';

    delete(outputname);
    fid = fopen(outputname,'at');
    fprintf(fid,'%s\t\%s\t\%s\t\%s\t\%s\t\%s\t\%s\t\%s\t\%s\n\r\','Subject','Ritika_RH_1','Ritika_LH_1','Ritika_RH_2','Ritika_LH_2','Victoria_RH_1','Victoria_LH_1','Victoria_RH_2','Victoria_LH_2');
    fclose(fid);
    dlmwrite(outputname,' ','-append');
    

for subject = 1:length(subjects)
    subjid = subjects{subject};
    disp(subjid)
    
    Ritika_RH_1 = 'None';
    Ritika_LH_1 = 'None';
    try gunzip(['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subjid '/Struct/FIRST_Hip_seg/' subjid '_RC.nii.gz']); catch; end
    try    
        data = load_nii(['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subjid '/Struct/FIRST_Hip_seg/' subjid '_RC.nii']);
        Ritika_RH_1 = num2str(length(find(data.img==max(data.img(find(data.img))))));
        Ritika_LH_1 = num2str(length(find(data.img==min(data.img(find(data.img))))));
    catch
    end
    clear data

    Ritika_RH_2 = 'None';
    Ritika_LH_2 = 'None';
    try gunzip(['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subjid '/Struct/FIRST_Hip_seg/' subjid '_RC2.nii.gz']); catch; end
    try    
        data = load_nii(['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subjid '/Struct/FIRST_Hip_seg/' subjid '_RC2.nii']);
        Ritika_RH_2 = num2str(length(find(data.img==max(data.img(find(data.img))))));
        Ritika_LH_2 = num2str(length(find(data.img==min(data.img(find(data.img))))));
    catch
    end
    clear data
    
    Victoria_RH_1 = 'None';
    Victoria_LH_1 = 'None';
    try gunzip(['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subjid '/Struct/FIRST_Hip_seg/' subjid '_VT.nii.gz']); catch; end
    try    
        data = load_nii(['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subjid '/Struct/FIRST_Hip_seg/' subjid '_VT.nii']);
        Victoria_RH_1 = num2str(length(find(data.img==max(data.img(find(data.img))))));
        Victoria_LH_1 = num2str(length(find(data.img==min(data.img(find(data.img))))));
    catch
    end
    clear data

    Victoria_RH_2 = 'None';
    Victoria_LH_2 = 'None';
    try gunzip(['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subjid '/Struct/FIRST_Hip_seg/' subjid '_VT2.nii.gz']); catch; end
    try    
        data = load_nii(['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subjid '/Struct/FIRST_Hip_seg/' subjid '_VT2.nii']);
        Victoria_RH_2 = num2str(length(find(data.img==max(data.img(find(data.img))))));
        Victoria_LH_2 = num2str(length(find(data.img==min(data.img(find(data.img))))));
    catch
    end
    clear data
    
 
    
    texttowrite = [subjid '  ' Ritika_RH_1 '  ' Ritika_LH_1 '  ' Ritika_RH_2 '  ' Ritika_LH_2 '  ' Victoria_RH_1 '  ' Victoria_LH_1 '  ' Victoria_RH_2 '  ' Victoria_LH_2];
    dlmwrite(outputname,texttowrite,'-append','delimiter','');
    
%     try data = load_nii(['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subjid '/Struct/FIRST_Hip_seg/' subjid '_RC2.nii']);
%     catch
%         gunzip(['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subjid '/Struct/FIRST_Hip_seg/' subjid '_RC2.nii.gz']);
%         data = load_nii(['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subjid '/Struct/FIRST_Hip_seg/' subjid '_RC2.nii']);
%     end
%     
%     RCdata(subject,2) = length(find(data.img));
%     clear data
%     
%     
%     try data = load_nii(['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subjid '/Struct/FIRST_Hip_seg/' subjid '_VT.nii']);
%     catch
%         gunzip(['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subjid '/Struct/FIRST_Hip_seg/' subjid '_VT.nii.gz']);
%         data = load_nii(['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subjid '/Struct/FIRST_Hip_seg/' subjid '_VT.nii']);
%     end
%     
%     VTdata(subject,1) = length(find(data.img));
%     clear data
%     
%     
%     try data = load_nii(['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subjid '/Struct/FIRST_Hip_seg/' subjid '_VT2.nii']);
%     catch
%         gunzip(['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subjid '/Struct/FIRST_Hip_seg/' subjid '_VT2.nii.gz']);
%         data = load_nii(['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subjid '/Struct/FIRST_Hip_seg/' subjid '_VT2.nii']);
%     end
%     
%     VTdata(subject,2) = length(find(data.img));
%     clear data
    
end

% interrater = corr2(mean(RCdata,2),mean(VTdata,2));
% disp(['Interrater reliability: r = ' num2str(interrater)])
% 
% intraraterRC = corr2(RCdata(:,1),RCdata(:,2));
% disp(['Ritika''s Intrarater reliability: r = ' num2str(intraraterRC)])
% 
% intraraterVT = corr2(VTdata(:,1),VTdata(:,2));
% disp(['Victoria''s Intrarater reliability: r = ' num2str(intraraterVT)])

