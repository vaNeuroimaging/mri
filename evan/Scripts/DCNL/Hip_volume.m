subjects = {'113','118','120','122','125','126'};
%{'101','102','113','118','120','122','125','127','132','138','147','150','151','154','156','159','160','161','162','166','172','181','182','187','202','207','211','214','215','221','225','229','232','233','242','250','254','255','272','274','112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','281','283','292','301','307','322','327','247','295','297','300','305','309','334','339','340','343','359','362','383','395','396','374','400','401','402','406','407','410','412','415','416','417','420','199','269','189','110'};
%


outputname = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/Hip_Tracing/volumes.txt';

    delete(outputname);
    fid = fopen(outputname,'at');
    fprintf(fid,'%s\t\%s\t\%s\n\r\','Subject','RH','LH','TCV');
    fclose(fid);
    dlmwrite(outputname,' ','-append');
    

for subject = 1:length(subjects)
    subjid = subjects{subject};
    disp(subjid)
    
    RH = 'None'; LH = 'None';
    
    try gunzip(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/Hip_Tracing/' subjid '_RC2.nii.gz']); 
    catch
        gunzip(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/Hip_Tracing/' subjid '_RC.nii.gz']);
    end
    try    
        data = load_nii(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/Hip_Tracing/' subjid '_RC2.nii']);
    catch
        data = load_nii(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/Hip_Tracing/' subjid '_RC.nii']);
    end
    RH = num2str(length(find(data.img==max(data.img(find(data.img))))));
    LH = num2str(length(find(data.img==min(data.img(find(data.img))))));
    clear data
    
    c1file = dir(['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subjid '/Struct/c1D*.img']);
    c2file = dir(['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subjid '/Struct/c2D*.img']);
    c3file = dir(['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subjid '/Struct/c3D*.img']);
    if isempty(c1file)
        c1file = dir(['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subjid '/Struct/c1S*.nii']);
        c2file = dir(['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subjid '/Struct/c2S*.nii']);
        c3file = dir(['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subjid '/Struct/c3S*.nii']);
    end

    
    t = VBM_get_TCV(['/fmri/data3/Evan/Gene-Rest-Nback/Data/' subjid '/Struct/' c1file(1).name ',1'; '/fmri/data3/Evan/Gene-Rest-Nback/Data/' subjid '/Struct/' c2file(1).name ',1'; '/fmri/data3/Evan/Gene-Rest-Nback/Data/' subjid '/Struct/' c3file(1).name ',1'],.3);
        
    TCV = num2str(sum(t));

    texttowrite = [subjid '  ' RH '  ' LH '  ' TCV];
    dlmwrite(outputname,texttowrite,'-append','delimiter','');
    
    
end


