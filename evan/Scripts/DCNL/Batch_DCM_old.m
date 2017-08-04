

warning off

directory = pwd;

outputfilename = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/DCM/DCMOutput.txt';

subjects = {'101'};
    %,'102','113','118','120','122','125','127','132','138','147','150','151','154','156','159','160','161','162','166','172','187','202','207','211','214','215','221','225','229','232','233','242','250','254','255','272','274'};

regions = {'vmPFC','PCC','r_IPL','r_ilPFC','r_dlPFC','ACC','l_aIns','r_precun'};
%{'rilPFC','PCC','VS'};%caudate ventralstriatum putamen
centers = {[-4 56 -8] [4 -60 24] [-36 -60 38] [32 60 0] [48 28 24] [4 32 20] [-44 12 -4] [12 -64 56]};
%{[32 60 0],[4,-60,24],[6 16 -8]};% 9,9,-8 25,8,6
models = {'ForFig'};


delete([outputfilename]);
fid = fopen([outputfilename],'at');
fprintf(fid,'%s\t\%s\t\%s\t\%s\t\%s\n\r\','Subject','ROI_connection','ModulatedBy','Parameter','Model');
fclose(fid);
dlmwrite([outputfilename],' ','-append');

for modelnum = 1:length(models)
    model = models{modelnum};
    
    for subject = 1:length(subjects)
        
        disp(subjects{subject})
        
        [SPM, xSPM] = spm_getSPM_Evan(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subjects{subject} '/SPM.mat']);
        
        for region = 1:length(regions)
            
            xY.xyz = centers{region};
            xY.name = regions{region};
            xY.Ic = 0;
            xY.Sess = 1;
            xY.def = 'sphere';
            xY.spec = 6;
            
            cd(directory)
            
            [Y,xY] = spm_regions_Evan(xSPM,SPM,xY);
            
            clear xY
            
        end
        
        
        load(['DCM_' model '.mat'])
        for region = 1:length(regions)
            load(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subjects{subject} '/VOI_' regions{region} '_1.mat']);
            
            DCM.Y.X0 = xY.X0;
            DCM.Y.y(:,region) = xY.u;
            DCM.Y.name{region} = regions{region};
            DCM.xY(region) = xY;
            
        end
        
        save(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subjects{subject} '/DCM_' model],'DCM');
        
        clear DCM
        
        spm_dcm_estimate(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subjects{subject} '/DCM_' model]);
        
        load(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subjects{subject} '/DCM_' model]);
        
        
        %
            texttowrite = [subjects{subject},'   ',[regions{2} '_to_' regions{1}],'   ',regions{3},'   ',num2str(DCM.D(1,2,3)),'   ',model];
            dlmwrite([outputfilename],texttowrite,'-append','delimiter','');
        
        %     texttowrite = [subjects{subject},'   ',[regions{2} '_to_' regions{1}],'   ','none','   ',num2str(DCM.A(1,2)),'   ',model];
        %     dlmwrite([outputfilename],texttowrite,'-append','delimiter','');
        
%         texttowrite = [subjects{subject},'   ',[regions{1} '_to_' regions{3}],'   ','none','   ',num2str(DCM.A(3,1)),'   ',model];
%         dlmwrite([outputfilename],texttowrite,'-append','delimiter','');
%         
%         texttowrite = [subjects{subject},'   ',[regions{2} '_to_' regions{1}],'   ','none','   ',num2str(DCM.A(1,2)),'   ',model];
%         dlmwrite([outputfilename],texttowrite,'-append','delimiter','');
%         
%         texttowrite = [subjects{subject},'   ',[regions{3} '_to_' regions{4}],'   ','none','   ',num2str(DCM.A(4,3)),'   ',model];
%         dlmwrite([outputfilename],texttowrite,'-append','delimiter','');
%         
%         texttowrite = [subjects{subject},'   ',[regions{4} '_to_' regions{2}],'   ','none','   ',num2str(DCM.A(2,4)),'   ',model];
%         dlmwrite([outputfilename],texttowrite,'-append','delimiter','');
%         
        
%P(subject,:) = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subjects{subject} '/DCM_' model '.mat'];


    end
    
%     cd /fmri/data3/Evan/Gene-Rest-Nback/Analysis/SPM_group/DCM/
%     spm_dcm_average(1,P,['Group' model '.mat']);
%     clear P
%     cd(directory)
    
end





