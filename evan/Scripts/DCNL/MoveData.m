subjects = {'374'};
    %'112','126','133','137','208','222','227','251','253','256','258','260','261','264','270','279','283','281','292','301','307','322','327','247','295','297','300','305','309','334','339','340','343','359','362','383','395','396'};
%{'101','102','113','118','120','122','125','127','132','138','147','150','151','154','156','159','160','161','162','166','172','187','202','207','211','214','215','221','225','229','232','233','242','250','254','255','272','274'};

 rootdir = '/fmri/data3/Evan/Gene-Rest-Nback/Data/';

for sub = 1:length(subjects);
    
    disp(['Subject ' subjects{sub}])

%     mkdir([rootdir subjects{sub} '/SPM8/']);
%     mkdir([rootdir subjects{sub} '/SPM8/Nback/']);
%     mkdir([rootdir subjects{sub} '/SPM8/Rest/']);
%     mkdir([rootdir subjects{sub} '/SPM8/FirstRest/']);

     copyfile([rootdir subjects{sub} '/avw/*Nbac*'], [rootdir subjects{sub} '/SPM8/Nback/']);
     delete([rootdir subjects{sub} '/SPM8/Nback/*00001*']); delete([rootdir subjects{sub} '/SPM8/Nback/*00002*']);
     
     Restfiles = dir([rootdir subjects{sub} '/avw/*Rest*']);
     for i = 1:(length(Restfiles)/2)
         copyfile([rootdir subjects{sub} '/avw/' Restfiles(i).name], [rootdir subjects{sub} '/SPM8/FirstRest/']);
     end
     delete([rootdir subjects{sub} '/SPM8/FirstRest/*00001*']); delete([rootdir subjects{sub} '/SPM8/FirstRest/*00002*']);
     
     for i = (length(Restfiles)/2)+1:length(Restfiles)
         copyfile([rootdir subjects{sub} '/avw/' Restfiles(i).name], [rootdir subjects{sub} '/SPM8/Rest/']);
     end
     delete([rootdir subjects{sub} '/SPM8/Rest/*00001*']); delete([rootdir subjects{sub} '/SPM8/Rest/*00002*']);
     

     copyfile([rootdir subjects{sub} '/avw/*MPRAGE*'], [rootdir subjects{sub} '/Struct/']);    
     copyfile([rootdir subjects{sub} '/avw/*' subjects{sub} '_run*.SER'], [rootdir subjects{sub} '/DTI/']);
end
