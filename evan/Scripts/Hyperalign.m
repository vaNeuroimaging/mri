subjects = {''};

maskfilename = '/home/usr/fidl/lib/glm_atlas_mask_222.4dfp.img';

conditions = {''};

Scratchdir = '';

maskdata = read_4dpfimg(maskfilename);


%Pass 1
for subject = 1:length(subjects)
    
    concatenateddata = [];
    
    for condition = 1:length(conditions)
        
        filename = [Scratchdir subjects{subject} dataprefix conditions{condition}];
    
        thiscondition = read_4dpfimg([filename '.4dfp.img']);
        concatenateddata = cat(2,concatenateddata,thiscondition);
    end
    
    data{subject} = concatenateddata(find(maskdata(:,1)),:);
    
    
    
    if subject==1
        meandata = data{subject};
    
    else
        
        [D, transformeddata] = PROCRUSTES(meandata, data{subject});
        
        meandata = (meandata .* (subject-1) + transformeddata) ./ subject;
        
    end
    
    clear transformeddata
    
end

%Pass2
for subject = 1:length(subjects)
    
    [D, transformeddata(:,:,subject)] = PROCRUSTES(meandata, data{subject});
    
end

finalmeandata = mean(transformeddata,3);

clear meandata


%Pass3
for subject = 1:length(subjects)
    
    [D, transformeddata, hyperparams{subject}] = PROCRUSTES(finalmeandata, data{subject});
    
end

clear 
    
    


