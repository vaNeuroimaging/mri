subs = {'161' '166' '102' '113' '172' '187' '207' '211' '214' '229' '232' '233' '242' '250' '254' '255' '272' '274' '118' '120' '122' '125' '127' '132' '138' '147' '150' '151' '154' '156' '159' '160' '162' '181' '182' '202' '215' '221' '225'};
%'101'
%'189','110' '269' '199' 

rootdir = '/fmri/data3/Evan/Gene-Rest-Nback/Data/';
destdir = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ICA_lowthresh30dim/';

warning off

unixcodefile = '/fmri/data3/Evan/Gene-Rest-Nback/Scripts/unixcode1';

for sub = 1:length(subs)
    disp(subs{sub})
    
    
    
    
%     try rmdir(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subs{sub} '/DualRegression_lowdim'],'s'); catch; end
%     try rmdir(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subs{sub} '/DualRegression_Nback'],'s'); catch; end
%     try rmdir(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subs{sub} '/DualRegression_Nback_byload'],'s'); catch; end
%     try rmdir(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subs{sub} '/DualRegression_Nback_lowdim'],'s'); catch; end
%     try rmdir(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subs{sub} '/Nback_rdlPFC_DR'],'s'); catch; end
%     try rmdir(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subs{sub} '/Rest_rdlPFC_DR'],'s'); catch; end
%     try rmdir(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subs{sub} '/PCC_seed_GSR'],'s'); catch; end
%     try rmdir(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subs{sub} '/Nback_caudate'],'s'); catch; end
%     try rmdir(['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subs{sub} '/Rest_caudate'],'s'); catch; end

%     try
%         
%         
%         
%         movefile([rootdir subs{sub} '/Rest/Rest_ICA_lowthresh30dim.ica'], [destdir subs{sub} '.ica']);
%         %rmdir([rootdir subs{sub} '/Rest/Rest_lowthresh30dim.ica'],'s');
%     catch
%                  disp('failed to move')
%     end
    
    fid = fopen(unixcodefile,'w+');
    
    files = dir([destdir subs{sub} '.ica/filtered_func_data.ica/stats/thresh_zstat*.nii.gz']);
    for file = 1:length(files)
        command = ['fslchfiletype ANALYZE ' destdir subs{sub} '.ica/filtered_func_data.ica/stats/' files(file).name];
        fprintf(fid,'%s\n',command);
    end
    
    flag=fclose(fid); 
    
    
    
    
    eval(['!chmod u+rwx ' unixcodefile])
    
    pause(.5)
    
    eval(['!' unixcodefile])
    
    delete(unixcodefile);
 
end