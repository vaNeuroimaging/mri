%setname = 'Motion-free';

%runs = {'Motor','Rest'};
%{'Switching','Motor','Rest'};

%for run = 1:length(runs)
    directory = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/GIFT_40/'];
    filestounzip = dir([directory '*.zip']);
    for file = 1:length(filestounzip)
        disp(filestounzip(file).name)
        mkdir([directory filestounzip(file).name(1:end-4)]);
        unzip([directory filestounzip(file).name],[directory filestounzip(file).name(1:end-4)]);
        
        cd([directory filestounzip(file).name(1:end-4)]);
        components = dir([filestounzip(file).name(1:end-4) '*.img']);
        files = [];
        for comp = 1:length(components)
            files = [files; [directory filestounzip(file).name(1:end-4) '/' components(comp).name ',1']];
        end
        addpath /fmri/spm5/toolbox/gift/icatb_helper_functions
        icatb_convert_to_z_shift(files);
        
    end
%end

