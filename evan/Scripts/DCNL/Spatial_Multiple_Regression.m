%Created by E. Gordon 06/16/10

warning off

%USER MUST INPUT DATA
%--------------------------------------------------------------------------

%Define the participants whom you want to analyze
subjects = {'101','102','113','118','120','122','125','127','132','138','147','150','151','154','156','159','160','161','162','166','172','181','182','187','202','207','211','214','215','221','225','229','232','233','242','250','254','255','272','274'};


conditions = {'Allconds','1Back','2Back','3Back'};
conditionfiles = {'Cond/con_0007.hdr','Cond/con_0001.hdr','Cond/con_0002.hdr','Cond/con_0003.hdr'};


%Location of group ICA component images
GroupICAdir = '/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ICA_lowthresh20dim.gica/groupmelodic.ica/stats/';

%Non-artifactual group component numbers
validgroupcomponents = [1 3 6 7 8 10 12 13 14 16 19];

%ICA components begin with
filetype = 'thresh_zstat';

%Pearson's r value minimum to count as a valid match 
goodmatchthreshold = .25;

%Location of Marsbar whole-brain mask.  I recommend making this out of the mask image created by the group ICA procedure
brainroipath = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ICA_lowthresh20dim.gica/mask_roi.mat'];
brainroi = maroi(brainroipath);

%Location and name of text file to write out
outputfilename = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/Task_made_of_components/Task_made_of_components.txt'];

%END USER INPUT
%----------------------------------------------------------------

%Delete any previous output file, then create a new one and write its header
try delete(outputfilename); catch; end
fid = fopen(outputfilename,'at');
fprintf(fid,'%s\t\%s\t\%s\t\%s\t\%s\t\%s\n\r\','Subject','Condition','GroupComp','IndComp','GroupIndCorrelation','Beta','R-sq','pval');
fclose(fid);
dlmwrite(outputfilename,' ','-append');

%Subject loop
for subnum = 1:length(subjects)

    subj = subjects{subnum};
    disp(subj)

    %USER MUST INPUT DATA
    %--------------------------------------------------------------------------

    %Location of this subject's activation map
    activation_pattern = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/' conditionfiles{1}];


    %Location of this subject's individual ICA components
    ICAdir = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ICA_lowthresh20dim/' subj '.ica/filtered_func_data.ica/'];

    %Number of timepoints and TR of the rest scans
    ntime = 148;
    TR=2;

    %Location of t and f .txt files (should be in /report/)
    name_t = [ICAdir 'report/t'];
    name_f = [ICAdir 'report/f'];

    %END USER INPUT
    %----------------------------------------------------------------

    %How many components does this subject have?
    comp_count = dir([ICAdir 'stats/' filetype '*.img']);
    numica = length(comp_count);
    time_data = zeros(numica, ntime);

    %Read in this subject's components, in the same space as the activation map
    P = sprintf('%-130s',activation_pattern);
    for compnum = 1:numica
        P = [P;sprintf('%-130s',[ICAdir 'stats/' filetype num2str(compnum) '.hdr'])];
    end
    Y = getdata(brainroi,P)';
    TestingComponentdata = Y(:,2:end);

    %Get nyquist freq
    nfft = ntime/2;
    freq_data = zeros(numica, nfft);

    extension = '.txt';

    %Get frequency domain for all components
    for compnum = 1:numica

        % Acquiring time data
        file_name = [name_t, int2str(compnum), extension];
        data = textread(file_name);
        ntime = size(data,1);
        time_data(compnum,:) = (data(:,1))';

        % FFT of time data
        temp = abs(fft(time_data(compnum,:), ntime));
        freq_data(compnum,:) = temp(1:ntime/2);

    end

    %Initialize some vars
    max_value = zeros(numica,1);
    total_energy = zeros(numica,1);
    high_freq_noise = zeros(numica,1);
    energy_percent = zeros(numica, 1);
    valid_component = zeros(numica,1);
    number_of_valid_components = 0;
    len = size(freq_data,2); %nfft

    %Define thresholds for considering this a noise component: if high frquency noise energy is more than highfreq_thresh% of the total signal energy then the component is considered a noise component.
    %Threshold for % of energy
    highfreq_thresh = 50;
    %Threshold for being considered high frequency
    highfreq_hz_thresh = .1; %in Hz

    %Create some values for this procedure
    hf_init_index = (nfft*highfreq_hz_thresh)/(1/(2*TR));
    hf_init_index = round(hf_init_index);
    valid_component_index = [];


    %Determine validity of components using power spectra
    for num_component = 1:numica

        max_value(num_component) = max(freq_data(num_component,:));

        % Determining high frequency noise energy
        for column = hf_init_index:len
            high_freq_noise(num_component) = high_freq_noise(num_component) + freq_data(num_component, column)^2;
        end

        % take low freq power?
        % low_freq_power = sum(freq_data(num_component,1:hf_init_index-1).^2);
        for column = 1:hf_init_index-1
            total_energy(num_component) = total_energy(num_component) + freq_data(num_component, column)^2;
        end

        % Determining total energy
        total_energy(num_component) = total_energy(num_component) + high_freq_noise(num_component);

        %Determining energy percent
        energy_percent(num_component) = 100 * high_freq_noise(num_component)/total_energy(num_component);

        if energy_percent(num_component) <= highfreq_thresh;

            %define as a valid component
            valid_component_index(length(valid_component_index)+1) = num_component;

        end
    end

    %Restrict components to those identified as valid above
    TestingComponentdata = TestingComponentdata(:,valid_component_index);

    %Read in the group components, in the same space as the activation map
    P = sprintf('%-130s',activation_pattern);
    for compnum = 1:length(validgroupcomponents)
        P = [P;sprintf('%-130s',[GroupICAdir filetype num2str(validgroupcomponents(compnum)) '.hdr'])];
    end
    Y = getdata(brainroi,P)';
    Groupcomponentdata = Y(:,2:end);

    %Initialize some variables
    final_valid_components = [];
    bestmatch = zeros(1,length(validgroupcomponents));

    %Find best matches for group components
    for num_group_component = 1:length(validgroupcomponents)

        %Correlate all group vs all individual components
        Correlvals = corrcoef([Groupcomponentdata TestingComponentdata]);

        %Restrict the resulting correlation matrix to only group vs individual correlations
        Correlvals = Correlvals(1:length(validgroupcomponents),length(validgroupcomponents)+1:end);

        %Find the maximum correlation
        maxcorrelation = max(max(Correlvals));

        %find which group/individual pairing resulted in this maximum correlation
        [groupcomp, indcomp] = find(Correlvals == maxcorrelation);

        %Label this individual component as the best match for this group component
        bestmatch(groupcomp) = valid_component_index(indcomp);
        final_valid_components = [final_valid_components valid_component_index(indcomp)];
        %end

        %Make sure that the above group and individual components are no longer considered for subsequent matches
        Groupcomponentdata(:,groupcomp) = 0;
        TestingComponentdata(:,indcomp) = 0;

        %Save the correlation value
        maxcorrelations(groupcomp) = maxcorrelation;

    end

    %Sort the valid individual components
    final_valid_components = sort(final_valid_components);


    for condition = 1:length(conditions)

        activation_pattern = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/' conditionfiles{condition}];

        %Read in the activation map and the valid individual components, all in the same space as the activation map
        P = sprintf('%-130s',activation_pattern);
        for compnum = final_valid_components
            P = [P;sprintf('%-130s',[ICAdir 'stats/' filetype num2str(compnum) '.hdr'])];
        end
        Y = getdata(brainroi,P)';
        Componentdata = Y(:,2:end);
        
        for num_group_component = 1:length(validgroupcomponents)
            if maxcorrelations(num_group_component)<goodmatchthreshold & ~isempty(find(final_valid_components==bestmatch(num_group_component)));
                indcompnum = find(final_valid_components==bestmatch(num_group_component));
                Componentdata(:,indcompnum) = zeros(size(Componentdata,1),1);
            end
        end
        
        Activationdata = Y(:,1);

        %Regress the individual components against the activation map
        [Betas, Bint, Resid, Rint, Stats] = regress(Activationdata, [Componentdata]); %ones(length(Activationdata),1)

        for num_group_component = 1:length(validgroupcomponents)

            %Save each Beta value and identify which group-level component it was associated with
            if bestmatch(num_group_component) > 0
                thisBeta = num2str(Betas(find(final_valid_components==bestmatch(num_group_component))));
            else
                thisBeta = '0';
            end

            %Write out the Beta values
            texttowrite = [subj '  ' conditions{condition} '   ' num2str(validgroupcomponents(num_group_component)) '  ' num2str(bestmatch(num_group_component)) '  ' num2str(maxcorrelations(num_group_component)) '  ' thisBeta '   ' num2str(Stats(1)) '   ' num2str(Stats(3))];
            dlmwrite(outputfilename,texttowrite,'-append','delimiter','');

        end
    end
end
