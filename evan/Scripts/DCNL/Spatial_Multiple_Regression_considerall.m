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

%Total group components
numbergroupcomponents = 20;

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
fprintf(fid,'%s\t\%s\t\%s\t\%s\t\%s\n\r\','Subject','Condition','GroupComp','IndComp','GroupIndCorrelation','Beta','R-sq');
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
    ICAdir = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/ICA_lowthresh30dim/' subj '.ica/filtered_func_data.ica/'];

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
    
    valid_component_index = 1:numica;

%     %Get nyquist freq
%     nfft = ntime/2;
%     freq_data = zeros(numica, nfft);
% 
%     extension = '.txt';
% 
%     %Get frequency domain for all components
%     for compnum = 1:numica
% 
%         % Acquiring time data
%         file_name = [name_t, int2str(compnum), extension];
%         data = textread(file_name);
%         ntime = size(data,1);
%         time_data(compnum,:) = (data(:,1))';
% 
%         % FFT of time data
%         temp = abs(fft(time_data(compnum,:), ntime));
%         freq_data(compnum,:) = temp(1:ntime/2);
% 
%     end
% 
%     %Initialize some vars
%     max_value = zeros(numica,1);
%     total_energy = zeros(numica,1);
%     high_freq_noise = zeros(numica,1);
%     energy_percent = zeros(numica, 1);
%     valid_component = zeros(numica,1);
%     number_of_valid_components = 0;
%     len = size(freq_data,2); %nfft
% 
%     %Define thresholds for considering this a noise component: if high frquency noise energy is more than highfreq_thresh% of the total signal energy then the component is considered a noise component.
%     %Threshold for % of energy
%     highfreq_thresh = 50;
%     %Threshold for being considered high frequency
%     highfreq_hz_thresh = .1; %in Hz
% 
%     %Create some values for this procedure
%     hf_init_index = (nfft*highfreq_hz_thresh)/(1/(2*TR));
%     hf_init_index = round(hf_init_index);
%     valid_component_index = [];
% 
% 
%     %Determine validity of components using power spectra
%     for num_component = 1:numica
% 
%         max_value(num_component) = max(freq_data(num_component,:));
% 
%         % Determining high frequency noise energy
%         for column = hf_init_index:len
%             high_freq_noise(num_component) = high_freq_noise(num_component) + freq_data(num_component, column)^2;
%         end
% 
%         % take low freq power?
%         % low_freq_power = sum(freq_data(num_component,1:hf_init_index-1).^2);
%         for column = 1:hf_init_index-1
%             total_energy(num_component) = total_energy(num_component) + freq_data(num_component, column)^2;
%         end
% 
%         % Determining total energy
%         total_energy(num_component) = total_energy(num_component) + high_freq_noise(num_component);
% 
%         %Determining energy percent
%         energy_percent(num_component) = 100 * high_freq_noise(num_component)/total_energy(num_component);
% 
%         if energy_percent(num_component) <= highfreq_thresh;
% 
%             %define as a valid component
%             valid_component_index(length(valid_component_index)+1) = num_component;
% 
%         end
%     end

    %Restrict components to those identified as valid above
    TestingComponentdata = TestingComponentdata(:,valid_component_index);

    %Read in the group components, in the same space as the activation map
    P = sprintf('%-130s',activation_pattern);
    for compnum = 1:numbergroupcomponents
        P = [P;sprintf('%-130s',[GroupICAdir filetype num2str(compnum) '.hdr'])];
    end
    Y = getdata(brainroi,P)';
    Groupcomponentdata = Y(:,2:end);
    
    Correlvals = corrcoef([Groupcomponentdata TestingComponentdata]);
    
    Correlvals = Correlvals(1:numbergroupcomponents,numbergroupcomponents+1:end);
    
    final_valid_components = [];
    bestmatch = [];
    maxcorrelations = [];
    best_of_multiple_matches = [];
    
    for num_component = 1:length(valid_component_index)
        
        matchinggroupcomp = find(Correlvals(:,num_component)==max(Correlvals(:,num_component)));
        
        if any(find(validgroupcomponents==matchinggroupcomp)) & max(Correlvals(:,num_component)) > goodmatchthreshold;
            
            final_valid_components = [final_valid_components num_component];
    
            bestmatch = [bestmatch matchinggroupcomp];
            
            maxcorrelations = [maxcorrelations max(Correlvals(:,num_component))];
            
        end
        
    end

    for num_component = 1:length(final_valid_components)
        
        if max(maxcorrelations(find(bestmatch==bestmatch(num_component))))==maxcorrelations(num_component)
            best_of_multiple_matches(num_component) = 1;
        else
            best_of_multiple_matches(num_component) = 0;
        end
    end

    

    for condition = 1:length(conditions)

        activation_pattern = ['/fmri/data3/Evan/Gene-Rest-Nback/Analysis/' subj '/' conditionfiles{condition}];

        %Read in the activation map and the valid individual components, all in the same space as the activation map
        P = sprintf('%-130s',activation_pattern);
        for compnum = final_valid_components
            P = [P;sprintf('%-130s',[ICAdir 'stats/' filetype num2str(compnum) '.hdr'])];
        end
        Y = getdata(brainroi,P)';
        Componentdata = Y(:,2:end);
        
        Activationdata = Y(:,1);

        %Regress the individual components against the activation map
        [Betas, Bint, Resid, Rint, Stats] = regress(Activationdata, [Componentdata(:,find(best_of_multiple_matches)) ones(length(Activationdata),1)]);

        used_comp_counter = 0;
        for num_component = 1:length(final_valid_components)
            
            if best_of_multiple_matches(num_component)
                used_comp_counter = used_comp_counter +1;
                ThisBeta = num2str(Betas(used_comp_counter));
            else
                ThisBeta = 'NotUsed';
            end
            %Write out the Beta values
            texttowrite = [subj '  ' conditions{condition} '   ' num2str(bestmatch(num_component)) '  ' num2str(final_valid_components(num_component)) '  ' num2str(maxcorrelations(num_component)) '  ' ThisBeta '   ' num2str(Stats(1))];
            dlmwrite(outputfilename,texttowrite,'-append','delimiter','');

        end
    end
end
