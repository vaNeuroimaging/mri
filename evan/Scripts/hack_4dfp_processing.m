
folder = '/data/cn3/steven/CoE/';
subjects = {'MAV006','MAV008','SMN'};


for s = 2%1:length(subjects)
    subject = subjects{s};
    cd([folder '/' subject])
    contents = dir;
    for sess = 3:length(contents);
        sessions{sess-2} = contents(sess).name;
    end
    
    mkdir atlas
    mprage_size = [];
    
    avgstr = ['avgmpr_4dfp '];
    
    boldcounter = 0;
    mpragecounter = 0;
    for sess = 1:length(sessions)
        
        mprage_folders = dir([sessions{sess} '/NIFTI/*T1*']);
        
        for mpragenum = 1:length(mprage_folders)
            
            filename = dir([sessions{sess} '/NIFTI/' mprage_folders(mpragenum).name '/*.nii.gz']);
            
            if ~isempty(filename)
                    temp = load_untouch_nii([sessions{sess} '/NIFTI/' mprage_folders(mpragenum).name '/' filename(1).name]);
                    this_mprage_size = size(temp.img);

                    if this_mprage_size(3)==160
                        mpragecounter = mpragecounter+1;
                        %system(['niftigz_4dfp -4 -N ' sessions{sess} '/NIFTI/' mprage_folders(mpragenum).name '/' filename(1).name ' atlas/sess_mpr' num2str(mpragecounter)]);
                        avgstr = [avgstr ' ' folder '/' subject '/atlas/sess_mpr' num2str(mpragecounter) '.4dfp.img'];
                    end
                
                
                
            end
        end
        
        boldfolders = dir([sessions{sess} '/NIFTI/*resting_state*']);
        
        
        
        for boldnum = 1:length(boldfolders)
            
            filename = dir([sessions{sess} '/NIFTI/' boldfolders(boldnum).name '/*.nii.gz']);
            
            if ~isempty(filename)
                
                temp = load_untouch_nii([sessions{sess} '/NIFTI/' boldfolders(boldnum).name '/' filename.name]);
                this_mprage_size = size(temp.img);
                
                if all(this_mprage_size(1:3)==[64 64 34])
                    
                    boldcounter = boldcounter + 1;
                                        
                    mkdir(['bold' num2str(boldcounter)])
                    
                    system(['niftigz_4dfp -4 -N ' sessions{sess} '/NIFTI/' boldfolders(boldnum).name '/' filename.name ' bold' num2str(boldcounter) '/' subject '_b' num2str(boldcounter)])
                    
                    [voxelsize frames I J K etype] = read_4dfpifh_HCP(['bold' num2str(boldcounter) '/' subject '_b' num2str(boldcounter) '.4dfp.ifh']);
                    
                    if ~all(str2num([I J K])'==[64 64 34])
                        system(['rm -R bold ' num2str(boldcounter)])
                        boldcounter = boldcounter - 1;
                    end
                    
                end
            end
        end
    end
    
    
    cd atlas
    avgstr = [avgstr ' ' folder '/' subject '/atlas/' subject '_mpr_n1.4dfp.img'];
    %system(avgstr)
    %system('t4_inv sess_mpr1_to_711-2B_t4')
    %system(['t4img_4dfp 711-2B_to_sess_mpr1_t4 ' folder '/' subject '/atlas/' subject '_mpr_n1_111_t88.4dfp.img  ' folder '/' subject '/atlas/' subject '_mpr1.4dfp.img -O' folder '/' subject '/atlas/sess_mpr1.4dfp.img'])
    
    cd ..
    outputfilename = [subject '_preprocess.prm'];
        
    delete([outputfilename]);
    fid = fopen([outputfilename],'at'); %open the output file for writing
    fclose(fid);
    
    linestowrite = {'1.0','TRIO_KY_NDC',subject,[folder '/' subject],[subject '.MR.HEAD_FRAN'],'VB17','0','1','0',num2str([1:boldcounter]),...
        ['1-' num2str(boldcounter)],[num2str(boldcounter) 'x2.5'],[num2str(boldcounter) 'x.078125'],'1','5','0','0.000590','even'};
    for line = 1:length(linestowrite)
        dlmwrite([outputfilename],linestowrite{line},'-append','delimiter','');
    end
    
    system(['/data/cn4/evan/Scripts/Processing_pipeline/preprocess_evan ' outputfilename])
    
end