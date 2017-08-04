function FD = Calc_FD(data_dir)


           %Find and load motion param file
           motionparamfile = dir([data_dir 'rp*.txt']);
           motionparams = textread([data_dir motionparamfile.name]);
           
           FD(1) = 0;
           
           
           for i = 2:size(motionparams,1)
            FD(i) = abs(motionparams(i,1) - motionparams(i-1,1)) + abs(motionparams(i,2) - motionparams(i-1,2)) + abs(motionparams(i,3) - motionparams(i-1,3)) + ...
                abs(sin(motionparams(i,4)/2)*100 - sin(motionparams(i-1,4)/2)*100) + abs(sin(motionparams(i,5)/2)*100 - sin(motionparams(i-1,5)/2)*100) + abs(sin(motionparams(i,6)/2)*100 - sin(motionparams(i-1,6)/2)*100);
            
           end
    

end

