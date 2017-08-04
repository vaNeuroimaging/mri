%function call_jdp_preprocess(subject)
subject = 'MAV011';

[status, message] = system(['tcsh /home/data/scripts/JDP/scripts/JDP_P1_preprocess.sh ' subject]);

logfolder = ['/home/data/subjects/' subject '/processing_logs/'];
logfile = [logfolder 'jdp_preproc_' datestr(now) '.txt']; tokens = tokenize(logfile,' '); logfile = [tokens{1} '_' tokens{2}];
dlmwrite(logfile,message,'delimiter','')

warninglogfile = [logfile(1:end-4) '_error.txt'];
longwarninglogfile = [logfile(1:end-4) '_longerror.txt'];

system(['grep -v ''If you are performing spatial transformations on an oblique dset,'' ' logfile ' | grep -n ''WARNING\|Warning\|ERROR\|Error''  > ' warninglogfile]);
system(['grep -v ''If you are performing spatial transformations on an oblique dset,'' ' logfile ' | grep -n -A 30 -B 30 ''WARNING\|Warning\|ERROR\|Error''  > ' longwarninglogfile]);

