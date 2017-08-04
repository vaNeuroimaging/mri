matlabversion = ['R' version('-release')];
warning off
delete(['~/.matlab/local_cluster_jobs/' matlabversion '/*'])
foldersleft = dir(['~/.matlab/local_cluster_jobs/' matlabversion '/']);
for f = 3:length(foldersleft)
    [success,message,mid] = rmdir(['~/.matlab/local_cluster_jobs/' matlabversion '/' foldersleft(f).name],'s');
    if ~success
        disp(message)
    end
end