function make_cifti_label(ciftifile)
% make_cifti_label(ciftifile)

cifti = ft_read_cifti_mod(ciftifile); cifti = cifti.data;
IDs = unique(cifti); IDs(IDs==0) = [];

colors = round(distinguishable_colors(length(IDs)) .* 255);

warning off
delete('labellist.txt');
fid = fopen('labellist.txt','at'); %open the output file for writing
fclose(fid);

for i = 1:length(IDs)
    dlmwrite('labellist.txt',['Label_' num2str(IDs(i))],'-append','delimiter','');
    dlmwrite('labellist.txt',num2str([IDs(i) colors(i,:) 1]),'-append','delimiter','');
end

file_parts = tokenize(ciftifile,'.');
outname = [];
for i = 1:(length(file_parts)-2)
    outname = [outname file_parts{i} '.'];
end
outname = [outname 'dlabel.nii'];

[~,~] = system(['wb_command -cifti-label-import ' ciftifile ' labellist.txt ' outname]);
