function cifti_write_wHDR(data,ciftitemplatefile,outnamestem)    
%cifti_write_wHDR(data,ciftitemplatefile,outnamestem)    

dir = pwd;
save(gifti(single(data)'),[outnamestem '.func.gii'],'ExternalFileBinary')

%Read template file and write out with correct header
bufsize = 524288;
ciftiheadertext = textread(ciftitemplatefile,'%s','delimiter','\r','bufsize',bufsize);
delete([outnamestem '.func.gii']);
fid = fopen([outnamestem '.func.gii'],'at'); %open the output file for writing
fclose(fid);

for row = 1:length(ciftiheadertext)
    thisline = ciftiheadertext{row};
    if length(thisline) >= 4 && strcmp(thisline(1:4),'Dim1')
        thisline = ['Dim1="' num2str(size(data,2)) '"'];
    elseif length(thisline) >= 16 && strcmp(thisline(1:16),'ExternalFileName')
        thisline = ['ExternalFileName="' outnamestem '.func.dat"'];
    elseif row>1 && strcmp(ciftiheadertext{row-1},'<Name>Provenance</Name>')
        thisline = '<Value>created with Matlab cifti_write function</Value>';
    end
    
    dlmwrite([outnamestem '.func.gii'],thisline,'-append','delimiter','');
    
end


system(['/data/cn4/laumannt/workbench/bin_linux64/wb_command -cifti-convert -from-gifti-ext ' dir '/' outnamestem '.func.gii ' dir '/' outnamestem '.dtseries.nii'])
delete([outnamestem '.func.*'])
