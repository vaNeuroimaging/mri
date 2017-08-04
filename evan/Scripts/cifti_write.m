function cifti_write(data,outnamestem)
%cifti_write(data,outnamestem)
%Write the contents of the 66697xN matrix "data" to the file
%"outnamestem.dtseries.nii"

ciftitemplatefile = '/data/cn4/evan/Scripts/cifti_template.func.gii';

bufsize = 524288;
ciftiheadertext = textread(ciftitemplatefile,'%s','delimiter','\r','bufsize',bufsize);

save(gifti(data'),[outnamestem '.func.gii'],'ExternalFileBinary')

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

evalc(['!wb_command -cifti-convert -from-gifti-ext ' outnamestem '.func.gii ' outnamestem '.dtseries.nii']);

delete([outnamestem '.func.gii'])
delete([outnamestem '.func.dat'])


