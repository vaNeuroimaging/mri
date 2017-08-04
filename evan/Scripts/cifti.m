function cifti(data,outnamestem)

ciftitemplatefile = '/data/cn4/evan/Scripts/cifti_template.func.gii';

bufsize = 524288;
ciftiheadertext = textread(ciftitemplate,'%s','delimiter','\r','bufsize',bufsize);

    save(gifti(data'),[outnamestem '.func.gii'],'ExternalFileBinary')
    
    delete([outnamestem '.func.gii']);
    fid = fopen([outnamestem '.func.gii'],'at'); %open the output file for writing
    fclose(fid);
    
    for row = 1:length(ciftiheadertext)
        thisline = ciftiheadertext{row};
        if strcmp(thisline(1:4),'Dim1')
            thisline = ['Dim1="' num2str(size(data,2)) '"'];
        elseif strcmp(thisline(1:16),'ExternalFileName')
            thisline = ['ExternalFileName="' outnamestem '.func.dat"'];
        elseif strcmp(ciftiheadertext{row-1},'<Name>Provenance</Name>')
            thisline = '<Value>created with Matlab cifti_write function</Value>';
        end
        
        dlmwrite([outnamestem '.func.gii'],thisline,'-append','delimiter','');
        
    end
    
    system(['wb_command -cifti-convert -from-gifti-ext 
        
    
    