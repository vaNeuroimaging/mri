function data = cifti_read(filename)
%function data = cifti_read(filename)

slashloc = strfind(filename,'/');
if isempty(slashloc)
    slashloc = 0;
end

result = evalc(['!wb_command -cifti-convert -to-gifti-ext ' filename ' ' filename(slashloc(end)+1:end-13) '.func.gii']);
    
try
    data = gifti([filename(slashloc(end)+1:end-13) '.func.gii']);
catch
    try
    data = gifti([filename(slashloc(end)+1:end-13) '.func.gii.data']);
    catch
        error(result)
    end
end
data = data.cdata;
delete([filename(slashloc(end)+1:end-13) '.func.gii'])
delete([filename(slashloc(end)+1:end-13) '.func.gii.data'])