function [etype] = endian_checker(imagename,varargin)

% this function checks endian type of a file. the user can either pass in
% the desired type and it will error out if the file isn't that type, or it
% will return the file type as etype
%
% jdp 6/8/10

if ~isempty(varargin)
    etype=varargin{1,1};
    switch etype
        case 'big'
        case 'little'
        otherwise
            error('Endian type selected was neither big nor little..');
    end
    filecommand = [ 'endian_4dfp ' imagename ];
    [ss,rr]=system(filecommand);
    pleasebetheregod=strfind(rr,etype);
    if isempty(pleasebetheregod)
        error('%s not %sendian, have to do this the oldfashioned way.',imagename,etype);
    end
else
    filecommand = [ 'endian_4dfp ' imagename ];
    [ss,rr]=system(filecommand);
    
    % start fishing for the endian type
    etype='big';
    pleasebetheregod=strfind(rr,etype);
    if isempty(pleasebetheregod)
        clear etype; etype='little';
        pleasebetheregod=strfind(rr,etype);
        if isempty(pleasebetheregod)
            error('%s appears to be neither big nor little endian',imagename);
        end
    end
end
