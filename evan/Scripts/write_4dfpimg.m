function write_4dfpimg(datamat,imgname,etype)


% what is endian type
switch etype
    case 'bigendian'
        etypespec=['ieee-be'];
    case 'littleendian'
        etypespec=['ieee-le'];
    otherwise
        error('Endian type selected was neither bigendian nor littleendian..');
end

fid=fopen(imgname,'w',etypespec);
fwrite(fid,datamat,'float32',0,etypespec);
fclose(fid);