function encode(File,File2)
%ENCODE  Converting ASCII text to Binary numbers 
%   DECODE(FILE,FILE2) where FILE is an input text file consisting of ASCII text and
%       FILE2 is the output file.
%       The code key file (code.m) should be included in the same directory in order
%       for this program to run.
%
%   Copyright 2004 Fahad Al Mahmood. 
%   $Revision: 1.00 $  $Date: 08-Jan-2004 09:05:10
%   $Revision: 1.50 $  $Date: 09-Feb-2004 15:21:40

wpl=0;
[spc_ent,C0]=textread('code.m','%s %s',2);
[L,C]=textread('code.m','%s %s','headerlines',2);
L=char(L);
FID = fopen(File,'r');
OUT = fopen(File2,'w');
while 1
    tline = fgetl(FID);
    if ~ischar(tline), break, end
    for i=1:size(tline,2)
        x=find(L==tline(i));        
        if isempty(x)==1
            fprintf(OUT,'%s',char(C0(1)));
            wpl=wpl+1;
        else
            fprintf(OUT,'%s',char(C(x)));
            wpl=wpl+1;
        end
        if wpl==10, fprintf(OUT,'\n');, wpl=0;, end
    end
    fprintf(OUT,char(C0(2)));
    wpl=wpl+1;
end
fclose('all');