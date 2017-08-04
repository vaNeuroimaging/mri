function decode(File,File2)
%DECODE  Converting Binary numbers to ASCII text 
%   DECODE(FILE,FILE2) where FILE is an input text file consisting of binary numbers
%       and FILE2 is the output ASCII text file.
%       The code key file (code.m) should be included in the same directory in order
%       for this program to run.
%
%   Copyright 2004 Fahad Al Mahmood. 
%   $Revision: 1.00 $  $Date: 08-Jan-2004 09:05:10
%   $Revision: 1.50 $  $Date: 09-Feb-2004 15:21:40


t=0;
[spc_ent,C0]=textread('code.m','%s %s',2);
[L,C]=textread('code.m','%s %s','headerlines',2);
L=char(L);
bit = length(char(C(1)));
FID = fopen(File,'r');
OUT = fopen(File2,'w');
while 1
    tline = fgetl(FID);
    if ~ischar(tline), break, end
    for i=1:size(tline,2)/bit
        for j=1:93
            x=isequal(char(C(j,:)),tline(8*i-7:8*i));
            if x==1
                t=j;
            end
        end
        if t~=0
            fprintf(OUT,'%s',L(t));
        elseif tline(8*i-7:8*i)==char(C0(2))
            fprintf(OUT,'\n');   
        elseif tline(8*i-7:8*i)==char(C0(1))
            fprintf(OUT,' ');
        end
        t=0;
    end
end
fclose('all');