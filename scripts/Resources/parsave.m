function parsave(filename,varargin)

v7p3 = false;
if strcmp(varargin{end},'-v7.3')
    v7p3 = true;
    varargin(end) = [];
end

string = ['save(''' filename ''','];
for i = 1:(length(varargin))
    eval([inputname(i+1) '=varargin{i};'])
    varargin{i} = [];
    string = [string '''' inputname(i+1) ''','];
end

if v7p3
    string = [string '''-v7.3'')'];
else
    string = [string(1:end-1) ')'];
end
    
eval(string);
