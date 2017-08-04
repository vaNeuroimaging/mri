%newclrs: regularized pattern file 
%b: number of thresholds
%varargin: subjective adjustment to the cutoff

function Final = Haoxin_MakingConsensus(newclrs, b, varargin)
%find the cutoff
for i = 1:b
    S(1,i) = sum([1:i]);
    S(2,i) = sum([b:-1:i-1]);
    if S(1,i)>S(2,i)
        if ~isempty(varargin)
            B = i+varargin{1};
        else
            B = i+1;
        end
        break
    end
end
% adjust the -1s
newclrs1 = newclrs;
for i = 1:size(newclrs,1)
    clear tp
    tp = find(newclrs(i,:)==-1);
    if ~isempty(tp)
        if tp(1)>b/3
            newclrs1(i,tp(1):end) = newclrs(i,tp(1)-1);
        else
            newclrs1(i,:) = -1;
        end
    end
end
% making the consensus 
Communities = unique(newclrs1)';
Communities(1) = [];
Final = zeros(size(newclrs1,1), 1);
for i = 1:length(Communities)
    Temp = zeros(size(newclrs));
    Temp(newclrs1 == Communities(i)) = Communities(i);
    firstPass = sum(logical(Temp)*(1:1:b)',2);
    firstPass(firstPass<=sum(B:b)) = 0;
    Final(find(firstPass)) = 0;
    Final = Final+Communities(i).*logical(firstPass);
end
Final(Final==0) = -1;