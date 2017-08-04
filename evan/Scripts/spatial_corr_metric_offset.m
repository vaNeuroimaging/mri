function spatial_corr_metric_offset(metric1name,metric2name,neighdist,alloffsetdist,outputdir,outputname)

for maxoffsetdist = 1:alloffsetdist

bufsize=16384;

a = dlmread(['/data/cn4/evan/fsaverage_LR32k/node_neighbors' num2str(neighdist) '.txt']);
maxneighbors = size(a,2);

outputstring = '[';
for i = 1:maxneighbors
    outputstring = [outputstring 'neighbors(:,' num2str(i) ') '];
end
outputstring = [outputstring(1:end-1) ']'];

eval([outputstring ' = textread([''/data/cn4/evan/fsaverage_LR32k/node_neighbors'' num2str(neighdist) ''.txt''],repmat(''%u '',[1 maxneighbors]),''delimiter'','' '',''bufsize'',bufsize,''emptyvalue'',NaN);']);

neighbors = neighbors+1;

[nextneighbors(:,1) nextneighbors(:,2) nextneighbors(:,3) nextneighbors(:,4)...
    nextneighbors(:,5) nextneighbors(:,6) nextneighbors(:,7)] = ...
    textread('/data/cn4/evan/fsaverage_LR32k/node_neighbors1.txt','%u %u %u %u %u %u %u',...
    'delimiter',' ','bufsize',bufsize,'emptyvalue',NaN);

nextneighbors = nextneighbors + 1;


metric1 = gifti(metric1name);
metric1 = metric1.cdata;

metric2 = gifti(metric2name);
metric2 = metric2.cdata;

offsetcorr_metric = zeros(size(metric1));


for i = 1:length(metric1)
    
    string{i} = ['Node ' num2str(i)];
    if i==1; fprintf('%s',string{i}); else fprintf([repmat('\b',1,length(string{i-1})) '%s'],string{i}); end
    
    nodeneigh = neighbors(i,:);
    nodeneigh(isnan(nodeneigh)) = [];
    
    metric1_val = metric1(nodeneigh);
    metric2_val = metric2(nodeneigh);
    
    nooffset_corrval = paircorr_mod(metric1_val,metric2_val);
    
    offsetcorr_metric(i) = 0;
    
    for direction = 2:7
        newnode = i;
        for offsetdist = 1:maxoffsetdist
            
            
            if isnan(nextneighbors(newnode,direction))
                newnode = nextneighbors(newnode,2);
            else
                newnode = nextneighbors(newnode,direction);
            end
            
            offsetnodeneigh = neighbors(newnode,:);
            
            offsetnodeneigh(isnan(offsetnodeneigh)) = [];
            
            offset_metric2_val = metric2(offsetnodeneigh);
            
            if length(nodeneigh) == length(offsetnodeneigh)
            
            offsetcorr = paircorr_mod(metric1_val,offset_metric2_val);
            
            if ((offsetcorr - nooffset_corrval) > offsetcorr_metric(i)) && (offsetcorr > .1) && (offsetdist == maxoffsetdist);
                offsetcorr_metric(i) = (offsetcorr - nooffset_corrval);
            end
            end
        end
    end
       
end

disp(' ')

  save(gifti(single(offsetcorr_metric)),[outputdir '/' outputname num2str(maxoffsetdist) '.func.gii']);
  
end
 
    
            

