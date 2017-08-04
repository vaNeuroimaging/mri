%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Inputs    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

data = textread('/data/cn4/evan/Occipitotemporal/Timecourses_for_clustering.txt');        % text file where rows are region timecourses and columns are datapoints (r x d)
names = {'lDA_iOT' 'lDA_IPS1' 'lDA_IPS2' 'lDA_IPS3' 'lDA_ParOcc' 'lDA_SMA' 'lDA_sOT' 'lDA_SPL' 'lFPC_aIFG' 'lFPC_midIFG' 'lFPC_OT' 'lFPC_peanut' 'lFPC2_aIFG' 'lFPC2_dmPFC' 'lFPC2_IPL' 'lFPC2_pMFG' 'lFPC2_pTemp' 'rDA_iOT' 'rDA_IPS1' 'rDA_IPS2' 'rDA_ParOcc' 'rDA_SMA' 'rDA_sOT' 'rDA2_IPS1' 'rDA2_IPS2' 'rDA2_SPL' 'rFPC_dmPFC' 'rFPC_IFG1' 'rFPC_IFG2' 'rFPC_IFG3' 'rFPC_IFG4' 'rFPC_IPL1' 'rFPC_IPL2' 'rFPC_MFG1' 'rFPC_MFG2' 'rFPC_MFG3' 'rFPC_peanut' 'rFPC_pTemp'}; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    Output   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Y = pdist(data,'correlation');           % 'pdist' converts the square adjacency matrix to a 
                                         %  1 x n matrix so that the function linkage can construct the tree
                                      
clustering = linkage(Y, 'average');      % 'linkage' computes the data to construct the tree
										 % 'average' refers to the UPGMA algorithm


[H,T,perm] = dendrogram(clustering, 0, 'orientation','left','labels', names, 'colorthreshold', .4);         % 'dendrogram' creates the tree

clusters = cluster(clustering, 'Cutoff', 0.5, 'Criterion', 'distance');										% 'cluster' reorders the regions as they 																												  appear on the dendrogram

orient landscape;                                                         							        % orients the dendrogram to 
																											% either landscape (as shown) or portrait

cophenetic_r = cophenet(clustering, Y)   % calculates the cophenetic correlation coefficient 
                                         % which is a measure of the 'goodness of fit' of the data
                                         % to the dendrogram
                                         
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     End    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%