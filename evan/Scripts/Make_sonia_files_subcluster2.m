load /data/cn4/evan/RestingState/FC_Mapping_120/Subject_clustering_allparcels/Parcel_networkIDs_subjectvariability_UPGMA_maxQ_minsize_20_L_164_alt.mat


%parcels = gifti('/data/cn4/evan/RestingState/FC_Mapping_120/120_L_crossthresh_watershedmerge_164.func.gii'); parcels = parcels.cdata;
%parcelIDs = unique(parcels); parcelIDs(parcelIDs==0) = [];

badparcels = [21440 21790 22092 22502 22718 26419 22260 22481 21837 21780 21824 22111 21459 21477 31421 12 10200 19941 10166 20269 20548 20853 20855 20793 21187 19811 28113 28256 23644 24266 23679];
goodparcelindex = ones(length(parcelIDs),1);
for parcelnum = 1:length(parcelIDs)
    if any(badparcels==parcelIDs(parcelnum));
        goodparcelindex(parcelnum) = 0;
    end
end
parcelIDs = parcelIDs(logical(goodparcelindex));
clusterconnections = clusterconnections(logical(goodparcelindex),:);
clustersizes = clustersizes(logical(goodparcelindex),:);

%parcelnetworks = gifti('/data/cn4/evan/RestingState/FC_Mapping_120/Subject_clustering_allparcels/Parcel_networkIDs_subjectvariability_UPGMA_maxQ_orderedclusters_minsize_20_L_164_alt.func.gii'); parcelnetworks = parcelnetworks.cdata;
%parcelnetworks_alt = gifti('/data/cn4/evan/RestingState/FC_Mapping_120/Subject_clustering_allparcels/Parcel_networkIDs_subjectvariability_UPGMA_maxQ_orderedclusters_minsize_10_L_164.func.gii'); parcelnetworks_alt = parcelnetworks_alt.cdata;

networkcolormap = [1 0 0;0 0 .6;1 1 0;1 .7 .4;0 .8 0;1 .6 1;0 .6 .6;0 0 0;.3 0 .6;.2 1 1;1 .5 0;.6 .2 1;0 .2 .4;.2 1 .2;0 0 1;1 1 .8; 0 .4 0;.25 .25 .25];
networknames = {'DMN','Vis','FPC1','FPC2','DA1','DA2','VA','Sal','CO','SM-H','SM-M','Aud','aMT:','pMTL','MemN','MemR','rPFC',' '};

prmfilename = '/data/cn4/evan/RestingState/FC_Mapping_120/120_L_crossthresh_infomap_goodparcels_forsonia/prmfile2.txt';
clrfilename = '/data/cn4/evan/RestingState/FC_Mapping_120/120_L_crossthresh_infomap_goodparcels_forsonia/watershed_bothhem_Tk0005to01in0005_S1to1_xd20_BI_INFMAP/rawassn_minsize3.txt';
%load('/data/cn4/evan/RestingState/FC_Mapping_120/120_L_crossthresh_infomap_temp/corrmat.mat')

nodecolor = zeros(length(parcelIDs),20,3);
clear color2_names
for parcelnum = 1:length(parcelIDs)
    color1 = clusterconnections(parcelnum,1);
    if color1==0
        color1 = 18;
    end
    
    for colornum = 1:3
        nodecolor(parcelnum,:,colornum) = networkcolormap(color1,colornum);
    end
    
    if clusterconnections(parcelnum,2)
        
        color2_names{parcelnum,1} = [];
        
        for i = 2:nnz(clusterconnections(parcelnum,:))
            
            color2 = clusterconnections(parcelnum,i);
            
            color2_names{parcelnum,1} = [color2_names{parcelnum,1} networknames{color2} '(' num2str(clustersizes(parcelnum,i)) ') '];
            
        end
    else
        color2_names{parcelnum,1} = ' ';
    end
    
end


M_visuals_evan(prmfilename,'thr',clrfilename,1,nodecolor.*255,color2_names,[],20,0,0)








