
function wrap_NP1(isub,caser)

for i=isub
    
    switch i
        case 1
            pdir=['/misc/data75/powerjd/ABIDE'];
            figdir=[pdir '/FIGS'];
            mat=[figdir '/P1.mat'];
            s=[1 187];
            callfn(pdir,'JDP_P1_preprocess.sh',s(1),s(2),caser,mat);
        case 2
            pdir=['/misc/data48/powerjd/NIH'];
            figdir=[pdir '/FIGS'];
            mat=[figdir '/P1.mat'];
            s=[1 91];
            callfn(pdir,'JDP_P1_preprocess.sh',s(1),s(2),caser,mat);
        case 3
            pdir=['/misc/data48/powerjd/WASHU'];
            figdir=[pdir '/FIGS'];
            mat=[figdir '/P1.mat'];
            s=[1 120];
            callfn(pdir,'JDP_P1_preprocess.sh',s(1),s(2),caser,mat);
        case 4
            pdir=['/misc/data48/powerjd/RANDY'];
            figdir=[pdir '/FIGS'];
            mat=[figdir '/P1.mat'];
            s=[1 231];
            callfn(pdir,'JDP_P1_preprocess.sh',s(1),s(2),caser,mat);
        case 5
            pdir=['/misc/data75/powerjd/ME'];
            figdir=[pdir '/FIGS'];
            mat=[figdir '/P1.mat'];
            s=[1 89];
            callfn(pdir,'JDP_P1_preprocess.sh',s(1),s(2),caser,mat);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function callfn(v1,v2,v3,v4,v5,mat)

load(mat);
JDP_P1_gather(v1,v2,v3,v4,v5,QC);
clear;


