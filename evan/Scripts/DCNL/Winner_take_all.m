contrast = '5'; run = 'Nback'; group = '10_10';

DA = load_nii([run '_DorsAttn_DAT/spmT_000' contrast '.hdr']); DA = DA.img;
DMN = load_nii([run '_Whole_DMN_DAT/spmT_000' contrast '.hdr']); DMN = DMN.img;
FPC = load_nii([run '_Whole_FPC_DAT/spmT_000' contrast '.hdr']); FPC = FPC.img;
Lang = load_nii([run '_Language_DAT/spmT_000' contrast '.hdr']); Lang = Lang.img;
SM = load_nii([run '_SetMaint_DAT/spmT_000' contrast '.hdr']); SM = SM.img;
All(:,:,:,1) = DA; All(:,:,:,2) = DMN; All(:,:,:,3) = FPC; All(:,:,:,4) = Lang; All(:,:,:,5) = SM;
Maxes = max(All,[],4);
Maxes(find(Maxes==0)) = NaN;
DAbest = Template; DAbest(find(Maxes==DA)) = 1;
DMNbest = Template; DMNbest(find(Maxes==DMN)) = 1;
FPCbest = Template; FPCbest(find(Maxes==FPC)) = 1;
Langbest = Template; Langbest(find(Maxes==Lang)) = 1;
SMbest = Template; SMbest(find(Maxes==SM)) = 1;

cd([run '_winners/'])
DAwrite = NIItemplate; DAwrite.img = DAbest; save_nii(DAwrite,[group '_DorsAttn'])
DMNwrite = NIItemplate; DMNwrite.img = DMNbest; save_nii(DMNwrite,[group '_DMN'])
FPCwrite = NIItemplate; FPCwrite.img = FPCbest; save_nii(FPCwrite,[group '_FPC'])
Langwrite = NIItemplate; Langwrite.img = Langbest; save_nii(Langwrite,[group '_Lang'])
SMwrite = NIItemplate; SMwrite.img = SMbest; save_nii(SMwrite,[group '_SM'])
cd ..