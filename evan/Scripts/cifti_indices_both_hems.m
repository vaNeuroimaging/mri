%%
surf_L = 1:29696;
surf_R = 29697:59412;
acc_L = 59413:59436;
acc_R = 59437:59471;
amyg_L = 59472:59555;
amyg_R = 59556:59631;
caud_L = 59632:59812;
caud_R = 59813:59982;
cbll_L = 59983:62451;
cbll_R = 62452:64949;
hipp_L = 64950:65165;
hipp_R = 65166:65371;
pall_L = 65372:65452;
pall_R = 65453:65518;
puta_L = 65519:65781;
puta_R = 65782:66031;
thal_L = 66032:66361;
thal_R = 66362:66697;

%%

cifti_coords_L_contra = cifti_coords([surf_L acc_L amyg_L caud_L cbll_R hipp_L pall_L puta_L thal_L],:);
dlmwrite(['/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/modified_cifti_network/LEFT_contra_cbllR/cifti_coords.txt'],cifti_coords_L_contra,'delimiter',' ')
quickroifile(cifti_coords_L_contra,['/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/modified_cifti_network/LEFT_contra_cbllR/cifti_coords.roi'])

cifti_coords_R_contra = cifti_coords([surf_R acc_R amyg_R caud_R cbll_L hipp_R pall_R puta_R thal_R],:);
dlmwrite(['/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/modified_cifti_network/RIGHT_contra_cbllL/cifti_coords.txt'],cifti_coords_R_contra,'delimiter',' ')
quickroifile(cifti_coords_R_contra,['/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/modified_cifti_network/RIGHT_contra_cbllL/cifti_coords.roi'])


%%


avgcrosscorr_L_contra = avgcrosscorr.avgcrosscorr([surf_L acc_L amyg_L caud_L cbll_R hipp_L pall_L puta_L thal_L],[surf_L acc_L amyg_L caud_L cbll_R hipp_L pall_L puta_L thal_L]);
save(['/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/modified_cifti_network/LEFT_contra_cbllR/cifti_avgcrosscorr.mat'],'avgcrosscorr_L_contra','-v7.3')

avgcrosscorr_R_contra = avgcrosscorr.avgcrosscorr([surf_R acc_R amyg_R caud_R cbll_L hipp_R pall_R puta_R thal_R],[surf_R acc_R amyg_R caud_R cbll_L hipp_R pall_R puta_R thal_R]);
save(['/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/modified_cifti_network/RIGHT_contra_cbllL/cifti_avgcrosscorr.mat'],'avgcrosscorr_R_contra','-v7.3')

%%

HEMS = {'L';'R'};
hemname = {'LEFT';'RIGHT'};
cd('/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/modified_cifti_network/LEFT_contra_cbllR')
graphcluster('prmfile.txt','thr',20,1,0,'infomap');
cd('/data/cn4/laumannt/fcMapping_redux/FCPROCESS_NEW_042213/modified_cifti_network/RIGHT_contra_cbllL')
graphcluster('prmfile.txt','thr',20,1,0,'infomap');
