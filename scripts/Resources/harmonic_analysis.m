function harmonic_analysis(datain)

data = load_untouch_nii_2D(datain);
img_out = log(abs(fftshift(fft(double(data.img),[],2),2)));
data.img = img_out;

save_untouch_nii_2D(data,[datain(1:end-7) '_harmonic.nii.gz'])