clear all
%%Load Subject IDs
Joe_a = fopen('Joe_subj.txt');
Maital_a = fopen('Maital_subj.txt');
Joe_b = textscan(Joe_a,'%s');
Joe = Joe_b{1,1};
Maital_b = textscan(Maital_a,'%s');
Maital = Maital_b{1,1};
fclose (Joe_a);
fclose (Maital_a);
clear Joe_a Joe_b Maital_a Maital_b

Joe{13,1} = 'vc34307';
Maital{13,1} = 'vc33543';
%n=1;
%x=13;
x = length(Joe);
for n=1:x
    n=13
if n<13
disp ('start');
disp (Maital{n,1});
%% load matrices

cd ACRN/subjects/

AC_Errors = read_4dfpimg(strcat(Maital{n,1},'_10runs_Error10frames_avgtc_AC_Errors_1_2_3_4_5_6_7_8_9_10_333_t88.4dfp.img'));
AC_Ambiguous = read_4dfpimg(strcat(Maital{n,1},'_10runs_Error10frames_avgtc_AC_Ambiguous_1_2_3_4_5_6_7_8_9_10_333_t88.4dfp.img'));
AC_Clear = read_4dfpimg(strcat(Maital{n,1},'_10runs_Error10frames_avgtc_Abstract+Concrete_1_2_3_4_5_6_7_8_9_10_333_t88.4dfp.img'));

RN_Errors = read_4dfpimg(strcat(Maital{n,1},'_10runs_Error10frames_avgtc_RN_Errors_1_2_3_4_5_6_7_8_9_10_333_t88.4dfp.img'));
RN_Ambiguous = read_4dfpimg(strcat(Maital{n,1},'_10runs_Error10frames_avgtc_RN_Ambiguous_1_2_3_4_5_6_7_8_9_10_333_t88.4dfp.img'));
RN_Clear = read_4dfpimg(strcat(Maital{n,1},'_10runs_Error10frames_avgtc_Rhyme+NoRhyme_1_2_3_4_5_6_7_8_9_10_333_t88.4dfp.img'));

cd ../../

cd MRTNV/subjects/

GP_Correct = read_4dfpimg(strcat(Joe{n,1},'_Scrubbed_f06_NameShortened_X2_TrialsExtended_SustainedDelayed_avgtc_Glass_0_Coh+Glass_0.125_Coh+Glass_0.25_Coh+Glass_0.5_Coh_1_2_3_4_5_6_7_8_9_10_333_t88.4dfp.img'));
GP_Errors = read_4dfpimg(strcat(Joe{n,1},'_Scrubbed_f06_NameShortened_X2_TrialsExtended_SustainedDelayed_avgtc_Glass_Incorrect_1_2_3_4_5_6_7_8_9_10_333_t88.4dfp.img'));

MR_Errors = read_4dfpimg(strcat(Joe{n,1},'_Scrubbed_f06_NameShortened_X2_TrialsExtended_SustainedDelayed_avgtc_MRT_Incorrect_1_2_3_4_5_6_7_8_9_10_333_t88.4dfp.img'));
MR_Correct = read_4dfpimg(strcat(Joe{n,1},'_Scrubbed_f06_NameShortened_X2_TrialsExtended_SustainedDelayed_avgtc_MRT_Same_Small+MRT_Same_Mid+MRT_Same_Large+MRT_Mirror_Small+MRT_Mirror_Mid+MRT_Mirror_Large_1_2_3_4_5_6_7_8_9_10_333_t88.4dfp.img'));

NV_Errors = read_4dfpimg(strcat(Joe{n,1},'_Scrubbed_f06_NameShortened_X2_TrialsExtended_SustainedDelayed_avgtc_Verbal_Incorrect_1_2_3_4_5_6_7_8_9_10_333_t88.4dfp.img'));
NV_Correct = read_4dfpimg(strcat(Joe{n,1},'_Scrubbed_f06_NameShortened_X2_TrialsExtended_SustainedDelayed_avgtc_Verbal_Noun+Verbal_Verb_1_2_3_4_5_6_7_8_9_10_333_t88.4dfp.img'));

%%GP_Correct = read_4dfpimg(strcat(Joe{n,1},'_Scrubbed_f06_NameShortened_X2_TrialsExtended_SustainedDelayed_Glass_0_Coh+Glass_0.125_Coh+Glass_0.25_Coh+Glass_0.5_Coh_1_2_3_4_5_6_7_8_9_10_333_t88.4dfp.img'));
%%GP_Errors = read_4dfpimg(strcat(Joe{n,1},'_Scrubbed_f06_NameShortened_X2_TrialsExtended_SustainedDelayed_Glass_Incorrect_1_2_3_4_5_6_7_8_9_10_333_t88.4dfp.img'));

cd ../../

test_image=zeros(size(AC_Clear,1),1);


%% make corrmat from new averaged timepoint matrices
Conc = zeros(length(test_image),120);
Conc(:,1:10) = AC_Errors;
Conc(:,11:20) = AC_Ambiguous;
Conc(:,21:30) = AC_Clear;
Conc(:,31:40) = RN_Errors;
Conc(:,41:50) = RN_Ambiguous;
Conc(:,51:60) = RN_Clear;
Conc(:,61:70) = GP_Errors;
Conc(:,71:80) = GP_Correct;
Conc(:,81:90) = MR_Errors;
Conc(:,91:100) = MR_Correct;
Conc(:,101:110) = NV_Errors;
Conc(:,111:120) = NV_Correct;

%% put into volume
cd MergedSubjects/
file_out = strcat('Subject_',Maital{n,1},'_All','.4dfp.img');
write_4dfpimg(Conc,file_out,'bigendian');
write_4dfpifh(file_out ,120,'bigendian');

cd 2Tasks/
Short = zeros(length(test_image),40);
Short(:,1:10) = AC_Errors;
Short(:,11:20) = AC_Clear;
Short(:,21:30) = GP_Errors;
Short(:,31:40) = GP_Correct;

file_out2 = strcat('Subject_',Maital{n,1},'_AC_GP','.4dfp.img');
write_4dfpimg(Short,file_out2,'bigendian');
write_4dfpifh(file_out2 ,120,'bigendian');
cd ../../

disp ('done');
disp(n);


else
cd ACRN/subjects/

AC_Errors = read_4dfpimg(strcat(Maital{n,1},'_10runs_Error10frames_avgtc_AC_Errors_1_2_3_4_5_6_7_8_9_10_333_t88.4dfp.img'));
AC_Ambiguous = read_4dfpimg(strcat(Maital{n,1},'_10runs_Error10frames_avgtc_AC_Ambiguous_1_2_3_4_5_6_7_8_9_10_333_t88.4dfp.img'));
AC_Clear = read_4dfpimg(strcat(Maital{n,1},'_10runs_Error10frames_avgtc_Abstract+Concrete_1_2_3_4_5_6_7_8_9_10_333_t88.4dfp.img'));

cd ../../

cd MRTNV/subjects/

GP_Correct = read_4dfpimg(strcat(Joe{n,1},'_Scrubbed_f06_NameShortened_X2_TrialsExtended_SustainedDelayed_Glass_0_Coh+Glass_0.125_Coh+Glass_0.25_Coh+Glass_0.5_Coh_1_2_3_4_5_6_7_8_9_10_333_t88.4dfp.img'));
GP_Errors = read_4dfpimg(strcat(Joe{n,1},'_Scrubbed_f06_NameShortened_X2_TrialsExtended_SustainedDelayed_Glass_Incorrect_1_2_3_4_5_6_7_8_9_10_333_t88.4dfp.img'));

cd ../../

test_image=zeros(size(AC_Clear,1),1);

%% put into volume
cd MergedSubjects/
cd 2Tasks/
Short = zeros(length(test_image),40);
Short(:,1:10) = AC_Errors;
Short(:,11:20) = AC_Clear;
Short(:,21:30) = GP_Errors;
Short(:,31:40) = GP_Correct;

file_out2 = strcat('Subject_',Maital{n,1},'_AC_GP','.4dfp.img');
write_4dfpimg(Short,file_out2,'bigendian');
write_4dfpifh(file_out2 ,120,'bigendian');
cd ../../

disp ('done');
disp(n);  


end
end







