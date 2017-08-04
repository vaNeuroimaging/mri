delete outliers.log
system('R --no-save --slave < outliers.R > outliers.log')

addpath /home/data/Analysis/ENIGMA/ENIGMA_QC_scripts/
subjects = '/home/data/subjects/processing_list_060717.txt';
if iscell(subjects)

elseif exist(subjects,'file')
    
    subjects = textread(subjects,'%s');
    
elseif ischar(subjects)
    
    subjects = {subjects};
    
end
for s = 1:length(subjects)
    func_make_corticalpngs_ENIGMA_QC('/home/data/Analysis/ENIGMA/QC/',subjects{s},'/mri/orig.mgz','/mri/aparc+aseg.mgz')
end
%system('/home/data/Analysis/ENIGMA/ENIGMA_QC_scripts/make_ENIGMA_QC_webpage.sh QC/')

%%
% cd QC
% system('echo "<html>" 												>> ENIGMA_Cortical_QC.html; echo "<head>"                                                   >> ENIGMA_Cortical_QC.html')
% system('echo "<style type=\text/css\>"								>> ENIGMA_Cortical_QC.html')
% system('echo "*"                                                        >> ENIGMA_Cortical_QC.html; echo "{"														>> ENIGMA_Cortical_QC.html; echo "margin: 0px;"												>> ENIGMA_Cortical_QC.html; echo "padding: 0px;"											>> ENIGMA_Cortical_QC.html; echo "}"														>> ENIGMA_Cortical_QC.html; echo "html, body"												>> ENIGMA_Cortical_QC.html; echo "{"														>> ENIGMA_Cortical_QC.html; echo "height: 100%;"											>> ENIGMA_Cortical_QC.html; echo "}"														>> ENIGMA_Cortical_QC.html; echo "</style>"													>> ENIGMA_Cortical_QC.html; echo "</head>"													>> ENIGMA_Cortical_QC.html; echo "<body>" 													>> ENIGMA_Cortical_QC.html;')
% for s = 1:length(subjects)
%     
%     system(['sub=' subjects(s).name '; echo "<table cellspacing=\"1\" style=\"width:100%;background-color:#000;\">"				>> ENIGMA_Cortical_QC.html; echo "<tr>"																					>> ENIGMA_Cortical_QC.html; echo "<td> <FONT COLOR=WHITE FACE=\"Geneva, Arial\" SIZE=5> $sub </FONT> </td>"             >> ENIGMA_Cortical_QC.html; echo "</tr>"                                                                                >> ENIGMA_Cortical_QC.html; echo "<tr>"                                                                                 >> ENIGMA_Cortical_QC.html; echo "<td><a href=\"file:"$sub"/Cortical_set_Axial_20.png\"><img src=\""$sub"/Cortical_set_Axial_70_20.png\" width=\"100%\" ></a></td>"	>> ENIGMA_Cortical_QC.html; echo "<td><a href=\"file:"$sub"/Cortical_set_Axial_40.png\"><img src=\""$sub"/Cortical_set_Axial_70_40.png\" width=\"100%\" ></a></td>"	>> ENIGMA_Cortical_QC.html; echo "<td><a href=\"file:"$sub"/Cortical_set_Axial_50.png\"><img src=\""$sub"/Cortical_set_Axial_70_50.png\" width=\"100%\" ></a></td>"	>> ENIGMA_Cortical_QC.html; echo "<td><a href=\"file:"$sub"/Cortical_set_Axial_60.png\"><img src=\""$sub"/Cortical_set_Axial_70_60.png\" width=\"100%\" ></a></td>"	>> ENIGMA_Cortical_QC.html; echo "<td><a href=\"file:"$sub"/Cortical_set_Axial_80.png\"><img src=\""$sub"/Cortical_set_Axial_70_80.png\" width=\"100%\" ></a></td>"	>> ENIGMA_Cortical_QC.html; echo "</tr>"																				>> ENIGMA_Cortical_QC.html; echo "<tr>"																					>> ENIGMA_Cortical_QC.html; echo "<td><a href=\"file:"$sub"/Cortical_set_Coronal_20.png\"><img src=\""$sub"/Cortical_set_Coronal_70_20.png\" width=\"100%\" ></a></td>"	>> ENIGMA_Cortical_QC.html; echo "<td><a href=\"file:"$sub"/Cortical_set_Coronal_40.png\"><img src=\""$sub"/Cortical_set_Coronal_70_40.png\" width=\"100%\" ></a></td>"	>> ENIGMA_Cortical_QC.html; echo "<td><a href=\"file:"$sub"/Cortical_set_Coronal_50.png\"><img src=\""$sub"/Cortical_set_Coronal_70_50.png\" width=\"100%\" ></a></td>"	>> ENIGMA_Cortical_QC.html; echo "<td><a href=\"file:"$sub"/Cortical_set_Coronal_60.png\"><img src=\""$sub"/Cortical_set_Coronal_70_60.png\" width=\"100%\" ></a></td>"	>> ENIGMA_Cortical_QC.html; echo "<td><a href=\"file:"$sub"/Cortical_set_Coronal_80.png\"><img src=\""$sub"/Cortical_set_Coronal_70_80.png\" width=\"100%\" ></a></td>"	>> ENIGMA_Cortical_QC.html; echo "</tr>"																				>> ENIGMA_Cortical_QC.html; echo "</table>"                                                                             >> ENIGMA_Cortical_QC.html;'])
% end
% system('echo "</body>"                                                                              >> ENIGMA_Cortical_QC.html; echo "</html>"                                                                              >> ENIGMA_Cortical_QC.html')

