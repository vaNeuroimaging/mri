list = '/data/cn3/steven/NP847/glmsLists/28subs_16tps_Ret_list_of_glms.glm_list';
%'/data/cn3/steven/NP847/glmsLists/28subs_16tps_Prim_list_of_glms.glm_list';

outfolder = '/data/cn4/evan/Task_parcellation/SteveSubjects';
oldpath = '\/data\/nil-external\/mcl\/Nelson\/';
newpath = '\/data\/cn3\/steven\/';

[front entries] = textread(list,'%s%s','delimiter',':');

numsubjects = str2num(entries{2});
glmfiles = entries(3:2+numsubjects);
T4files = entries(4+numsubjects:end);

for glmnum = 1:length(glmfiles)
    
    filename = dir(glmfiles{glmnum}); filename = filename(1).name;
    copyfile(glmfiles{glmnum},[outfolder '/temp_' filename])
    
    evalc(['!sed ''s/' oldpath '/' newpath '/g'' temp_' filename '>' filename]);
    
    delete([outfolder '/temp_' filename])
    
%     fid = fopen([outfolder '/' filename],'r+');
%     header = textscan(fid,'%s',132,'Delimiter','%\n');
%     header = header{1};
%     fclose(fid);
%     
%     for headerline = 1:length(header)
%         thisheader = header{headerline};
%         pathloc = strfind(thisheader,oldpath);
%         
%         if ~isempty(pathloc)
%             
%             newheaderline = [thisheader(1:pathloc-1) newpath thisheader(pathloc+length(oldpath) : end) repmat(' ',1,(length(oldpath)-length(newpath)))];
%             
%             header{headerline} = newheaderline;
%         end
%     end
%     
%     fid = fopen([outfolder '/' filename],'r+');
%     frewind(fid)
%     for headerline = 1:length(header)
% %       if headerline<length(header)
%             fprintf(fid, '%s\n', header{headerline});
% %         else
% %             fprintf(fid, '%s\n', header{headerline});
% %         end
%     end
%     fclose(fid)
end
