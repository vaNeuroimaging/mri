function FSL_from_matlab(fsl_command)
% FSL_from_matlab(fsl_command)
% 
% Takes any string and executes it as a command on the terminal.  Any
% terminal output text is printed onto the Matlab command screen.
%
% Technique developed by X. You; script created by E. Gordon, Jan 2011

commandfilename = 'FSLCommand';
delete(commandfilename);
fid = fopen(commandfilename,'at');
fprintf(fid,fsl_command);
fclose(fid);
eval(['!chmod a+rwx ' commandfilename]);
eval(['!' pwd '/' commandfilename]);
delete(commandfilename);