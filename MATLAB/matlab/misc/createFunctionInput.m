function createFunctionInput( time,value,name,number )
% utility for SierraSD text input files
% writes to a file of given 'name' a SierraSD function of given
% 'number'.  The function is a table of 'time' and 'value' pairs.
% time and value must have the same dimension    
fid = fopen(name,'wt');    
fprintf(fid,'FUNCTION ');
fprintf(fid,'%d\n',number);
fprintf(fid,'  name "');
fprintf(fid,name);
fprintf(fid,'" \n');
fprintf(fid,'  type LINEAR \n');    
fprintf(fid,'  data   %15.10g\t%15.10g\n',[time(:) value(:)]');        
fprintf(fid,'END');    
fclose(fid);

