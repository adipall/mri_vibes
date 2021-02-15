function [mpcVector,mpcOrder]=readMpc(filename)
% [mpcVector,mpcOrder]=readMpc(filename)
mpcVector=[];
mpcOrder=[];
mpcnumber = 0;
index = 0;
fid=fopen(filename,'r');
finfo = dir(filename);% make sure the file is not empty
fsize = finfo.bytes;
active=false;
if fsize > 0 
    tline = fgetl(fid);
    while ischar(tline)    % !done
       iscomment = strncmp(tline,'//',2);
       if ~iscomment,
           isstart = strncmp(tline,'MPC',3);
           if active && isstart,
               disp('illegal file syntax: start mpc inside an mpc');
           end
           if ~active && isstart,
               active = true;
           end
           isEnd = strncmp(tline,'END',3);
           if ~active && isEnd,
                disp('illegal file syntax: end mpc outside mpc');
           end
           if active && isEnd,
               active = false;
           end
           if active && ~isEnd,               
               [a,count] = fscanf(fid,'%u %s %f\n');
               count = count/3;
               if count > 0
                   mpcnumber = mpcnumber + 1;
                   mpcOrder(mpcnumber,1)=count;
                   node = a(1:3:end);
                   dof = a(2:3:end);
                   dof = dof - 120;
                   weight = a(3:3:end);
                   mpcVector(index+1:index+count,:)=[node,dof,weight];
               end % count > 0
               index = index + count;
           end % active & !end
       end % !comment
       tline = fgetl(fid);
    end % while
else
    disp('warning: empty file');
end
st = fclose(fid);
if st ~= 0, 
    disp('close failed');
end
%a active count dof fid finfo fsize filename index
%iscomment isEnd isstart mpcnumber node st tline weight
