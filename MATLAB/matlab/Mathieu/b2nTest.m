% Gutierrez-Vega, Table 4.1. 
% Example 1  b_{2n+2} modes of angular Mathieu functions
clear
parity = 0;
category = 3;
[modes, parameter] = bEvenModes();
number = 16; 
[values,mv,vt]=eig_Spm(category,parameter,number);
if( norm(values(1:8)- modes,inf) < 1.e-9 )
     disp('pass')
else
    disp('fail');
end