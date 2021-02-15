clear
[zhangModes, parameter] = first3Modes();
numberModes = 3;
modes = zeros(numberModes*4,1);
last = 0;
for category = 1:4,
    numberCoefficients = 2*numberModes;
    [a,~,~]=eig_Spm(category,parameter,numberCoefficients);
    first = 1 + last;
    last = last + numberModes;
    modes(first:last) = a(1:numberModes);
end
if( norm(modes-zhangModes,inf) < 1.e-7 )
     disp('pass')
else
    disp('fail');
end