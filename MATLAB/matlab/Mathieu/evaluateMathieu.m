%% function [axial,radial] = evaluateMathieu(math,angles,radii)
%
function [axial,radial] = evaluateMathieu(math,angles,radii)
[~,coefficients,indices]=eig_Spm(math.category,math.parameter,math.numberTerms);
if math.category == 1 || math.category == 2,
    trigonometry = 'cosine';
else
    trigonometry = 'sine';
end
axial = series(angles,coefficients,indices,trigonometry); % numberAxial by numberTerms
numberRadial = size(radii,1);
radial = zeros(math.numberTerms,numberRadial);
for i=1:numberRadial,
     radial(:,i)=Jpm(math.category,radii(i),math.parameter,coefficients,math.numberTerms);
end
radial = radial';% numberRadial by numberTerms
