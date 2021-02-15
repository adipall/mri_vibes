clear
%% elliptic cylinder functions
%  What is this?  This is a cleaned up version
% of exampleMathieu.  Once this is tested,
% replace exampleMathieu with this.
s=struct('parameter',.1,'numberTerms',8,'category',1);
numberAxial = 11;
numberRadial = 20;
angles = linspace(0,pi/2,numberAxial)';
% 2.03, 2.2, 2.5, 2.7
radii = linspace(0,4,numberRadial)';
[axial,radial] = evaluateMathieu(s,angles,radii);

numberModes = s.numberTerms - 3;
% subplot(1,2,1);
%set(gca,'XLim',[0 1.6])

% plot( axial(:,1:numberModes) );
% hold off; xlabel('\eta');

% subplot(1,2,2);
for i=1:numberModes,
    plot(radii, radial(:,i),'--'); hold on;
end
hold off; xlabel('\xi');
