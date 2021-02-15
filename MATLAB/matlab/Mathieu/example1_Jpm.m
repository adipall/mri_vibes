% example evaluate radial Mathieu function Jpm at nmax and orders
% [( nmax >= max(n) ) & (nmax <= 25)]
ellipticalParameter=1;
orders=[3 6 9];
numberOrders = length(orders);
nmax=9;
nCoeffs = getNumberCoefficients();
if ~safeNumberCoefficients(orders, nmax, nCoeffs ),
    error('nmax is illegal');
end
vu=0:5e-02:2.5;            % coordinate u
numberCoordinates = length(vu); % 51
for category=1:4
    % characteristic values expansion coefficients
    [va,mv,vt]=eig_Spm(category,ellipticalParameter,nCoeffs);
    my= zeros(numberCoordinates,numberOrders);
    for ku=1:numberCoordinates;
        u=vu(ku); 
        vy=Jpm(category,u,ellipticalParameter,mv,nmax);
        my(ku,:)=vy(orders);% extract values at orders n=[3 6 9]
    end
    subplot(2,2,category);
    plot(vu,my(:,1),'k.-',vu,my(:,2),'r.-',vu,my(:,3),'b.-'); hold on
    set(gca,'XLim',[0 2.5]) 
    xlabel('u')
    ylabel('Jpm')
    title([' ellipticalParameter 1; orders=3 6 9; KF=',num2str(category)])
end
