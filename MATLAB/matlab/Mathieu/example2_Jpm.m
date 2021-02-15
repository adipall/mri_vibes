% example radial Mathieu function Jpm at n, nmax, elliptical parameter q
numberParameters = 5;
parameters=1:numberParameters;
nmax=1;
nCoeffs = getNumberCoefficients();
if ~safeNumberCoefficients(nmax, nmax, nCoeffs ),
    error('nmax is illegal');
end
coordinates=0:5e-02:2.5;
numberCoordinates = length(coordinates); % 51 
for category=1:4
    my = zeros( numberCoordinates, numberParameters);
    for kq=1:numberParameters
        q=parameters(kq);
        % compute characteristic values and expansion coefficients
        [va,mv,vt]=eig_Spm(category,q,nCoeffs);
        y=zeros(numberCoordinates,1);
        for ku=1:numberCoordinates;
            u=coordinates(ku);
            vy=Jpm(category,u,q,mv,nmax);% nmax orders;
            y(ku)=vy(nmax);
        end
        my(:,kq)=y;
    end
    subplot(2,2,category);
    plot(coordinates,my(:,1),'k.-',coordinates,my(:,2),'r.-',...
         coordinates,my(:,3),'b.-',coordinates,my(:,4),'c.-',...
         coordinates,my(:,5),'m.-'); hold on
    set(gca,'XLim',[0 2.5])  
    xlabel('u')
    ylabel('Jpm')
    title([' n=1; q=1:5; KF=',num2str(category)])
end
