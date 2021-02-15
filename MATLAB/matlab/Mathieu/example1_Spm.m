% example angular Mathieu function at given: q > 0, nmax, and order n
% [( nmax >= max(n) ) & (nmax <= 25)]
% Results are shown in Fig.2
q=1;                      % elliptical parameter 
n=3:5;                    % order n
numberOrders = length(n); % 3
nmax=5;
nCoeffs = getNumberCoefficients();
if ~safeNumberCoefficients(n, nmax, nCoeffs ),
    error('nmax is illegal');
end
vv=0:pi/100:pi/2;         % angles v in radians
numberCoordinates = length(vv); % 51
% Specify the function code k:
for k=1:4
    % characteristic values and expansion coefficients
    [va,mv,vt]=eig_Spm(k,q,nCoeffs);
    my=zeros(numberCoordinates,numberOrders);
    for kv=1:numberCoordinates;
        v=vv(kv);                % take a value of angle v
        vy=Spm(k,v,mv,nmax);     % size [1 nmax]
        my(kv,:)=vy(n);% extract values of Spm corresponding to orders                                                   
    end
    subplot(2,2,k);
    plot(vv/pi,my(:,1),'k.-',vv/pi,my(:,2),'r.-',vv/pi,my(:,3),'b.-'); hold on
    set(gca,'XLim',[0 0.5],'YLim',[-1 1]) 
    xlabel('v/pi')
    ylabel('Spm')
    title([' q=1; n=3:5; KF=',num2str(k)])
end

%
