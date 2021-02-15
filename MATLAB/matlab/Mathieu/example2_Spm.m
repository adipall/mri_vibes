% example angular Mathieu function Spm at order n, nmax, q
% [( nmax >= n )& (nmax <= 25)]
% Results are shown in Fig.3
vq=1:5;                % 5 values of elliptical parameter q
n=1;                   % single value of order n
nmax=1;
nCoeffs = getNumberCoefficients();
if ~safeNumberCoefficients(n, nmax, nCoeffs ),
    error('nmax is illegal');
end
vv=0:pi/20:2*pi;       % 41 values of angle v in radians 
for k=1:4  % Specify the function code k:
    my = zeros(41,5);
    for kq=1:5  
        q=vq(kq);             % take a value of q 
        % compute characteristic values and expansion coefficients 
        [va,mv,vt]=eig_Spm(k,q,nCoeffs);
        y = zeros(41,1);
        for kv=1:length(vv);
            v=vv(kv);
            vy=Spm(k,v,mv,nmax);% Spm for nmax orders;
            yn=vy(n); % extract order n
            y(kv) =yn;  % column vector of Spm values at different
                        % angles v, at order n, and a value q;
        end
        my(:,k)=y; % matrix of Spm v and q;
    end
    subplot(2,2,k);
    plot(vv/pi,my(:,1),'k.-',vv/pi,my(:,2),'r.-',vv/pi,my(:,3),'b.-', ...
         vv/pi,my(:,4),'c.-',vv/pi,my(:,5),'m.-'); hold on
    xlabel('v/pi')
    ylabel('Spm')
    title([' n=1; q=1:5; KF=',num2str(k)])
end
