%function [Ke]=soli8Kf(ex,ey,ez,ep,fl_rho,c)
%
%-------------------------------------------------------------
% PURPOSE
%  Calculate the stiffness matrix for a 8 node (brick)
%  fluid element.
%
% INPUT:   ex = [x1 x2 x3 ... x8]
%          ey = [y1 y2 y3 ... y8]  element coordinates
%          ez = [z1 z2 z3 ... z8]
%
%          ep = [ir]               ir integration rule
%
%          fl_rho: density of fluid 
%          c: speed of sound of fluid
%
% OUTPUT: Ke : element stiffness matrix
%         fe : equivalent nodal forces 
%
% LAST MODIFIED: M Ross   16 May 2005
%-------------------------------------------------------------

function [Ke]=soli8Kf(ex,ey,ez,ep,fl_rho,c)

%-------------------------------------------------------------

  ir=ep(1);  ngp=ir*ir*ir;

%--------- gauss points --------------------------------------
%   if ir==1
%     g1=0.0; w1=2.0;
%     %gp=[ g1 g1 ];  w=[ w1 w1 ];
%     gp=[ g1 g1 g1 ];  w=[ w1 w1 w1 ];
  if ir==2
    g1=0.577350269189626; w1=1;
    gp(:,1)=[-1; 1; 1;-1;-1; 1; 1;-1]*g1; w(:,1)=[ 1; 1; 1; 1; 1; 1; 1; 1]*w1;
    gp(:,2)=[-1;-1; 1; 1;-1;-1; 1; 1]*g1; w(:,2)=[ 1; 1; 1; 1; 1; 1; 1; 1]*w1;
    gp(:,3)=[-1;-1;-1;-1; 1; 1; 1; 1]*g1; w(:,3)=[ 1; 1; 1; 1; 1; 1; 1; 1]*w1;
  elseif ir==3
    g1=0.774596669241483; g2=0.;
    w1=0.555555555555555; w2=0.888888888888888;

    I1=[-1; 0; 1;-1; 0; 1;-1; 0; 1]'; I2=[ 0;-1; 0; 0; 1; 0; 0; 1; 0]';
    gp(:,1)=[I1 I1 I1]'*g1;           gp(:,1)=[I2 I2 I2]'*g2+gp(:,1);
    I1=abs(I1); I2=abs(I2);
    w(:,1)=[I1 I1 I1]'*w1;
    w(:,1)=[I2 I2 I2]'*w2+w(:,1);
    I1=[-1;-1;-1; 0; 0; 0; 1; 1; 1]'; I2=[ 0; 0; 0; 1; 1; 1; 0; 0; 0]';
    gp(:,2)=[I1 I1 I1]'*g1;           gp(:,2)=[I2 I2 I2]'*g2+gp(:,2);
    I1=abs(I1); I2=abs(I2);
    w(:,2)=[I1 I1 I1]'*w1;
    w(:,2)=[I2 I2 I2]'*w2+w(:,2);
    I1=[-1;-1;-1;-1;-1;-1;-1;-1;-1]'; I2=[ 0; 0; 0; 0; 0; 0; 0; 0; 0]';
    I3=abs(I1);
    gp(:,3)=[I1 I2 I3]'*g1;           gp(:,3)=[I2 I3 I2]'*g2+gp(:,3);
    w(:,3)=[I3 I2 I3]'*w1;
    w(:,3)=[I2 I3 I2]'*w2+w(:,3);
  else
    disp('Used number of integration points not implemented');
    return
  end;
  
%   if ir==1
%       wp=w(:,1).*w(:,2);
%   else
%       wp=w(:,1).*w(:,2).*w(:,3);
%   end
  wp=w(:,1).*w(:,2).*w(:,3);
  xsi=gp(:,1);  eta=gp(:,2); zet=gp(:,3);  r2=ngp*3;

%--------- shape functions -----------------------------------
  N(:,1)=(1-xsi).*(1-eta).*(1-zet)/8;  N(:,5)=(1-xsi).*(1-eta).*(1+zet)/8;
  N(:,2)=(1+xsi).*(1-eta).*(1-zet)/8;  N(:,6)=(1+xsi).*(1-eta).*(1+zet)/8;
  N(:,3)=(1+xsi).*(1+eta).*(1-zet)/8;  N(:,7)=(1+xsi).*(1+eta).*(1+zet)/8;
  N(:,4)=(1-xsi).*(1+eta).*(1-zet)/8;  N(:,8)=(1-xsi).*(1+eta).*(1+zet)/8;

  dNr(1:3:r2,1)=-(1-eta).*(1-zet);    dNr(1:3:r2,2)= (1-eta).*(1-zet);
  dNr(1:3:r2,3)= (1+eta).*(1-zet);    dNr(1:3:r2,4)=-(1+eta).*(1-zet);
  dNr(1:3:r2,5)=-(1-eta).*(1+zet);    dNr(1:3:r2,6)= (1-eta).*(1+zet);
  dNr(1:3:r2,7)= (1+eta).*(1+zet);    dNr(1:3:r2,8)=-(1+eta).*(1+zet);
  dNr(2:3:r2+1,1)=-(1-xsi).*(1-zet);  dNr(2:3:r2+1,2)=-(1+xsi).*(1-zet);
  dNr(2:3:r2+1,3)= (1+xsi).*(1-zet);  dNr(2:3:r2+1,4)= (1-xsi).*(1-zet);
  dNr(2:3:r2+1,5)=-(1-xsi).*(1+zet);  dNr(2:3:r2+1,6)=-(1+xsi).*(1+zet);
  dNr(2:3:r2+1,7)= (1+xsi).*(1+zet);  dNr(2:3:r2+1,8)= (1-xsi).*(1+zet);
  dNr(3:3:r2+2,1)=-(1-xsi).*(1-eta);  dNr(3:3:r2+2,2)=-(1+xsi).*(1-eta);
  dNr(3:3:r2+2,3)=-(1+xsi).*(1+eta);  dNr(3:3:r2+2,4)=-(1-xsi).*(1+eta);
  dNr(3:3:r2+2,5)= (1-xsi).*(1-eta);  dNr(3:3:r2+2,6)= (1+xsi).*(1-eta);
  dNr(3:3:r2+2,7)= (1+xsi).*(1+eta);  dNr(3:3:r2+2,8)= (1-xsi).*(1+eta);
  dNr=dNr/8.;


  Ke=zeros(8,8);  
  JT=dNr*[ex;ey;ez]';

%--------- three dimensional case ----------------------------
  for i=1:ngp
    indx=[ 3*i-2; 3*i-1; 3*i ];
    detJ=det(JT(indx,:));
    if detJ<10*eps
      disp('Jacobideterminant equal or less than zero!')
    end
    JTinv=inv(JT(indx,:));
    dNx=JTinv*dNr(indx,:);

    B(1,1:3:24-2)=dNx(1,:);
    B(1,2:3:24-1)=dNx(2,:);
    B(1,3:3:24)  =dNx(3,:);
%     B(4,1:3:24-2)=dNx(2,:);
%     B(4,2:3:24-1)=dNx(1,:);
%     B(5,1:3:24-2)=dNx(3,:);
%     B(5,3:3:24)  =dNx(1,:);
%     B(6,2:3:24-1)=dNx(3,:);
%     B(6,3:3:24)  =dNx(2,:);


    % Ke=Ke+fl_rho*c^2*B'*B*detJ*wp(i);
    Ke=Ke+(dNx(1,:)'*dNx(1,:)+dNx(2,:)'*dNx(2,:)+dNx(3,:)'*dNx(3,:))*detJ*wp(i);

  end
%--------------------------end--------------------------------
