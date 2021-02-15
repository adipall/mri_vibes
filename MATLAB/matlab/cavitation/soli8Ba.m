% function [Bele] = soli8Ba(ex,ey,ez,ep,pvapor,dvel_pot,penalty)
%
%--------------------------------------------------------------------------
% This is a function to calculate the penalty matrix of an acoustic brick element. 
%
% Input:   ex = [x1 x2 x3 ... x8]
%          ey = [y1 y2 y3 ... y8]  element coordinates
%          ez = [z1 z2 z3 ... z8]
%
%          ep = [ir]               ir integration rule
%         
%          rho = mass density
%
%          c = speed of sound
%
%          penalty = penalty parameter
%
% Ouput: The Element Mass Matrix (Me)
%
% Created: M. Ross 16May05
% Modified: M. Ross 1/3/2017
%
%--------------------------------------------------------------------------

function [Bele] = soli8Ba(ex,ey,ez,ep,pvapor,dvel_pot,penalty)

%--------------------------------------------------------------------------

ir=ep(1);  ngp=ir*ir*ir;

%--------- gauss points --------------------------------------
%   if ir==1
%     g1=0.0; w1=2.0;
%     gp=[ g1 g1 ];  w=[ w1 w1 ];
  if ir==2
    g1=0.577350269189626; w1=1;
    gp(:,1)=[-1; 1; 1;-1;-1; 1; 1;-1]*g1; 
    w(:,1)=[ 1; 1; 1; 1; 1; 1; 1; 1]*w1;
    gp(:,2)=[-1;-1; 1; 1;-1;-1; 1; 1]*g1; 
    w(:,2)=[ 1; 1; 1; 1; 1; 1; 1; 1]*w1;
    gp(:,3)=[-1;-1;-1;-1; 1; 1; 1; 1]*g1; 
    w(:,3)=[ 1; 1; 1; 1; 1; 1; 1; 1]*w1;
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
  dNr(2:3:r2,1)=-(1-xsi).*(1-zet);  dNr(2:3:r2,2)=-(1+xsi).*(1-zet);%1
  dNr(2:3:r2,3)= (1+xsi).*(1-zet);  dNr(2:3:r2,4)= (1-xsi).*(1-zet);
  dNr(2:3:r2,5)=-(1-xsi).*(1+zet);  dNr(2:3:r2,6)=-(1+xsi).*(1+zet);
  dNr(2:3:r2,7)= (1+xsi).*(1+zet);  dNr(2:3:r2,8)= (1-xsi).*(1+zet);
  dNr(3:3:r2,1)=-(1-xsi).*(1-eta);  dNr(3:3:r2,2)=-(1+xsi).*(1-eta);%2
  dNr(3:3:r2,3)=-(1+xsi).*(1+eta);  dNr(3:3:r2,4)=-(1-xsi).*(1+eta);
  dNr(3:3:r2,5)= (1-xsi).*(1-eta);  dNr(3:3:r2,6)= (1+xsi).*(1-eta);
  dNr(3:3:r2,7)= (1+xsi).*(1+eta);  dNr(3:3:r2,8)= (1-xsi).*(1+eta);
  dNr=dNr/8.0;

  Bele = zeros(8,8);
  Bz = zeros(8,8);
  JT=dNr*[ex;ey;ez]';
  %NtN = zeros(8,8);
  
  for i=1:ngp
    pg = N(i,:)*dvel_pot';
    if (pg<pvapor)
        
        indx=[ 3*i-2; 3*i-1; 3*i ];
        detJ=det(JT(indx,:));
        if detJ<10*eps
            disp('Jacobideterminant equal or less than zero!')
        end
        
        Bele = Bele + N(i,:)'*N(i,:)*detJ*wp(i)*(penalty);
        
    end
        
    %NtN = NtN + N'*N*detJ*wp(i)*rho;
  end