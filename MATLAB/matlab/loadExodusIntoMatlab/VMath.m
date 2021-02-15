classdef VMath
    properties
    end
    methods (Static)
        %% General, up to 3x3 [A11 A12 A13 A21 A22 A23 A31 A32 A33]
        function [z,detA,A]=Solve33(A,z,dim)
            %
            %  A is stored as A=(A11 A12 A13 A21 A22 A23 A31 A32 A33)
            %          if dim==2, A is stored the same way but only
            %          A11 A12, A21 and A22 are used.
            %
            %  z is a vector stored as z(:,:,i)=(z1 z2 z3)
            %
            % We'll use gauss elimination with partial pivioting to solve
            %
            %  the hold(idx) lines are slowing it down.
            if nargin>2 && dim<3,
                switch dim
                    case 2
                        detA=A(:,1).*A(:,5)-A(:,2).*A(:,4);
                        A(:,1:4)=[A(:,5) -A(:,2) -A(:,4) A(:,1)]./detA(:,[1 1 1 1]);
                        z=[A(:,1).*z(:,1)+A(:,2).*z(:,2) A(:,3).*z(:,1)+A(:,4).*z(:,2) zeros(size(z,1),1)];
                    case 1
                        detA=A(:,1);
                        z=[z(:,1)./A(:,1) zeros(size(z,1),2)];
                end
                return
            end

            [n,m,nlhs]=size(z);
            
            nrows=size(A,1);
            [ju,midx1]=max(abs(A(:,[1 4 7])),[],2);
            pividx=[1 2 3 4 5 6 7 8 9
                4 5 6 1 2 3 7 8 9
                7 8 9 4 5 6 1 2 3];
            idx=pividx(midx1,:)*nrows-repmat(((nrows-1):-1:0)',1,9);

            A=A(idx);

            pividx=[1 2 3
                2 1 3
                3 2 1];
            sgn=[1;-1;-1];
            sign=sgn(midx1);
            idx=pividx(midx1,:)*nrows-repmat(((nrows-1):-1:0)',1,3);
            
            %Pivioting moved to next for loop
            %             for i=1:nlhs
            %                 hold=z(:,:,i);
            %                 z(:,:,i)=hold(idx);
            %             end
            
            if any(~A(:,1)),
                error('Some matrices are singular');
            end
            const=A(:,4)./A(:,1);
            A(:,[5 6])=A(:,[5 6])-const(:,[1 1]).*A(:,[2 3]);
            
            
            for i=1:nlhs
                hold=z(:,:,i);
                z(:,:,i)=hold(idx);
                z(:,2,i)=z(:,2,i)-const.*z(:,1,i);
            end
            
            
            const=A(:,7)./A(:,1);
            A(:,[8 9])=A(:,[8 9])-const(:,[1 1]).*A(:,[2 3]);
            
            for i=1:nlhs
                z(:,3,i)=z(:,3,i)-const.*z(:,1,i);
            end
            
            %%%
            % Pivot again
            %%%
            [ju,midx2]=max(abs(A(:,[5 8])),[],2);
            pividx=[1 2 3 4 5 6 7 8 9
                1 2 3 7 8 9 4 5 6];
            
            idx=pividx(midx2,:)*nrows-repmat(((nrows-1):-1:0)',1,9);
            A=A(idx);
            
            pividx=[1 2 3
                1 3 2];
            sgn=[1;-1];
            sign=sign.*sgn(midx2);
            
            idx=pividx(midx2,:)*nrows-repmat(((nrows-1):-1:0)',1,3);
            
            %% pivioting moved to next for loop.
            %             for i=1:nlhs
            %                 hold=z(:,:,i);
            %                 z(:,:,i)=hold(idx);
            %             end
            
            if any(~A(:,5)),
                error('Some matrices are singular');
            end
            
            const=A(:,8)./A(:,5);
            A(:,9)=A(:,9)-const.*A(:,6);
            
            for i=1:nlhs
                hold=z(:,:,i);
                z(:,:,i)=hold(idx);
                z(:,3,i)=z(:,3,i)-const.*z(:,2,i);
            end
            
            %%
            % Back substitute
            %%
            
            for i=1:nlhs
                z(:,3,i)=z(:,3,i)./A(:,9);
                z(:,2,i)=(z(:,2,i)-z(:,3,i).*A(:,6))./A(:,5);
                z(:,1,i)=(z(:,1,i)-z(:,3,i).*A(:,3)-z(:,2,i).*A(:,2))./A(:,1);
            end
            if nargout>1,
                detA=sign.*A(:,1).*A(:,5).*A(:,9);
            end
        end
        function [z,detA,A]=Solve33Loop(A,z,dim)
            if nargout>1,
                detA=zeros(size(A,1),1);
            end
            for i=1:size(A,1),
                z(i,:)=([A(i,1:3);A(i,4:6);A(i,7:9)]\z(i,:)')';
                if nargout>1,
                    detA(i)=det([A(i,1:3);A(i,4:6);A(i,7:9)]);
                end
            end
        end
        function [invJ,detJ]=Inv33(J1,dim)
            if nargin==1,
                dim=3;
            end
            switch dim
                case 3
                    nz=size(J1,1);
                    z=zeros(nz,3,3);
                    on=ones(nz,1);
                    zr=zeros(nz,1);
                    z(:,:,1)=[on zr zr];
                    z(:,:,2)=[zr on zr];
                    z(:,:,3)=[zr zr on];
                    [z,detJ]=VMath.Solve33(J1,z);
                    invJ=[squeeze(z(:,:,1)) squeeze(z(:,:,2)) squeeze(z(:,:,3))];
                    invJ=VMath.Tran33(invJ);
                case 2
                    detJ=J1(:,1).*J1(:,4)-J1(:,2).*J1(:,3);
                    if any(detJ<1e-14),
                        error('A Jacobian is nearly singular')
                    end
                    invJ=[J1(:,4) -J1(:,2) -J1(:,3) J1(:,3)]./detJ(:,ones(1,4));
                case 1
                    invJ=1./J1(:,1);
                    detJ=J1(:,1);
                otherwise
                    error('Unknown dimension')
            end
        end
        function [detJ]=Det33(J1,dim)
            if nargin==1,
                dim=3;
            end
            switch dim
                case 3
                    detJ=J1(:,1).*(J1(:,5).*J1(:,9)-J1(:,6).*J1(:,8)) ...
                        -J1(:,2).*(J1(:,4).*J1(:,9)-J1(:,6).*J1(:,7)) ...
                        +J1(:,3).*(J1(:,4).*J1(:,8)-J1(:,5).*J1(:,7));
                case 2
                    detJ=J1(:,1).*J1(:,5)-J1(:,2).*J1(:,4);
                case 1
                    detJ=J1(:,1);
                otherwise
                    error('Unknown Dimension')
            end
        end
        function [A]=Tran33(A)
            %% Calculates the transpose of a 3x3 matrix
            A=A(:,[1 4 7 2 5 8 3 6 9]);
        end
        function [b]=Ax33(A,x)
            b=[A(:,1).*x(:,1)+A(:,2).*x(:,2)+A(:,3).*x(:,3) ...
                A(:,4).*x(:,1)+A(:,5).*x(:,2)+A(:,6).*x(:,3) ...
                A(:,7).*x(:,1)+A(:,8).*x(:,2)+A(:,9).*x(:,3)];
        end
        %%%% 2x2 [A11 A12 A21 A22]
        function [invJ,detJ]=Inv22(J1)
            detJ=J1(:,1).*J1(:,4)-J1(:,2).*J1(:,3);
            if any(detJ<1e-14),
                error('A Matrix is nearly singular')
            end
            invJ=[J1(:,4) -J1(:,2) -J1(:,3) J1(:,3)]./detJ(:,ones(1,4));
        end
        function [z,detA,A]=Solve22(A,z)
            %
            %  A is stored as A=(A11 A12 A21 A22 )
            %         
            %  z is a vector stored as z(:,:,i)=(z1 z2)
            %
            [n,m,nlhs]=size(z);
            detA=A(:,1).*A(:,4)-A(:,2).*A(:,3);
            A=[A(:,4) -A(:,2) -A(:,3) A(:,1)]./detA(:,[1 1 1 1]);
            for i=1:nlhs,
                z(:,:,i)=[A(:,1).*z(:,1,i)+A(:,2).*z(:,2,i) A(:,3).*z(:,1,i)+A(:,4).*z(:,2,i)];
            end
        end
        function [detJ]=Det22(J1)
            detJ=J1(:,1).*J1(:,4)-J1(:,2).*J1(:,3);
        end
        function [A]=Tran22(A)
            A=A(:,[1 3 2 4]);
        end
        function [b]=Ax22(A,x)
            b=[A(:,1).*x(:,1)+A(:,2).*x(:,2) A(:,3).*x(:,1)+A(:,4).*x(:,2)];
        end
        
        %%%%% symmetric 3x3 [A11 A22 A33 A23 A13 A12]
        function [z,detJ]=SSolve33(A,z)
            %
            %  A is symmetric and stored as A=(A11 A22 A33 A23 A31 A12)
            %  z is a vector stored as z=(z1 z2 z3)
            %
            %
            %
            % We'll use gauss elimination to solve
            %
            if nargout==2,
                [z,detJ]=VMath.Solve33(A(:,[1 6 5 6 2 4 5 4 3]),z);
            else
                [z]=VMath.Solve33(A(:,[1 6 5 6 2 4 5 4 3]),z);
            end
        end
        function D=RtDR(R,D)
            %
            %  calculates d=R'*D*R
            %
            %
            % R=(R11 R12 R13 R21 R22 R23 R31 R32 R33)
            % D=(D11 D22 D33 D23 D31 D12)
            % for R*S*R' => R=R(:,[1 4 7 2 5 8 3 6 9]);
            %
            A=zeros(size(R,1),9);
            
            A(:,1)=R(:,1).*D(:,1)+R(:,4).*D(:,4)+R(:,7).*D(:,6);
            A(:,2)=R(:,1).*D(:,4)+R(:,4).*D(:,2)+R(:,7).*D(:,5);
            A(:,3)=R(:,1).*D(:,6)+R(:,4).*D(:,5)+R(:,7).*D(:,3);
            %
            A(:,4)=R(:,2).*D(:,1)+R(:,5).*D(:,4)+R(:,8).*D(:,6);
            A(:,5)=R(:,2).*D(:,4)+R(:,5).*D(:,2)+R(:,8).*D(:,5);
            A(:,6)=R(:,2).*D(:,6)+R(:,5).*D(:,5)+R(:,8).*D(:,3);
            %
            A(:,7)=R(:,3).*D(:,1)+R(:,6).*D(:,4)+R(:,9).*D(:,6);
            A(:,8)=R(:,3).*D(:,4)+R(:,6).*D(:,2)+R(:,9).*D(:,5);
            A(:,9)=R(:,3).*D(:,6)+R(:,6).*D(:,5)+R(:,9).*D(:,3);
            %%
            
            D(:,1)=A(:,1).*R(:,1)+A(:,2).*R(:,4)+A(:,3).*R(:,7);
            D(:,2)=A(:,4).*R(:,2)+A(:,5).*R(:,5)+A(:,6).*R(:,8);
            D(:,3)=A(:,7).*R(:,3)+A(:,8).*R(:,6)+A(:,9).*R(:,9);
            %%
            D(:,4)=A(:,1).*R(:,2)+A(:,2).*R(:,5)+A(:,3).*R(:,8);
            D(:,5)=A(:,4).*R(:,3)+A(:,5).*R(:,6)+A(:,6).*R(:,9);
            D(:,6)=A(:,1).*R(:,3)+A(:,2).*R(:,6)+A(:,3).*R(:,9);
            %
            
        end
        function [invJ,detJ]=SInv33(J1,dim)
            %% A better way to do this would be to call Solve33 with z the
            %% identity matrix.
            if nargin==1,
                dim=3;
            end
            switch dim
                case 3
                    nz=size(J1,1);
                    z=zeros(nz,3,3);
                    on=ones(nz,1);
                    zr=zeros(nz,1);
                    z(:,:,1)=[on zr zr];
                    z(:,:,2)=[zr on zr];
                    z(:,:,3)=[zr zr on];
                    [z,detJ]=VMath.SSolve33(J1,z);
                    invJ=[squeeze(z(:,1,1)) squeeze(z(:,2,2)) squeeze(z(:,3,3)) ...
                        squeeze(z(:,2,3)) squeeze(z(:,1,3)) squeeze(z(:,1,2))];
                case 2
                    detJ=J1(:,1).*J1(:,2)-J1(:,3).*J1(:,3);
                    if any(detJ<1e-14),
                        warning('A matrix is nearly singular')
                    end
                    invJ=[J1(:,2) -J1(:,3) J1(:,1)]./detJ(:,ones(1,3));
                case 1
                    invJ=1./J1(:,1);
                    detJ=J1(:,1);
                otherwise
                    error('Unknown dimension')
            end
        end
        function b=SAx33(A,x)
            b=[A(:,1).*x(:,1)+A(:,6).*x(:,2)+A(:,5).*x(:,3) ...
                A(:,6).*x(:,1)+A(:,2).*x(:,2)+A(:,4).*x(:,3) ...
                A(:,5).*x(:,1)+A(:,4).*x(:,2)+A(:,3).*x(:,3)];
        end
        %%%%% symmetric 2x2 [A11 A22 A12]
        function [Cinv,detC]=SInv22(C)
            detC=C(:,2).*C(:,1)-C(:,3).*C(:,3);
            Cinv=[C(:,2) C(:,1) -C(:,3)]./detC(:,ones(3,1));
        end
        function b=SAx22(A,x)
            b=[A(:,1).*x(:,1)+A(:,3).*x(:,2) A(:,3).*x(:,1)+A(:,2).*x(:,2)];
        end
        function z=SSolve22(A,z)
            [ju,idx]=max(abs(A(:,[1 3])),[],2);
            idx1=find(idx==1);
            if ~isempty(idx1),
                z(idx1,2)=(z(idx1,2)-A(idx1,3).*z(idx1,1)./A(idx1,1))./(A(idx1,2)-A(idx1,3).*A(idx1,3)./A(idx1,1));
                z(idx1,1)=(z(idx1,1)-A(idx1,3).*z(idx1,2))./A(idx1,1);
            end
            idx2=find(idx==2);
            if ~isempty(idx2),
                z(idx2,:)=z(idx2,[2 1]);
                z(idx2,2)=(z(idx2,2)-A(idx2,1).*z(idx2,1)./A(idx2,3))./(A(idx2,3)-A(idx2,1).*A(idx2,2)./A(idx2,3));
                z(idx2,1)=(z(idx2,1)-A(idx2,2).*z(idx2,2))./A(idx2,3);
            end
        end
        %%%%% Overdetermined 3x2 [A11 A21 A31 A21 A22 A23]
        function z=PInv32(A,z)
            ATA=[sum(A(:,1:3).*A(:,1:3),2) sum(A(:,4:6).*A(:,4:6),2) sum(A(:,1:3).*A(:,4:6),2)];
            z=[sum(A(:,1:3).*z,2) sum(A(:,4:6).*z,2)];
            z=VMath.SSolve22(ATA,z);
        end
        %%%%% Overdetermined 3x1 [A11 A21 A31]
        function z=PInv31(A,z)
            ATA=sum(A(:,1:3).*A(:,1:3),2);
            z=sum(A(:,1:3).*z,2);
            z=z./ATA;
        end
        %%%%% Overdetermined 2x1 [A11 A21];
        function z=PInv21(A,z)
            ATA=sum(A(:,1:2).*A(:,1:2),2);
            z=sum(A(:,1:2).*z,2);
            z=z./ATA;
        end
        %%%%% Vector tools
        function [c]=Dot(a,b)
            c=sum(a.*b,2);
        end
        function [c]=Cross(a,b)
            c=[a(:,2).*b(:,3) - a(:,3).*b(:,2) ...
                a(:,3).*b(:,1) - a(:,1).*b(:,3) ...
                a(:,1).*b(:,2) - a(:,2).*b(:,1)];
        end
        function [c]=Norm(a)
            c=sqrt(sum(a.*a,2));
        end
        %%%%% Debugging tools
        function A=VMath2Matrix(A,i)
            if nargin==1,
                i=1;
            end
            switch size(A,2),
                case 9
                    A=[A(i,1:3);A(i,4:6);A(i,7:9)];
                case 6
                    A=[A(i,[1 6 5]);A(i,[6 2 4]);A(i,[5 4 3])];
                case 4
                    A=[A(i,1:2);A(i,3:4)];
                case 3
                    A=[A(i,[1 3]);A(i,[3 2])];
                otherwise
                    error('Unknown dimension')
            end
        end
    end
end