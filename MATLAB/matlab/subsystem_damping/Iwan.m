classdef Iwan < Joint
    %
    % obj=Iwan(params);
    %
    %  params=[ chi, phi_max,R, S]
    %  y0 = slider positions (nsliders+1) (+1 for the macroslip)
    %  k = current value of tangent stiffness
    %  f = force in slider
    
    % Written by Todd Simmermacher 3/2009
    %

    properties (Access=public)
    end
    
    properties (Access=private)
        SliderDisplacements % current slider displacements
        %% Iwan parameters
        Chi     % slope
        PhiMax % Break-free force of slider
        R       % Iwan parameter
        S       % Iwan parameter
        %%
        NSliders = 51% number of sliders
        Alpha=1.2 % Internal variable that describes the Geometric distribution of the sliders
        %% Constants
        Phi
        Phi1
        Phi2
        %%
        % Variables to hold during iterations
        SDhold  % hold slider displacements
        %Fhold   % hold force
    end
    
    methods
        function obj=Iwan(params)
            if nargin>0
                %%params=[ chi, phi_max,R, S]
                obj.Chi=params(1);
                obj.PhiMax=params(2);
                obj.R=params(3);
                obj.S=params(4);
                
                obj.SliderDisplacements =zeros(obj.NSliders,1);
                obj.SDhold=zeros(obj.NSliders,1);
                % Set up constants
                dphi=params(2)*(obj.Alpha-1)/(obj.Alpha^(obj.NSliders-1) -1);
                ju=cumsum([0 dphi*1.2.^(0:obj.NSliders-2)]);
                phi=ju(1:end-1)+diff(ju)/2;
                
                obj.Phi=[phi(:);obj.PhiMax];
                
                obj.Phi1=ju(:).^(obj.Chi+1);
                obj.Phi2=ju(:).^(obj.Chi+2);
                
                % Initialize stiffness
                obj=obj.updateStiffness(0,0);
                obj.K0=obj.Kt;
            end
        end
        function [obj,force]=updateForce(obj,x,v)
            %
            % Differential Displacement (x)
            % Differential Velocity (v) is unused in Iwan
            %
            
            % define some commonly used quantities
            dif=x-obj.SliderDisplacements;
            dif=dif(:);
            
            % figure out what sliders are sliding
            Il=find((abs(dif(1:end-1))-obj.Phi(1:end-1))<0); %fabs(u-x)<phi
            Ig=find((abs(dif(1:end-1))-obj.Phi(1:end-1))>=0); %fabs(u-x)>=phi
            
            % some more useful quantities
            dif(~dif)=1e-30;
            dir=sign(dif);
            
            %
            % calculate the force in the element
            f1=obj.R*sum(dif(Il).*(obj.Phi1(Il+1)-obj.Phi1(Il)))/(obj.Chi+1); %not sliding contribution
            f2=obj.R*sum(dir(Ig).*(obj.Phi2(Ig+1)-obj.Phi2(Ig)))/(obj.Chi+2); %sliding contribution
            obj.Force=f1+f2;
            
            % add in the macroslip singularity
            if abs(dif(end))<obj.Phi(end),
                obj.Force=obj.Force+obj.S*dif(end);
            else
                obj.Force=obj.Force+dir(end)*obj.S*obj.PhiMax;
            end
            
            % update states
            obj.SDhold=obj.SliderDisplacements;
            obj.SDhold(Ig)=x-dir(Ig).*obj.Phi(Ig);
            if abs(dif(end)) > obj.Phi(end),
                obj.SDhold(end)=x-dir(end)*obj.Phi(end);
            end
            
            %% just so we can return the force
            force=obj.Force;
            
        end
        function [obj,kt]=updateStiffness(obj,x,v)
            kt=obj.calcStiffness(x,v);
            obj.Kt=kt;
        end
        function [kt]=calcStiffness(obj,x,v)
            %
            % Differential Displacement (x)
            % Differential Velocity (v) is unused in Iwan
            %
            dif=x-obj.SliderDisplacements;
            dif=dif(:);
            Il=find((abs(dif(1:end-1))-obj.Phi(1:end-1))<0); %fabs(u-x)<phi
            
            if ~isempty(Il),
                f1=obj.R*sum((obj.Phi1(Il+1) - obj.Phi1(Il)))/(obj.Chi+1);
                if abs(dif(end))<obj.PhiMax,
                    kt=f1+obj.S;
                else
                    kt=f1;
                end
            end
        end
        function obj=updateStates(obj)
            obj.SliderDisplacements=obj.SDhold;
        end
    end
end
